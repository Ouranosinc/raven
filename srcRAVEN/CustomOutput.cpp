/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team, Ayman Khedr
  ----------------------------------------------------------------*/

#include "CustomOutput.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the CCustomOutput constructor
/// \param variable [in] Output variable assessed in custom output file
/// \param sv [in] State variable output type (if output is a SV)
/// \param sv_index [in] State variable index (if output is a state variable)
/// \param force_string [in] Forcing function name (if output is a forcing function)
/// \param stat [in] Time aggregate statistic (average, maximum, range, etc.)
/// \param time_aggregation [in] How frequently data is aggregated (monthly, daily, hourly, etc.)
/// \param space_aggregation [in] How data is spatially aggregated (by HRU, by basin, etc.)
/// \param *pMod [in] Pointer to model
/// \param &Options [in] Global model options information
//
CCustomOutput::CCustomOutput( const diagnostic    variable,
                              const sv_type       sv,
                              const int           sv_index,
                              const int           sv_index2,
                              const string        force_string,
                              const agg_stat      stat,
                              const time_agg      time_aggregation,
                              const spatial_agg   space_aggregation,
                              const string        filename_spec,
                              const int           kk,
                              const CModel       *pMod,
                              const optStruct    &Options)
{
// no point in going further if we don't have a model
  ExitGracefullyIf(pMod==NULL,"CCustomOutput Constructor: NULL model",BAD_DATA);

// set up the objects member variables
  _netcdf_ID    = -9; //doesn't exist
  _var      =variable;                                      // forcing variable, state variable, flux, etc.
  _svtype   =sv;                                            // state variable type (if output var is a SV)
  _svind    =sv_index;                                      // state variable index (if output var is a SV or flux)
  _svind2   =sv_index2;                                     // state variable target index (if output var is a flux between)
  _force_str=force_string;                                  // name of forcing variable (if output var is a forcing function)

  _spaceAgg  =space_aggregation;                            // how we're aggregation the data (ByHRU,BySubbasin, etc.)
  _timeAgg   =time_aggregation;                             // yearly, monthly, etc.
  _aggstat   =stat;                                         // the statistic we're calculating

  pModel     =pMod;                                         // pointer to the main model


  _varName        ="UNKNOWN";
  _varUnits       ="none";
  _timeAggStr     ="UNKNOWN";
  _statStr        ="UNKNOWN";
  _spaceAggStr    ="UNKNOWN";

  num_data      =1;
  data          =NULL;  // pointer to data storage for aggregation

  _hist_min     =0;
  _hist_max     =10;
  _nBins        =1; //Arbitrary default

  count         =0;

  kk_only       =kk;

  ExitGracefullyIf((kk_only==DOESNT_EXIST) && (_spaceAgg==BY_SELECT_HRUS),
                   "CCustomOutput Constructor: invalid HRU group index for Select HRU Aggregation. Undefined HRU group?",BAD_DATA);
  _time_index=0;

  //-------------------------------------------------------------
  // figure out filename here
  // RunName if available
  ostrstream FILENAME;
  if (Options.run_name!=""){FILENAME<<Options.output_dir<<Options.run_name<<"_";}
  else                     {FILENAME<<Options.output_dir;                       }
      
  // Get the variable name (either the state variable name of the forcing variable name)
  // and the units
  switch(_var)
  {
  case (VAR_STATE_VAR):
  {
    sv_type typ=pModel->GetStateVarType(_svind);
    int     ind=pModel->GetStateVarLayer(_svind);
    _varName  = CStateVariable::SVTypeToString(typ,ind);
    _varUnits = CStateVariable::GetStateVarUnits(typ);
    if (pModel->GetStateVarType(_svind)==CONSTITUENT){ _varUnits = "mg/L"; } //overridden for concentrations
    break;
  }
  case (VAR_FORCING_FUNCTION):
  {
    _varName = _force_str;
    forcing_type typ=GetForcingTypeFromString(_force_str);
    ExitGracefullyIf(typ==F_UNRECOGNIZED,
                     "CCustomOutput Constructor: invalid forcing type string",BAD_DATA);
    _varUnits = GetForcingTypeUnits(typ);
    break;
  }
  case (VAR_TO_FLUX) : //cumulative flux
  {
    sv_type typ=pModel->GetStateVarType(_svind);
    int     ind=pModel->GetStateVarLayer(_svind);
    _varName  = "TO_"+CStateVariable::SVTypeToString(typ,ind);
    _varUnits = CStateVariable::GetStateVarUnits(typ);
    break;
  }
  case (VAR_FROM_FLUX) : //cumulative flux
  {
    sv_type typ=pModel->GetStateVarType(_svind);
    int     ind=pModel->GetStateVarLayer(_svind);
    _varName  = "FROM_"+CStateVariable::SVTypeToString(typ,ind);
    _varUnits = CStateVariable::GetStateVarUnits(typ);
    break;
  }
  case (VAR_BETWEEN_FLUX) :
  {
    sv_type typ=pModel->GetStateVarType(_svind);
    int     ind=pModel->GetStateVarLayer(_svind);
    ExitGracefullyIf(_svind2==DOESNT_EXIST,"Invalid second SV index in :CustomOutput :Between command",BAD_DATA_WARN);
    sv_type typ2=pModel->GetStateVarType(_svind2);
    int     ind2=pModel->GetStateVarLayer(_svind2);
    _varName  = "BETWEEN_"+CStateVariable::SVTypeToString(typ,ind)+"_AND_"+CStateVariable::SVTypeToString(typ2,ind2);
    _varUnits = CStateVariable::GetStateVarUnits(typ)+"/d";
    break;
  }
  default:
  {
    _varName = "UNKNOWN";
    _varUnits = "none";
    break;
  }
  }
  FILENAME<<_varName<<"_";

// temporal aggregation
  switch(_timeAgg)
  {
  case YEARLY:                    _timeAggStr="Yearly"; break;
  case MONTHLY:                   _timeAggStr="Monthly"; break;
  case DAILY:                     _timeAggStr="Daily"; break;
  case WATER_YEARLY:              _timeAggStr="WYearly";break;
  case EVERY_TSTEP:               _timeAggStr="Continuous"; break;
  }
  FILENAME<<_timeAggStr<<"_";

// statistic
  switch(_aggstat)
  {
  case AGG_AVERAGE:               _statStr="Average"; break;
  case AGG_MAXIMUM:               _statStr="Maximum"; break;
  case AGG_MINIMUM:               _statStr="Minimum"; break;
  case AGG_RANGE:                 _statStr="Range"; break;
  case AGG_MEDIAN:                _statStr="Median"; break;
  case AGG_95CI:                  _statStr="95%"; break;
  case AGG_QUARTILES:             _statStr="Quartiles"; break;
  case AGG_HISTOGRAM:             _statStr="Histogram"; break;
  }
  FILENAME<<_statStr<<"_";

// spatial aggregation
  switch(_spaceAgg)
  {
  case BY_HRU:                    _spaceAggStr="ByHRU"; break;
  case BY_BASIN:                  _spaceAggStr="BySubbasin"; break;
  case BY_WSHED:                  _spaceAggStr="ByWatershed"; break;
  case BY_HRU_GROUP:              _spaceAggStr="ByHRUGroup"; break;
  case BY_SELECT_HRUS:            _spaceAggStr="ByHRU"; break;
  }
  FILENAME<<_spaceAggStr;

  // filename extension
  switch(Options.output_format)
  {
  case OUTPUT_ENSIM:      FILENAME<<".tb0"<<ends; break;
  case OUTPUT_NETCDF:     FILENAME<<".nc" <<ends; break;
  case OUTPUT_STANDARD:
  default:                FILENAME<<".csv"<<ends; break;
  }

  if (filename_spec==""){_filename=FILENAME.str();}
  else                  {
    if (Options.run_name!=""){_filename=Options.output_dir+Options.run_name+"_"+filename_spec;}
    else                     {_filename=Options.output_dir+filename_spec; } // \todo [QA/QC]: should check for proper extension of filename_spec
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CCustomOutput::~CCustomOutput()
{
  delete [] data; data=NULL;
  CloseFiles();
}

//////////////////////////////////////////////////////////////////
/// \brief Set histogram-related parameters
/// \param minv [in] Histogram minimum
/// \param maxv [in] Histogram maximum
/// \param numBins [in] Histogram number of bins
//
void CCustomOutput::SetHistogramParams(const double minv,const double maxv, const int numBins)
{
  _hist_min=minv;
  _hist_max=maxv;
  _nBins=numBins;
}

///////////////////////////////////////////////////////////////////
/// \brief Allocates memory and initialize data storage of a CCustomOutput object
/// \remarks Called prior to simulation. Determines size of and allocates memory for (member) data[][] array needed in statistical calculations
/// \param &Options [in] Global model options information
//
void CCustomOutput::InitializeCustomOutput(const optStruct &Options)
{
// figure out how much space we need
  if      (_spaceAgg==BY_HRU        ){num_data=pModel->GetNumHRUs();}
  else if (_spaceAgg==BY_BASIN      ){num_data=pModel->GetNumSubBasins();}
  else if (_spaceAgg==BY_WSHED      ){num_data=1;}
  else if (_spaceAgg==BY_HRU_GROUP  ){num_data=pModel->GetNumHRUGroups();}
  else if (_spaceAgg==BY_SELECT_HRUS){num_data=pModel->GetHRUGroup(kk_only)->GetNumHRUs();}

  if      (_aggstat==AGG_AVERAGE){num_store=1;}
  else if (_aggstat==AGG_MAXIMUM){num_store=1;}
  else if (_aggstat==AGG_MINIMUM){num_store=1;}
  else if (_aggstat==AGG_RANGE)  {num_store=2;}
  else if ((_aggstat==AGG_MEDIAN   ) ||
           (_aggstat==AGG_95CI     ) ||
           (_aggstat==AGG_QUARTILES) ||
           (_aggstat==AGG_HISTOGRAM))
  {
    if      (_timeAgg==YEARLY      ){num_store=(int)ceil(366/Options.timestep)+1;}
    else if (_timeAgg==WATER_YEARLY){num_store=(int)ceil(366/Options.timestep)+1;}
    else if (_timeAgg==MONTHLY     ){num_store=(int)ceil( 31/Options.timestep)+1;}
    else if (_timeAgg==DAILY       ){num_store=(int)ceil(  1/Options.timestep)+1;}
    else if (_timeAgg==EVERY_TSTEP ){num_store=1;}
  }
  else{
    ExitGracefully("CCustomOutput::InitializeCustomOutput(): bad aggregator",BAD_DATA);
  }

  // allocate space for the aggregation
  data =new double *[num_data];
  for (int k=0;k<num_data;k++)
  {
    data[k]=NULL;
    data[k]=new double [num_store];
    ExitGracefullyIf(data[k]==NULL,"CCustomOutput constructor",OUT_OF_MEMORY);
    for (int a=0;a<num_store;a++){data[k][a]=0.0;}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Open a stream to the file and write header info
//
void CCustomOutput::WriteFileHeader(const optStruct &Options)
{
  _CUSTOM.open(_filename.c_str());
  if (_CUSTOM.fail()){
    WriteWarning("CCustomOutput::WriteFileHeader: Unable to create file "+_filename,true);
  }

  // write the appropriately formatted header
  switch(Options.output_format)
  {
  case OUTPUT_STANDARD:
  default:
    WriteCSVFileHeader();
    return;
    break;
  case OUTPUT_ENSIM:
    WriteEnSimFileHeader(Options);
    return;
    break;
  case OUTPUT_NETCDF:
    WriteNetCDFFileHeader(Options);
    return;
    break;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write header info to "STANDARD" CSV file
//
void CCustomOutput::WriteCSVFileHeader(void)
{
  // standard csv output
  //-Line 1-
  if      (_timeAgg==YEARLY      ){_CUSTOM<<",";}
  else if (_timeAgg==MONTHLY     ){_CUSTOM<<",";}
  else if (_timeAgg==DAILY       ){_CUSTOM<<",";}
  else if (_timeAgg==EVERY_TSTEP ){_CUSTOM<<",,";}
  else if (_timeAgg==WATER_YEARLY){_CUSTOM<<",";}

  if      (_spaceAgg==BY_HRU        ){_CUSTOM<<"HRU:,";}
  else if (_spaceAgg==BY_BASIN      ){_CUSTOM<<"SubBasin:,";}
  else if (_spaceAgg==BY_WSHED      ){_CUSTOM<<"Watershed:,";}
  else if (_spaceAgg==BY_HRU_GROUP  ){_CUSTOM<<"HRUGroup:,";}
  else if (_spaceAgg==BY_SELECT_HRUS){_CUSTOM<<"HRU:,";}

  for (int k=0;k<num_data;k++)
  {
    string title;
    ostrstream  TMP;
    if      (_spaceAgg==BY_HRU        ){TMP<<pModel->GetHydroUnit(k)->GetID() <<ends;}
    else if (_spaceAgg==BY_HRU_GROUP  ){TMP<<pModel->GetHRUGroup(k)->GetName()<<ends;}
    else if (_spaceAgg==BY_WSHED      ){TMP<<"Watershed"                      <<ends;}
    else if (_spaceAgg==BY_BASIN      ){TMP<<pModel->GetSubBasin(k)->GetID()  <<ends;}
    else if (_spaceAgg==BY_SELECT_HRUS){TMP<<pModel->GetHRUGroup(kk_only)->GetHRU(k)->GetID()<<ends;}

    title=TMP.str();

    if      (_aggstat == AGG_MEDIAN   ){_CUSTOM<<title<<",";}
    else if (_aggstat == AGG_95CI     ){_CUSTOM<<title<<",,";}
    else if (_aggstat == AGG_QUARTILES){_CUSTOM<<title<<",,,";}
    else if (_aggstat == AGG_HISTOGRAM){
      _CUSTOM<<title;
      for (int bin=0;bin<_nBins;bin++){_CUSTOM<<",";}
    }
    else    {
      _CUSTOM<<title<<",";for (int i=1;i<num_store;i++){_CUSTOM<<",";}
    }
  }
  _CUSTOM<<endl;

  //-Line 2-
  if      (_timeAgg==YEARLY      ){_CUSTOM<<"time,year,";}
  else if (_timeAgg==WATER_YEARLY){_CUSTOM<<"time,water year,";}
  else if (_timeAgg==MONTHLY     ){_CUSTOM<<"time,month,";}
  else if (_timeAgg==DAILY       ){_CUSTOM<<"time,day,";}
  else if (_timeAgg==EVERY_TSTEP ){_CUSTOM<<"time,day,hour,";}
  for (int k=0;k<num_data;k++)
  {
    if      (_aggstat==AGG_AVERAGE  ){_CUSTOM<<"mean,";}
    else if (_aggstat==AGG_MAXIMUM  ){_CUSTOM<<"maximum,";}
    else if (_aggstat==AGG_MINIMUM  ){_CUSTOM<<"minimum,";}
    else if (_aggstat==AGG_RANGE    ){_CUSTOM<<"minimum,maximum,";}
    else if (_aggstat==AGG_MEDIAN   ){_CUSTOM<<"median,";}
    else if (_aggstat==AGG_95CI     ){_CUSTOM<<"5% quantile,95% quantile,";}
    else if (_aggstat==AGG_QUARTILES){_CUSTOM<<"25% quartile,50% quartile,75% quartile,";}
    else if (_aggstat==AGG_HISTOGRAM){
      for (int bin=0;bin<_nBins;bin++){
        _CUSTOM<<_hist_min+bin*(_hist_max-_hist_min)/_nBins<<"-"<<_hist_min+(bin+1)*(_hist_max-_hist_min)/_nBins<<",";
      }
    }
  }
  _CUSTOM<<endl;
}

//////////////////////////////////////////////////////////////////
/// \brief Write header info to EnSim file
//
void CCustomOutput::WriteEnSimFileHeader(const optStruct &Options)
{
// Write the header
  _CUSTOM<<"#########################################################################"<<endl;
  _CUSTOM<<":FileType tb0 ASCII EnSim 1.0"<<endl;
  _CUSTOM<<"#"<<endl;
  _CUSTOM<<":Application   Raven"<<endl;
  if(!Options.benchmarking){
    _CUSTOM<<":Version       "<<Options.version <<endl;
    _CUSTOM<<":CreationDate  "<<GetCurrentTime()<<endl;
  }
  _CUSTOM<<"#"<<endl;
  _CUSTOM<<"#------------------------------------------------------------------------"<<endl;
  _CUSTOM<<"#"<<endl;

// some meta data information describing the contents of this file
  _CUSTOM<<":RunName              "<<Options.run_name <<endl;
  _CUSTOM<<"#"<<endl;
  if      ((_var == VAR_FORCING_FUNCTION) && (_timeAgg==EVERY_TSTEP)){ _CUSTOM<<":Format PeriodEnding "<<endl; }//period ending
  else if ((_var == VAR_FROM_FLUX       ) && (_timeAgg==EVERY_TSTEP)){ _CUSTOM<<":Format Instantaneous"<<endl;  }//snapshot
  else if ((_var == VAR_TO_FLUX         ) && (_timeAgg==EVERY_TSTEP)){ _CUSTOM<<":Format Instantaneous"<<endl;  }//snapshot
  else if ((_var == VAR_BETWEEN_FLUX    ) && (_timeAgg==EVERY_TSTEP)){ _CUSTOM<<":Format PeriodEnding"<<endl;  }//period ending
  else if ((_var == VAR_STATE_VAR       ) && (_timeAgg==EVERY_TSTEP)){ _CUSTOM<<":Format Instantaneous"<<endl;  }//snapshot
  else                                                               { _CUSTOM<<":Format PeriodStarting"<<endl;  }//period starting
  _CUSTOM<<"#"<<endl;
  _CUSTOM<<"# Custom output meta-data"<<endl;
  _CUSTOM<<"#"<<endl;
  _CUSTOM<<":Variable             "<<_varName    <<endl;
  _CUSTOM<<":TemporalAggregation  "<<_timeAggStr <<endl;
  _CUSTOM<<":Statistic            "<<_statStr    <<endl;
  _CUSTOM<<":SpatialAggregation   "<<_spaceAggStr<<endl;
  _CUSTOM<<"#"<<endl;

// Now define the column meta data
  _CUSTOM<<":ColumnMetaData"<<endl;

// First column is a year, month or date(the default)
// second column is model time (in days)
  string col1Name="Date", col1Units="date", col1Type="date";
  string colType="float"; // the default, (histogram however is integer)
  switch(_timeAgg)
  {
  case YEARLY: col1Name="Year"; col1Units="years"; col1Type="int"; break;
  case WATER_YEARLY: col1Name="WaterYear"; col1Units="years"; col1Type="int"; break;
  case MONTHLY: col1Name="MonthEnd"; col1Units="months"; col1Type="date"; break;
  default: break;
  }

  _CUSTOM<<"  :ColumnName "<<col1Name<<" time ";

// Build the rest of the column names
  int dataColumnCount = 0;
  for (int k=0;k<num_data;k++)
  {
    ostrstream  curSpaceIdentifier;
    switch(_spaceAgg)
    {
    case BY_HRU:                    curSpaceIdentifier<<"HRU_"      <<pModel->GetHydroUnit(k)->GetID() <<ends; break;
    case BY_HRU_GROUP:              curSpaceIdentifier<<"HRUGroup_" <<pModel->GetHRUGroup(k)->GetName()<<ends; break;
    case BY_WSHED:                  curSpaceIdentifier<<"Watershed_"<<"Watershed"                      <<ends; break;
    case BY_BASIN:                  curSpaceIdentifier<<"SubBasin_" <<pModel->GetSubBasin(k)->GetID()  <<ends; break;
    case BY_SELECT_HRUS:            curSpaceIdentifier<<"HRU_"      <<pModel->GetHRUGroup(kk_only)->GetHRU(k)->GetID()  <<ends; break;
    }
    string title=curSpaceIdentifier.str();

    switch(_aggstat)
    {
    case AGG_AVERAGE:               _CUSTOM<<title<<"_mean "; dataColumnCount++; break;
    case AGG_MAXIMUM:               _CUSTOM<<title<<"_max "; dataColumnCount++; break;
    case AGG_MINIMUM:               _CUSTOM<<title<<"_min "; dataColumnCount++; break;
    case AGG_RANGE:                 _CUSTOM<<title<<"_min "<<title<<"_max "; dataColumnCount+=2; break;
    case AGG_MEDIAN:                _CUSTOM<<title<<"_median "; dataColumnCount++; break;
    case AGG_95CI:                  _CUSTOM<<title<<"_5%_quantile "<<title<<"_95%_quantile "; dataColumnCount+=2; break;
    case AGG_QUARTILES:             _CUSTOM<<title<<"_25%_quartile "<<title<<"_50%_quartile "<<title<<"_75%_quartile "; dataColumnCount+=3; break;
    case AGG_HISTOGRAM:
    {
      for (int bin=0;bin<_nBins;bin++)
      {_CUSTOM<<title<<"_"<<_hist_min+bin*(_hist_max-_hist_min)/_nBins << "-" << _hist_min+(bin+1)*(_hist_max-_hist_min)/_nBins<<" ";}
      dataColumnCount+=_nBins;
      _varUnits="count";
      colType="int";
    }
    break;
    }
  }

  _CUSTOM<<endl;

  int colFormat;
  if      ((_var == VAR_FORCING_FUNCTION) && (_timeAgg==EVERY_TSTEP)){ colFormat = -1; }//period ending
  else if ((_var == VAR_FROM_FLUX       ) && (_timeAgg==EVERY_TSTEP)){ colFormat = 0;  }//snapshot (cumulative)
  else if ((_var == VAR_TO_FLUX         ) && (_timeAgg==EVERY_TSTEP)){ colFormat = 0;  }//snapshot (cumulative)
  else if ((_var == VAR_BETWEEN_FLUX    ) && (_timeAgg==EVERY_TSTEP)){ colFormat = -1; }//period ending
  else if ((_var == VAR_STATE_VAR       ) && (_timeAgg==EVERY_TSTEP)){ colFormat = 0;  }//snapshot
  else                                                               { colFormat = 1;  }//period starting

  _CUSTOM<<"  :ColumnUnits "<<col1Units<<" days";
  for(int i=0;i<dataColumnCount;i++) {_CUSTOM<<" "<<_varUnits;}
  _CUSTOM<<endl;

  _CUSTOM<<"  :ColumnType "<<col1Type<<" float";
  for(int i=0;i<dataColumnCount;i++) {_CUSTOM<<" "<<colType;}
  _CUSTOM<<endl;

  _CUSTOM << "  :ColumnFormat 0 0";
  for(int i=0;i<dataColumnCount;i++) {_CUSTOM<<" "<<colFormat;}
  _CUSTOM<<endl;

  _CUSTOM<<":EndColumnMetaData"<<endl;

  _CUSTOM<<"#"<<endl;
  _CUSTOM<<":EndHeader"<<endl;
}
//////////////////////////////////////////////////////////////////
/// \brief Write header info to netCDF file
//
void CCustomOutput::WriteNetCDFFileHeader(const optStruct &Options)
{

#ifdef _RVNETCDF_

  time_struct tt;                                    // start time structure
  const int   ndims1 = 1;
  const int   ndims2 = 2;
  int         dimids1[ndims1];                       // array which will contain all dimension ids for a variable
  int         dimids2[ndims2];                       // array which will contain all dimension ids for a variable
  int         time_dimid, varid_time;                // dimension ID (holds number of time steps) and variable ID (holds time values) for time
  int         ndata_dimid,varid_data,varid_grps;     // dimension ID, variable ID for simulated data and group (HRU/SB) info

  int         retval;                                // error value for NetCDF routines
  size_t      start[1], count[1];                    // determines where and how much will be written to NetCDF
  string      tmp,tmp2;

  bool cant_support=(_aggstat==AGG_RANGE || _aggstat==AGG_95CI || _aggstat==AGG_QUARTILES || _aggstat==AGG_HISTOGRAM);
  if(cant_support){
    WriteWarning("CCustomOutput::WriteNetCDFFileHeader: cannot support custom aggregation by range/quartile/histogram/95CI with netCDF output",Options.noisy);
  }

  // create file and obtain _netcdf_ID 
  _CUSTOM.close(); //because netcdf does the opening/closing  
  retval = nc_create(_filename.c_str(), NC_CLOBBER|NC_NETCDF4, &_netcdf_ID);       HandleNetCDFErrors(retval);

  // ---------------------------------------------------------- 
  // time                                                       
  // ---------------------------------------------------------- 
  // (a) Define the DIMENSIONS. NetCDF will hand back an ID for each. 
  retval = nc_def_dim(_netcdf_ID, "time", NC_UNLIMITED, &time_dimid);              HandleNetCDFErrors(retval);

  // (b) Define the time variable. 
  dimids1[0] = time_dimid;
  retval = nc_def_var(_netcdf_ID, "time", NC_DOUBLE, ndims1,dimids1, &varid_time); HandleNetCDFErrors(retval);

  // (c) Assign units attributes to the netCDF VARIABLES. 
  //     --> converts start day into "days since YYYY-MM-DD HH:MM:SS" 
  char  starttime[200]; // start time string in format 'days since YYY-MM-DD HH:MM:SS'
  JulianConvert( 0.0,Options.julian_start_day, Options.julian_start_year, tt);
  strcpy(starttime, "days since ") ;
  strcat(starttime, tt.date_string.c_str()) ;
  strcat(starttime, " 00:00:00");
  retval = nc_put_att_text(_netcdf_ID, varid_time, "units"   , strlen(starttime)  , starttime);   HandleNetCDFErrors(retval);
  retval = nc_put_att_text(_netcdf_ID, varid_time, "calendar", strlen("gregorian"), "gregorian"); HandleNetCDFErrors(retval);

  // ---------------------------------------------------------- 
  // custom data                                      
  // ---------------------------------------------------------- 
  if (num_data > 0)
  {
    // (a) create dimension "ndata"                          
    retval = nc_def_dim(_netcdf_ID, "ndata", num_data, &ndata_dimid);                             HandleNetCDFErrors(retval);

    // (b) create variable "HRUID" or "SBID" or...
    string group_name="HRU_ID";
    string long_name="ID of HRU";
    if      (_spaceAgg==BY_HRU        ){group_name="HRU_ID";    long_name="ID of HRU";}
    else if (_spaceAgg==BY_BASIN      ){group_name="SBID";      long_name="ID of subbasin";}
    else if (_spaceAgg==BY_WSHED      ){group_name="Watershed"; long_name="entire watershed";}
    else if (_spaceAgg==BY_HRU_GROUP  ){group_name="HRUGroup";  long_name="HRU group name";}
    else if (_spaceAgg==BY_SELECT_HRUS){group_name="HRU_ID";    long_name="ID of HRU";}

    dimids1[0] = ndata_dimid;
    retval = nc_def_var(_netcdf_ID, group_name.c_str(), NC_STRING, ndims1, dimids1, &varid_grps); HandleNetCDFErrors(retval);

    //(c) set some attributes to variable "HRUID" or "SBID"  
    tmp =long_name;
    tmp2="timeseries_ID";
    retval = nc_put_att_text(_netcdf_ID, varid_grps, "long_name",  tmp.length(), tmp.c_str());    HandleNetCDFErrors(retval);
    retval = nc_put_att_text(_netcdf_ID, varid_grps, "cf_role"  , tmp2.length(),tmp2.c_str());    HandleNetCDFErrors(retval);
    
    //(d) create variable "custom_data" which will contain at the end custom data; shape=(tt,num_data) 
    string netCDFtag=_statStr+"_"+_varName;
    dimids2[0] = time_dimid;
    dimids2[1] = ndata_dimid;
    retval = nc_def_var(_netcdf_ID, netCDFtag.c_str(), NC_DOUBLE, ndims2, dimids2, &varid_data);    HandleNetCDFErrors(retval);    
    
    //(f) set some attributes to variable _netCDFtag 
    tmp=_timeAggStr+" "+_statStr+" "+_varName+" "+_spaceAggStr;
    tmp2=_varUnits;
    static double fill_val[] = {NETCDF_BLANK_VALUE}; 
    static double miss_val[] = {NETCDF_BLANK_VALUE}; 
    retval = nc_put_att_text  (_netcdf_ID, varid_data, "long_name"    ,  tmp.length(), tmp.c_str());HandleNetCDFErrors(retval);
    retval = nc_put_att_text  (_netcdf_ID, varid_data, "units"        , tmp2.length(),tmp2.c_str());HandleNetCDFErrors(retval);
    retval = nc_put_att_double(_netcdf_ID, varid_data, "_FillValue"   , NC_DOUBLE,1, fill_val);     HandleNetCDFErrors(retval);
    retval = nc_put_att_double(_netcdf_ID, varid_data, "missing_value", NC_DOUBLE,1, miss_val);     HandleNetCDFErrors(retval);

  }// end if (num_data>0)

  // End define mode. This tells netCDF we are done defining metadata. 
  retval = nc_enddef(_netcdf_ID);  HandleNetCDFErrors(retval);

  // write values to NetCDF 
  // write HRU/subbasin/HRU group names to variable "HRUID" or "SBID" or... 
  int k = 0;
  const char *group_name[1];                         // HRU/watershed/basin name
  for(k=0; k<num_data; k++)
  {
    ostrstream  TMP;
    string temp;
    if      (_spaceAgg==BY_HRU        ){TMP<<pModel->GetHydroUnit(k)->GetID() <<ends;}
    else if (_spaceAgg==BY_HRU_GROUP  ){TMP<<pModel->GetHRUGroup(k)->GetName()<<ends;}
    else if (_spaceAgg==BY_WSHED      ){TMP<<"Watershed"                      <<ends;}
    else if (_spaceAgg==BY_BASIN      ){TMP<<pModel->GetSubBasin(k)->GetID()  <<ends;}
    else if (_spaceAgg==BY_SELECT_HRUS){TMP<<pModel->GetHRUGroup(kk_only)->GetHRU(k)->GetID()<<ends;}
    temp=TMP.str();
    group_name[0]=temp.c_str();
    start[0]=k;
    count[0]=1;
    retval = nc_put_vara_string(_netcdf_ID,varid_grps,start,count,&group_name[0]);  HandleNetCDFErrors(retval);
  }

  if(Options.noisy){ cout<<"netCDF file header written for "<<_filename<<endl; }
#endif   // end compilation if NetCDF library is available
}


//////////////////////////////////////////////////////////////////
/// \brief Calculates quartiles of passed 1D data array
/// \note values[] array (size 'size') must be pre-sorted
/// \author code developed by Ayman Khedr, University of Waterloo 3A, 2011
/// \param *values [in] Array of values whose quartiles are to be found
/// \param size [in] Size of values[] array
/// \param &Q1 [out] First quartile value
/// \param &Q2 [out] Second quartile value
/// \param &Q3 [out] Third Quartile value
//
void GetQuartiles(const double *values, const int size, double &Q1, double &Q2, double &Q3)
{
  if (size%2==0)
  {
    if ((size/2)%2==0)
    { // since values is sorted index value quartile values associated with quartiles of index
      Q1=(values[(int) ceil((double)(size-1)*0.25)]+values[(int)floor((double)(size-1)*0.25)])/2;
    }
    else
    {
      Q1=values[(int)floor((double)(size-1)*0.25)];
    }
    Q2=(values[(int)ceil((double)(size-1)*0.5)]+values[(int)floor((double)(size-1)*0.5)])/2;
    if ((size/2)%2==0)
    {
      Q3=(values[(int)ceil((double)(size-1)*0.75)]+values[(int)floor((double)(size-2)*0.75)])/2;
    }
    else{
      Q3=values[(int)ceil((double)(size-1)*0.75)];
    }
  }
  else
  {
    if ((int)floor((double)(size-1)/2)%2 == 0)
    {
      Q1=values[(int)floor((double)(size-1)*0.25)];
    }
    else
    {
      Q1=(values[(int)ceil((double)(size-1)/4)]+values[(int)floor((double)(size-1)/4)])/2;
    }
    Q2=values[((size-1)/2)];
    if ((int)floor((double)(size-1)/2)%2 == 0)
    {
      Q3=(values[(int)ceil((double)(size-1)*0.75)]+values[(int)floor((double)(size-2)*0.75)])/2;
    }
    else {
      Q3=values[(int)ceil((double)(size-1)*0.75)];
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Write custom output to file by querying the model and calculating diagnostics
/// \brief now handles .csv, .nc, and .tb0 (Ensim) formats
/// \param &tt [in] Current model time
/// \param &Options [in] Global model options information
//
void CCustomOutput::WriteCustomOutput(const time_struct &tt,
                                      const optStruct   &Options)
{
  //each custom output gets values of vals data[nHRUs] (max size)
  double val=0;
  double dday;//,jul_day;
  int    dmon,dyr;
  string thisdate, thishour, yesterday;
  bool   reset;

  double t=tt.model_time;

  time_struct yest;
  JulianConvert(t-1,Options.julian_start_day,Options.julian_start_year,yest); //get previous day, yest
  yesterday=yest.date_string;
  thisdate=tt.date_string;
  dday    =tt.day_of_month;
  dmon    =tt.month;
  dyr     =tt.year;

  thishour =DecDaysToHours(tt.julian_day);

  if (t==0){return;} //initial conditions should not be printed to custom output, only period data.

  double *output=NULL;
  if(Options.output_format==OUTPUT_NETCDF){output=new double[num_data]; }

  //Check to see if it is time to write to file
  //------------------------------------------------------------------------------
  reset=false;
  if      ((_timeAgg==YEARLY)  && (dday==1) && (dmon==1))    {reset=true;}//Jan 1 - print preceding year
  else if ((_timeAgg==MONTHLY) && (dday==1))                 {reset=true;}//first day of month - print preceding month info
  else if ((_timeAgg==DAILY)   &&
           ((fabs(floor(t)-t) <0.5*Options.timestep) ||
            (fabs( ceil(t)-t) <0.5*Options.timestep)))       {reset=true;}//start of day (hopefully 0:00!)- print preceding day
  else if (_timeAgg==EVERY_TSTEP)                            {reset=true;}//every timestep
  else if ((_timeAgg==WATER_YEARLY) && (dday==1) && (dmon==Options.wateryr_mo))    {reset=true;}//Oct 1 - print preceding year
  //cout <<t <<" ->"<<fabs(floor(t)-t)<<"   "<<(fabs(floor(t)-t) <0.5*Options.timestep)<<endl;

  if (reset)
  {
    if(Options.output_format==OUTPUT_STANDARD)//=================================================================
    {
      if      (_timeAgg==YEARLY      ){_CUSTOM<<t<<","<<yest.year<<",";}
      else if (_timeAgg==MONTHLY     ){_CUSTOM<<t<<","<<yesterday.substr(0,7)<<",";}//trims day, e.g., just return 2011-02,
      else if (_timeAgg==DAILY       ){_CUSTOM<<t<<","<<yesterday<<",";}
      else if (_timeAgg==EVERY_TSTEP ){_CUSTOM<<t<<","<<thisdate<<","<<thishour<<","; }//period ending for forcing data
      else if (_timeAgg==WATER_YEARLY){_CUSTOM<<t<<","<<yest.year<<"-"<<yest.year+1<<",";}
    }
    else if(Options.output_format==OUTPUT_ENSIM)//===============================================================
    {
      switch(_timeAgg)
      {
      case YEARLY:       _CUSTOM<<yest.year<<" "<<t<<" "; break;
      case MONTHLY:      _CUSTOM<<yest.year<<"-"<<yest.month<<"-"<<yest.day_of_month<<" "<<t<<" "; break;
      case DAILY:        _CUSTOM<<yest.year<<"-"<<yest.month<<"-"<<yest.day_of_month<<" "<<t<<" "; break;
      case EVERY_TSTEP:
      {
        char dQuote = '"';
        _CUSTOM<<dQuote<<tt.date_string<<" "<<DecDaysToHours(tt.julian_day)<<dQuote<<" "<<t<<" "; //period ending for forcings
      }
      break;
      case WATER_YEARLY: _CUSTOM<<t<<","<<yest.year<<"-"<<yest.year+1<<","; break;
      }
    }
    else if(Options.output_format==OUTPUT_NETCDF)//=============================================================
    {
#ifdef _RVNETCDF_
      int retval,time_id;
      size_t start1[1], count1[1];  // determines where and how much will be written to NetCDF; 1D variable (time)
      double current_time[1];       // current time in days since start time
       
      if      (_timeAgg==YEARLY      ){current_time[0]=yest.model_time-yest.julian_day;}
      else if (_timeAgg==MONTHLY     ){current_time[0]=yest.model_time-yest.day_of_month;}
      else if (_timeAgg==DAILY       ){current_time[0]=yest.model_time;}
      else if (_timeAgg==EVERY_TSTEP ){current_time[0]=tt.model_time-Options.timestep; }
      else if (_timeAgg==WATER_YEARLY){
        double days_to_first_of_month=0;
        for (int m=0; m<Options.wateryr_mo;m++){days_to_first_of_month+=DAYS_PER_MONTH[m];}
        if (IsLeapYear(yest.year) && (Options.wateryr_mo>2)){days_to_first_of_month+=1;}
        current_time[0]=yest.model_time-yest.julian_day+days_to_first_of_month;
      }
      
      //start1[0] = int(round(current_time[0]/Options.timestep));   // element of NetCDF array that will be written
      start1[0]=_time_index;
      count1[0] = 1;                                              // writes exactly one time step

      retval = nc_inq_varid (_netcdf_ID, "time",   &time_id);                             HandleNetCDFErrors(retval);
      retval = nc_put_vara_double(_netcdf_ID, time_id, start1, count1, &current_time[0]); HandleNetCDFErrors(retval);
#endif
    }
  }

  bool is_concentration=false;

  is_concentration = (_var == VAR_STATE_VAR) && (pModel->GetStateVarType(_svind)==CONSTITUENT);

  //Sift through HRUs, BASINs or watershed, updating aggregate statistics
  //--------------------------------------------------------------------------
  //numdata=1 if BY_WATERSHED, =nSubBasins if BY_BASIN, =nHRUs if BY_HRU...
  for (int k=0;k<num_data;k++)
  {
    //---access current diagnostic variable (from end of timestep)------------
    if (is_concentration){
      int m = pModel->GetStateVarLayer(_svind);
      int i_stor=pModel->GetTransportModel()->GetWaterStorIndexFromLayer(m);
      double conv=(MM_PER_METER/LITER_PER_M3);
      val=-9999;
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetStateVarValue(_svind)/pModel->GetHydroUnit(k)->GetStateVarValue(i_stor)*conv;}
      else {
        ExitGracefully("CustomOutput: cannot currently generate basin,watershed, or hru group based aggregate constituent concentrations",STUB);
      }
      /*else if (_spaceAgg==BY_BASIN      ){val=-9999;}// \todo [funct] pTransModel->GetSubBasinAvgConc(k,_svind)
        else if (_spaceAgg==BY_WSHED      ){val=-9999;}
        else if (_spaceAgg==BY_HRU_GROUP  ){val=-9999;}
        else if (_spaceAgg==BY_SELECT_HRUS){val=-9999;}*/
    }
    else if (_var==VAR_STATE_VAR){
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetStateVarValue(_svind);}
      else if (_spaceAgg==BY_BASIN      ){val=pModel->GetSubBasin (k)->GetAvgStateVar  (_svind);}
      else if (_spaceAgg==BY_WSHED      ){val=pModel->                 GetAvgStateVar  (_svind);}
      else if (_spaceAgg==BY_HRU_GROUP  ){val=pModel->GetHRUGroup (k)->GetAvgStateVar  (_svind);}
      else if (_spaceAgg==BY_SELECT_HRUS){val=pModel->GetHRUGroup (kk_only)->GetHRU(k)->GetStateVarValue(_svind);}
    }
    else if (_var==VAR_FORCING_FUNCTION){
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetForcing   (_force_str);}
      else if (_spaceAgg==BY_BASIN      ){val=pModel->GetSubBasin (k)->GetAvgForcing(_force_str);}
      else if (_spaceAgg==BY_WSHED      ){val=pModel->                 GetAvgForcing(_force_str);}
      else if (_spaceAgg==BY_HRU_GROUP  ){val=pModel->GetHRUGroup (k)->GetAvgForcing(_force_str);}
      else if (_spaceAgg==BY_SELECT_HRUS){val=pModel->GetHRUGroup (kk_only)->GetHRU(k)->GetForcing(_force_str);}
    }
    else if (_var == VAR_TO_FLUX){
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetCumulFlux   (_svind,true);}
      else if (_spaceAgg==BY_BASIN      ){val=pModel->GetSubBasin (k)->GetAvgCumulFlux(_svind,true);}
      else if (_spaceAgg==BY_WSHED      ){val=pModel->                 GetAvgCumulFlux(_svind,true);}
      else if (_spaceAgg==BY_HRU_GROUP  ){val=pModel->GetHRUGroup (k)->GetAvgCumulFlux(_svind,true);}
      else if (_spaceAgg==BY_SELECT_HRUS){val=pModel->GetHRUGroup (kk_only)->GetHRU(k)->GetCumulFlux(_svind,true);}
    }
    else if (_var == VAR_FROM_FLUX){
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetCumulFlux   (_svind,false);}
      else if (_spaceAgg==BY_BASIN      ){val=pModel->GetSubBasin (k)->GetAvgCumulFlux(_svind,false);}
      else if (_spaceAgg==BY_WSHED      ){val=pModel->                 GetAvgCumulFlux(_svind,false);}
      else if (_spaceAgg==BY_HRU_GROUP  ){val=pModel->GetHRUGroup (k)->GetAvgCumulFlux(_svind,false);}
      else if (_spaceAgg==BY_SELECT_HRUS){val=pModel->GetHRUGroup (kk_only)->GetHRU(k)->GetCumulFlux(_svind,false);}
    }
    else if (_var == VAR_BETWEEN_FLUX){
      if      (_spaceAgg==BY_HRU        ){val=pModel->GetHydroUnit(k)->GetCumulFluxBet(_svind,_svind2);}
      else if (_spaceAgg==BY_BASIN      ){val=pModel->GetSubBasin (k)->GetAvgCumulFluxBet(_svind,_svind2);}
      else if (_spaceAgg==BY_WSHED      ){val=pModel->                 GetAvgCumulFluxBet(_svind,_svind2);}
      else if (_spaceAgg==BY_HRU_GROUP  ){val=pModel->GetHRUGroup (k)->GetAvgCumulFluxBet(_svind,_svind2);}
      else if (_spaceAgg==BY_SELECT_HRUS){val=pModel->GetHRUGroup (kk_only)->GetHRU(k)->GetCumulFluxBet(_svind,_svind2);}
    }

    if (k==0){count++;}//increment number of data items stored

    //---Update diagnostics--------------------------------------------------
    //-----------------------------------------------------------------------
    if      (_aggstat==AGG_AVERAGE)
    {
      // \todo[funct] - should handle pointwise variables (e.g., state vars) differently from periodwise variables (e.g., forcings)
      if      (_var == VAR_STATE_VAR){
        data[k][0]=(double)(count-1)/(double)(count)*data[k][0]+val/count; //NOT CURRENTLY STRICTLY VALID
      }
      else {
        data[k][0]=(double)(count-1)/(double)(count)*data[k][0]+val/count;
      }
    }
    else if (_aggstat==AGG_MAXIMUM)
    {
      upperswap(data[k][0],val);
    }
    else if (_aggstat==AGG_MINIMUM)
    {
      lowerswap(data[k][0],val);
    }
    else if (_aggstat==AGG_RANGE)
    {
      lowerswap(data[k][0],val);
      upperswap(data[k][1],val);
    }
    else if ((_aggstat==AGG_MEDIAN)    ||
             (_aggstat==AGG_QUARTILES) ||
             (_aggstat==AGG_95CI)      ||
             (_aggstat==AGG_HISTOGRAM))
    {
      //populate data
      data [k][count-1] = val;
    }
    else
    {
      data[k][0]=val;//most recent value
    }

    //---Write output to file and reinitialize statistics, if needed---------
    //-----------------------------------------------------------------------
    if (reset)
    {
      if((Options.output_format==OUTPUT_STANDARD) || (Options.output_format==OUTPUT_ENSIM))
      {  
        string sep=",";
        if(Options.output_format==OUTPUT_ENSIM){ sep=" "; }//space separated
        //write to .csv or .tb0 file
        if     (_aggstat==AGG_AVERAGE){ _CUSTOM<<FormatDouble(data[k][0])      <<sep; }
        else if(_aggstat==AGG_MAXIMUM){ _CUSTOM<<FormatDouble(data[k][0])      <<sep; }
        else if(_aggstat==AGG_MINIMUM){ _CUSTOM<<FormatDouble(data[k][0])      <<sep; }
        else if(_aggstat==AGG_RANGE  ){ _CUSTOM<<FormatDouble(data[k][0])      <<sep<<FormatDouble(data[k][1])<<sep; }
        else if(_aggstat==AGG_MEDIAN)
        {
          double Q1,Q2,Q3;
          quickSort(data[k],0,count-1);
          GetQuartiles(data[k],count,Q1,Q2,Q3);
          _CUSTOM<<FormatDouble(Q2)<<sep;
        }
        else if(_aggstat==AGG_QUARTILES) //find lower quartile, median, then upper quartile
        {
          double Q1,Q2,Q3;
          quickSort(data[k],0,count-1);
          GetQuartiles(data[k],count,Q1,Q2,Q3);
          _CUSTOM<<FormatDouble(Q1)<<sep<<FormatDouble(Q2)<<sep<<FormatDouble(Q3)<<sep;
        }
        else if(_aggstat==AGG_95CI)
        {
          quickSort(data[k],0,count-1);//take floor and ceiling of lower and upper intervals to be conservative
          _CUSTOM  << FormatDouble(data[k][(int)floor((double)(count-1)*0.025)])<<sep;
          _CUSTOM  << FormatDouble(data[k][(int)ceil((double)(count-1)*0.975)])<<sep;
        }
        else if(_aggstat==AGG_HISTOGRAM)
        {
          quickSort(data[k],0,count-1);
          double binsize = (_hist_max-_hist_min)/_nBins;
          int bincount[MAX_HISTOGRAM_BINS];

          for(int bin=0;bin<_nBins;bin++){ bincount[bin]=0; }
          for(int a=0;a<count;a++)
          {
            for(int bin=0;bin<_nBins;bin++){
              if((data[k][a]>=(_hist_min+bin*binsize)) && (data[k][a]<(_hist_min+(bin+1)*binsize))){ bincount[bin]++; }
            }
          }

          for(int bin=0;bin<_nBins;bin++){
            _CUSTOM<<bincount[bin]<<sep;
          }
        }

      }
      else if(Options.output_format==OUTPUT_NETCDF)
      {
#ifdef _RVNETCDF_
        // Collect Data 
        double out;
        if     (_aggstat==AGG_AVERAGE){ out=data[k][0]; }
        else if(_aggstat==AGG_MAXIMUM){ out=data[k][0]; }
        else if(_aggstat==AGG_MINIMUM){ out=data[k][0]; }
        else if(_aggstat==AGG_MEDIAN)
        {
          double Q1,Q2,Q3;
          quickSort(data[k],0,count-1);
          GetQuartiles(data[k],count,Q1,Q2,Q3);
          out=Q2; 
        }
        output[k]=out;       
#endif
      }

      //-reset to initial conditions
      //-----------------------------------------------------------------------
      if      (_aggstat==AGG_AVERAGE){data[k][0]= 0.0;}
      else if (_aggstat==AGG_MAXIMUM){data[k][0]=-ALMOST_INF;}
      else if (_aggstat==AGG_MINIMUM){data[k][0]= ALMOST_INF;}
      else if (_aggstat==AGG_RANGE  )
      {
        data[k][0]= ALMOST_INF;
        data[k][1]=-ALMOST_INF;
      }
      else if ((_aggstat==AGG_MEDIAN)    ||
               (_aggstat==AGG_QUARTILES) ||
               (_aggstat==AGG_95CI)      ||
               (_aggstat==AGG_HISTOGRAM))
      {
        for(int j=0;j<num_store;j++){
          data[k][j]=0;
        }
      }
      //...more stats here
      if (k==num_data-1){count=0;}//don't reboot count until last HRU/basin is done
    }//end if (reset)
  }//end for (k=0;...

  if (reset){
  	if((Options.output_format==OUTPUT_STANDARD) || (Options.output_format==OUTPUT_ENSIM))
    { 
      _CUSTOM<<endl;//valid for ensim or .csv format
		}
		else if (Options.output_format==OUTPUT_NETCDF)
    {
#ifdef _RVNETCDF_
      // Write to NetCDF (done entire vector of data at once)
      if (num_data > 0)
      {
        int retval,data_id;
        size_t start2[2], count2[2]; 
        start2[0]=_time_index;  count2[0] = 1;         // writes exactly one time interval (yr/mo/day/hr)
        start2[1] = 0;          count2[1] = num_data;  // writes exactly num_data elements

        string netCDFtag=_statStr+"_"+_varName;
        retval = nc_inq_varid      (_netcdf_ID, netCDFtag.c_str(), &data_id);          HandleNetCDFErrors(retval);
        retval = nc_put_vara_double(_netcdf_ID, data_id, start2, count2, &output[0]);  HandleNetCDFErrors(retval);
      }
      delete [] output;
#endif
		}
    _time_index++;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Closes output stream after all information written to file
//
void CCustomOutput::CloseFiles()
{
  if (_CUSTOM.is_open()){_CUSTOM.close();}
}


