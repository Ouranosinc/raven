/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#include "Model.h"

/*****************************************************************
   Constructor/Destructor
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Implementation of the Model constructor
///
/// \param SM [in] Input soil model object
/// \param nsoillayers [in] Integer number of soil layers
//
CModel::CModel(const soil_model SM,
               const int        nsoillayers,
               const optStruct &Options)
{
  int i;
  _nSubBasins=0;    _pSubBasins=NULL;
  _nHydroUnits=0;   _pHydroUnits=NULL;
  _nHRUGroups=0;    _pHRUGroups=NULL;
  _nGauges=0;       _pGauges=NULL;
  _nForcingGrids=0; _pForcingGrids=NULL;
  _nProcesses=0;    _pProcesses=NULL;
  _nCustomOutputs=0;_pCustomOutputs=NULL;
  _nTransParams=0;  _pTransParams=NULL;
  _nClassChanges=0; _pClassChanges=NULL;
  _nObservedTS=0;   _pObservedTS=NULL; _pModeledTS=NULL; _aObsIndex=NULL;
  _nObsWeightTS =0; _pObsWeightTS=NULL;
  _nDiagnostics=0;  _pDiagnostics=NULL;

  _nTotalConnections=0;

  _WatershedArea=0;

  _aSubBasinOrder =NULL; _maxSubBasinOrder=0;
  _aOrderedSBind  =NULL;
  _aDownstreamInds=NULL;

  _pOptStruct = &Options;

  _HYDRO_ncid=-9; 
  _STORAGE_ncid=-9;
  _FORCINGS_ncid=-9; 

  ExitGracefullyIf(nsoillayers<1,
                   "CModel constructor::improper number of soil layers. SoilModel not specified?",BAD_DATA);

  //Initialize Lookup table for state variable indices
  for (int i=0;i<MAX_STATE_VARS;i++){
    for (int m=0;m<MAX_SV_LAYERS;m++){
      _aStateVarIndices[i][m]=DOESNT_EXIST;
    }
  }
  //determine first group of state variables based upon soil model
  //SW, atmosphere, atmos_precip always present, one for each soil layer, and 1 for GW (unless lumped)
  _nStateVars=4+nsoillayers;

  _aStateVarType =new sv_type [_nStateVars];
  _aStateVarLayer=new int     [_nStateVars];
  _nSoilVars     =nsoillayers;
  _nAquiferLayers=0;  //initially
  _nSnowLayers   =0;  //initially

  _aStateVarType[0]=SURFACE_WATER; _aStateVarLayer[0]=DOESNT_EXIST; _aStateVarIndices[(int)(SURFACE_WATER)][0]=0;
  _aStateVarType[1]=ATMOSPHERE;    _aStateVarLayer[1]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOSPHERE   )][0]=1;
  _aStateVarType[2]=ATMOS_PRECIP;  _aStateVarLayer[2]=DOESNT_EXIST; _aStateVarIndices[(int)(ATMOS_PRECIP )][0]=2;
  _aStateVarType[3]=PONDED_WATER;  _aStateVarLayer[3]=DOESNT_EXIST; _aStateVarIndices[(int)(PONDED_WATER )][0]=3;

  int count=0;
  for (i=4;i<4+_nSoilVars;i++)
  {
    _aStateVarType [i]=SOIL;
    _aStateVarLayer[i]=count;
    _aStateVarIndices[(int)(SOIL)][count]=i;
    count++;
  }
  _lake_sv=0; //by default, rain on lake goes direct to surface storage [0]


  CHydroProcessABC::SetModel(this);
  CLateralExchangeProcessABC::SetModel(this);

  _aGaugeWeights =NULL; //Initialized in Initialize
  _aGaugeWtTemp  =NULL;
  _aGaugeWtPrecip=NULL;

  _aCumulativeBal   =NULL;
  _aFlowBal         =NULL;
  _aCumulativeLatBal=NULL;
  _aFlowLatBal      =NULL;
  _CumulInput       =0.0;
  _CumulOutput      =0.0;
  _CumEnergyGain    =0.0;
  _CumEnergyLoss    =0.0;
  _initWater        =0.0;

  _UTM_zone=-1;

  _nOutputTimes=0;   _aOutputTimes=NULL;
  _currOutputTimeInd=0;
  _pOutputGroup=NULL;

  _aShouldApplyProcess=NULL; //Initialized in Initialize

  _pTransModel=new CTransportModel(this);
}

/////////////////////////////////////////////////////////////////
/// \brief Implementation of the default destructor
//
CModel::~CModel()
{
  if (DESTRUCTOR_DEBUG){cout<<"DELETING MODEL"<<endl;}
  int c,f,g,i,j,k,kk,p;

  CloseOutputStreams();

  for (p=0;p<_nSubBasins;    p++){delete _pSubBasins    [p];} delete [] _pSubBasins;    _pSubBasins=NULL;
  for (k=0;k<_nHydroUnits;   k++){delete _pHydroUnits   [k];} delete [] _pHydroUnits;   _pHydroUnits=NULL;
  for (g=0;g<_nGauges;       g++){delete _pGauges       [g];} delete [] _pGauges;       _pGauges=NULL;
  for (f=0;f<_nForcingGrids; f++){delete _pForcingGrids [f];} delete [] _pForcingGrids; _pForcingGrids=NULL;
  for (j=0;j<_nProcesses;    j++){delete _pProcesses    [j];} delete [] _pProcesses;    _pProcesses=NULL;
  for (c=0;c<_nCustomOutputs;c++){delete _pCustomOutputs[c];} delete [] _pCustomOutputs;_pCustomOutputs=NULL;
  for (i=0;i<_nObservedTS;   i++){delete _pObservedTS   [i];} delete [] _pObservedTS;   _pObservedTS=NULL;
  if (_pModeledTS != NULL){
    for (i = 0; i < _nObservedTS; i++){ delete _pModeledTS[i]; } delete[] _pModeledTS;    _pModeledTS = NULL;
  }
  for (i=0;i<_nObsWeightTS;  i++){delete _pObsWeightTS  [i];} delete [] _pObsWeightTS;  _pObsWeightTS=NULL;
  for (j=0;j<_nDiagnostics;  j++){delete _pDiagnostics  [j];} delete [] _pDiagnostics;  _pDiagnostics=NULL;

  if (_aCumulativeBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aCumulativeBal[k];} delete [] _aCumulativeBal; _aCumulativeBal=NULL;
  }
  if (_aFlowBal!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aFlowBal[k];      } delete [] _aFlowBal;       _aFlowBal=NULL;
  }
  if (_aCumulativeLatBal!=NULL){delete [] _aCumulativeLatBal; _aCumulativeLatBal=NULL;}
  if (_aFlowLatBal      !=NULL){delete [] _aFlowLatBal;       _aFlowLatBal=NULL;}
  if (_aGaugeWeights!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWeights[k]; } delete [] _aGaugeWeights;  _aGaugeWeights=NULL;
  }
  if (_aGaugeWtPrecip!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtPrecip[k]; } delete [] _aGaugeWtPrecip;  _aGaugeWtPrecip=NULL;
  }
  if (_aGaugeWtTemp!=NULL){
    for (k=0;k<_nHydroUnits;   k++){delete [] _aGaugeWtTemp[k]; } delete [] _aGaugeWtTemp;  _aGaugeWtTemp=NULL;
  }
  if (_aShouldApplyProcess!=NULL){
    for (k=0;k<_nProcesses;   k++){delete [] _aShouldApplyProcess[k]; } delete [] _aShouldApplyProcess;  _aShouldApplyProcess=NULL;
  }
  for (kk=0;kk<_nHRUGroups;kk++){delete _pHRUGroups[kk]; } delete [] _pHRUGroups;   _pHRUGroups  =NULL;
  for (j=0;j<_nTransParams;j++) {delete _pTransParams[j];} delete [] _pTransParams; _pTransParams=NULL;
  for (j=0;j<_nClassChanges;j++){delete _pClassChanges[j];} delete [] _pClassChanges; _pClassChanges=NULL;

  delete [] _aStateVarType;  _aStateVarType=NULL;
  delete [] _aStateVarLayer; _aStateVarLayer=NULL;
  delete [] _aSubBasinOrder; _aSubBasinOrder=NULL;
  delete [] _aOrderedSBind;  _aOrderedSBind=NULL;
  delete [] _aDownstreamInds;_aDownstreamInds=NULL;
  delete [] _aOutputTimes;   _aOutputTimes=NULL;
  delete [] _aObsIndex;      _aObsIndex=NULL;


  CSoilClass::      DestroyAllSoilClasses();
  CVegetationClass::DestroyAllVegClasses();
  CLandUseClass::   DestroyAllLUClasses();
  CTerrainClass::   DestroyAllTerrainClasses();
  CSoilProfile::    DestroyAllSoilProfiles();
  CAquiferStack::   DestroyAllAqStacks();
  CChannelXSect::   DestroyAllChannelXSections();

  delete _pTransModel;
}
/*****************************************************************
   Accessors
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns number of sub basins in model
///
/// \return Integer number of sub basins
//
int CModel::GetNumSubBasins   () const{return _nSubBasins;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of HRUs in model
///
/// \return Integer number of HRUs
//
int CModel::GetNumHRUs        () const{return _nHydroUnits;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of HRU groups
///
/// \return Integer number of HRU groups
//
int CModel::GetNumHRUGroups   () const{return _nHRUGroups;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of gauges in model
///
/// \return Integer number of gauges in model
//
int CModel::GetNumGauges      () const{return _nGauges;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of gridded forcings in model
///
/// \return Integer number of gridded forcings in model
//
int CModel::GetNumForcingGrids () const{return _nForcingGrids;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of state variables per HRU in model
///
/// \return Integer number of state variables per HRU in model
//
int CModel::GetNumStateVars   () const{return _nStateVars;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of soil layers
///
/// \return Integer number of soil layers
//
int CModel::GetNumSoilLayers  () const{return _nSoilVars;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of aquifer layers
///
/// \return Integer number of aquifer layers
//
int CModel::GetNumAquiferLayers  () const{return _nAquiferLayers;}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of hydrologic processes simulated by model
///
/// \return Integer number of hydrologic processes simulated by model
//
int CModel::GetNumProcesses   () const{return _nProcesses;}

//////////////////////////////////////////////////////////////////
/// \brief Returns total modeled watershed area
///
/// \return total modeled watershed area [km2]
//
double CModel::GetWatershedArea () const{return _WatershedArea;}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific hydrologic process denoted by parameter
///
/// \param j [in] Process index
/// \return pointer to hydrologic process corresponding to passed index j
//
CHydroProcessABC *CModel::GetProcess(const int j) const
{

#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetProcess::improper index",BAD_DATA);
#endif
  return _pProcesses[j];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific gauge denoted by index
///
/// \param g [in] Gauge index
/// \return pointer to gauge corresponding to passed index g
//
CGauge *CModel::GetGauge(const int g) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((g<0) || (g>=_nGauges),"CModel GetGauge::improper index",BAD_DATA);
#endif
  return _pGauges[g];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific forcing grid denoted by index
///
/// \param f [in] Forcing Grid index
/// \return pointer to gauge corresponding to passed index g
//
CForcingGrid *CModel::GetForcingGrid(const forcing_type &ftyp) const
{
  int f=GetForcingGridIndexFromType(ftyp);
#ifdef _STRICTCHECK_
  ExitGracefullyIf((f<0) || (f>=_nForcingGrids),"CModel GetForcingGrid::improper index",BAD_DATA);
#endif
  return _pForcingGrids[f];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU denoted by index k
///
/// \param k [in] HRU index
/// \return pointer to HRU corresponding to passed index k
//
CHydroUnit *CModel::GetHydroUnit(const int k) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel GetHydroUnit::improper index",BAD_DATA);
#endif
  return _pHydroUnits[k];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU with HRU identifier HRUID
///
/// \param HRUID [in] HRU identifier
/// \return pointer to HRU corresponding to passed ID HRUID, NULL if no such HRU exists
//
CHydroUnit *CModel::GetHRUByID(const int HRUID) const
{
  for (int k=0;k<_nHydroUnits;k++){
    if (HRUID==_pHydroUnits[k]->GetID()){return _pHydroUnits[k];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU group denoted by parameter kk
///
/// \param kk [in] HRU group index
/// \return pointer to HRU group corresponding to passed index kk
//
CHRUGroup  *CModel::GetHRUGroup(const int kk) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((kk<0) || (kk>=_nHRUGroups),"CModel GetHRUGroup::improper index",BAD_DATA);
#endif
  return _pHRUGroups[kk];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific HRU group denoted by string parameter
///
/// \param name [in] String name of HRU group
/// \return pointer to HRU group corresponding to passed name, or NULL if this group doesn't exist
//
CHRUGroup  *CModel::GetHRUGroup(const string name) const
{
  for (int kk=0;kk<_nHRUGroups;kk++){
    if (!name.compare(_pHRUGroups[kk]->GetName())){
      return _pHRUGroups[kk];
    }
  }
  return NULL;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns true if HRU with global index k is in specified HRU Group
///
/// \param k [in] HRU global index
/// \param HRUGroupName [in] String name of HRU group
/// \return true if HRU k is in HRU Group specified by HRUGroupName
//
bool              CModel::IsInHRUGroup(const int k, const string HRUGroupName) const{

  CHRUGroup *pGrp=NULL;
  pGrp=GetHRUGroup(HRUGroupName);
  if (pGrp == NULL){ return false; }//throw warning?

  int kk = pGrp->GetGlobalIndex();
  for (int k_loc=0; k_loc<_pHRUGroups[kk]->GetNumHRUs(); k_loc++)
  {
    if (_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex()==k){return true;}
  }
  return false;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns specific Sub basin denoted by index parameter
///
/// \param p [in] Sub basin index
/// \return pointer to the Sub basin object corresponding to passed index
//
CSubBasin  *CModel::GetSubBasin(const int p) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((p<0) || (p>=_nSubBasins),"CModel GetSubBasin::improper index",BAD_DATA);
#endif
  return _pSubBasins[p];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of subbasin downstream from subbasin referred to by index
///
/// \param p [in] List index for accessing subbasin
/// \return downstream subbasin index, if input index is valid; -1 if there is no downstream basin
//
int         CModel::GetDownstreamBasin(const int p) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((p<0) || (p>=_nSubBasins),"GetDownstreamBasin: Invalid index",BAD_DATA);
#endif
  return _aDownstreamInds[p];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns subbasin object corresponding to passed subbasin ID
///
/// \param SBID [in] Integer sub basin ID
/// \return pointer to Sub basin object corresponding to passed ID, if ID is valid
//
CSubBasin  *CModel::GetSubBasinByID(const long SBID) const
{ //could be quite slow...
  if (SBID<0){return NULL;}
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->GetID()==SBID){return _pSubBasins[p];}
  }
  return NULL;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns sub basin index corresponding to passed subbasin ID
///
/// \param SBID [in] Integer subbasin ID
/// \return Sub basin index corresponding to passed ID, if ID is valid
//
int         CModel::GetSubBasinIndex(const long SBID) const
{  //could be quite slow...
  if (SBID<0){return DOESNT_EXIST;}
  for (int p=0;p<_nSubBasins;p++){
    if (_pSubBasins[p]->GetID()==SBID){return p;}
  }
  return INDEX_NOT_FOUND;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns hydrologic process type corresponding to passed index
///
/// \param j [in] Integer index
/// \return Process type corresponding to process with passed index
//
process_type CModel::GetProcessType(const int j) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetProcessType::improper index",BAD_DATA);
#endif
  return _pProcesses[j]->GetProcessType();
}

//////////////////////////////////////////////////////////////////
/// \brief Returns number of connections of hydrological process associated with passed index
///
/// \param j [in] Integer index corresponding to a hydrological process
/// \return Number of connections associated with hydrological process symbolized by index
/// \note should be called only by solver
//
int         CModel::GetNumConnections (const int j) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) || (j>=_nProcesses),"CModel GetNumConnections::improper index",BAD_DATA);
#endif
  return _pProcesses[j]->GetNumConnections();
}

//////////////////////////////////////////////////////////////////
/// \brief Returns state variable type corresponding to passed state variable array index
///
/// \param i [in] state variable array index (>=0, <nStateVariables)
/// \return State variable type corresponding to passed state variable array index
//
sv_type     CModel::GetStateVarType(const int i) const
{
#ifdef _STRICTCHECK_
  string warn="CModel GetStateVarType::improper index ("+to_string(i)+")";
  ExitGracefullyIf((i<0) || (i>=_nStateVars),warn.c_str(),BAD_DATA);
#endif
  return _aStateVarType[i];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of state variable type passed
///
/// \param type [in] State variable type
/// \return Index which corresponds to the state variable type passed, if it exists; DOESNT_EXIST (-1) otherwise
/// \note  should only be used for state variable types without multiple levels; issues for soils, e.g.
//
int         CModel::GetStateVarIndex(sv_type type) const
{
  return _aStateVarIndices[(int)(type)][0];
}

//////////////////////////////////////////////////////////////////
/// \brief Returns index of state variable type passed (for repeated state variables)
///
/// \param type [in] State variable type
/// \param layer [in] Integer identifier of the layer of interest (or, possibly DOESNT_EXIST if variable doesn't have layers)
/// \return Index which corresponds to the state variable type passed, if it exists; DOESNT_EXIST (-1) otherwise
//
int         CModel::GetStateVarIndex(sv_type type, int layer) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((layer!=DOESNT_EXIST) && ((layer<0) || (layer>=MAX_SV_LAYERS)),
                   "CModel GetStateVarIndex::improper layer",BAD_DATA);
#endif
  if (layer==DOESNT_EXIST){return _aStateVarIndices[(int)(type)][0];    }
  else                    {return _aStateVarIndices[(int)(type)][layer];}
}

//////////////////////////////////////////////////////////////////
/// \brief Uses state variable type index to access index of layer to which it corresponds
///
/// \param ii [in] Index referencing a type of state variable
/// \return Index which corresponds to soil layer, or 0 for a non-layered state variable
//
int         CModel::GetStateVarLayer(const int ii) const
{
  int count=0;
  for (int i=0;i<ii;i++){
    if (_aStateVarType[i]==_aStateVarType[ii]){count++;}
  }
  return count;
}

//////////////////////////////////////////////////////////////////
/// \brief Checks if state variable passed exists in model
///
/// \param typ [in] State variable type
/// \return Boolean indicating whether state variable exists in model
//
bool        CModel::StateVarExists(sv_type typ) const
{
  return (GetStateVarIndex(typ)!=DOESNT_EXIST);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns lake storage variable index
///
/// \return Integer index of lake storage variable
//
int         CModel::GetLakeStorageIndex() const{return _lake_sv;}

//////////////////////////////////////////////////////////////////
/// \brief Returns gauge index of gauge with specified name
///
/// \return Integer index of gauge
/// \param name [in] specified name
//
int  CModel::GetGaugeIndexFromName (const string name) const{
  for (int g=0;g<_nGauges;g++){
    if (name==_pGauges[g]->GetName()){return g;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns forcing grid index of forcing grid with specified name
///
/// \return Integer index of forcing grid
/// \param name [in] specified name
//
/*int  CModel::GetForcingGridIndexFromName (const string name) const{
  for (int f=0;f<_nForcingGrids;f++){
    if (name==ForcingToString(_pForcingGrids[f]->GetName())){return f;}
  }
  return DOESNT_EXIST;
}*/
//////////////////////////////////////////////////////////////////
/// \brief Returns forcing grid index of forcing grid with specified type
///
/// \return Integer index of forcing grid
/// \param name [in] specified type
//
int  CModel::GetForcingGridIndexFromType (const forcing_type &typ) const{
  for (int f=0;f<_nForcingGrids;f++){
    if (typ==_pForcingGrids[f]->GetName()){return f;}
  }
  return DOESNT_EXIST;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flux (mm/d, MJ/m2/d, mg/m2/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param iFrom_test [in] index of "from" storage state variable
/// \param iTo_test [in] index of "to" storage state variable
/// \param &Options [in] Global model options information
/// \todo [optimize]
//
/*double CModel::GetFlux(const int k, const int iFrom_test, const int iTo_test, const optStruct &Options) const
{
  int q,j,js(0);
  double flow=0;
  const int *iFrom,*iTo;
  int nConnections;

  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetFlux: bad HRU index",RUNTIME_ERR);

  //Not Optimized
  //much more effective to store aFlowBal as aFlowRate[k][iFrom][iTo]?
  if (iFrom_test==iTo_test){return 0.0;}
  for (j=0;j<_nProcesses;j++)
  {
    iTo  =_pProcesses[j]->GetToIndices();
    iFrom=_pProcesses[j]->GetFromIndices();
    nConnections=_pProcesses[j]->GetNumConnections();
    for (q=0;q<nConnections;q++)
    {
      if ((iTo_test==iTo[q]) && (iFrom_test==iFrom[q]))
      {
        flow+=_aFlowBal[k][js]/Options.timestep;
      }
      if ((iTo_test==iFrom[q]) && (iFrom_test==iTo[q]))
      {
        flow-=_aFlowBal[k][js]/Options.timestep;
      }
      js++;
    }//end for q=0 to nConnections
  }//end for j=0 to nProcesses

  return flow;
}*/
//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flux (mm/d, MJ/m2/d, mg/m2/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param js [in] index of process connection (i.e., j*)
/// \param &Options [in] Global model options information
//
double CModel::GetFlux(const int k, const int js, const optStruct &Options) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((js<0) || (js>=_nTotalConnections),"CModel::GetFlux: bad connection index",RUNTIME_ERR);
#endif
  return _aFlowBal[k][js]/Options.timestep;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns current mass/energy flow (mm-m2/d, MJ/d, mg/d) between two storage compartments iFrom and iTo
/// \details required for advective transport processes
///
/// \param k [in] HRU index
/// \param qs [in] global index of process connection (i.e., q*)
/// \param &Options [in] Global model options information
//
double CModel::GetLatFlow(const int qs,const optStruct &Options) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((qs<0) || (qs>=_nTotalConnections),"CModel::GetFlux: bad connection index",RUNTIME_ERR);
#endif
  return _aFlowLatBal[qs]/Options.timestep;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative flux to/from storage unit i
///
/// \param k [in] index of HRU
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return cumulative flux to storage compartment i in hru K
//
double CModel::GetCumulativeFlux(const int k, const int i, const bool to) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetCumulativeFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((i<0) || (i>=_nStateVars),"CModel::GetCumulativeFlux: bad state var index",RUNTIME_ERR);
#endif
  int js=0;
  double sum=0;
  double area=_pHydroUnits[k]->GetArea();
  int jss=0;
  for(int j = 0; j < _nProcesses; j++)
  {
    for(int q = 0; q < _pProcesses[j]->GetNumConnections(); q++)//each process may have multiple connections
    {
      if(( to) && (_pProcesses[j]->GetToIndices()[q]   == i)){ sum+=_aCumulativeBal[k][js]; }
      if((!to) && (_pProcesses[j]->GetFromIndices()[q] == i)){ sum+=_aCumulativeBal[k][js]; }
      js++;
    }

    if(_pProcesses[j]->GetNumLatConnections()>0)
    {
      CLateralExchangeProcessABC *pProc=(CLateralExchangeProcessABC*)_pProcesses[j];
      for(int q = 0; q < _pProcesses[j]->GetNumLatConnections(); q++)//each process may have multiple connections
      {
        if(( to) && (pProc->GetToHRUIndices()[q]  ==k) && (pProc->GetLateralToIndices()[q]  ==i)){sum+=_aCumulativeLatBal[jss]/area; }
        if((!to) && (pProc->GetFromHRUIndices()[q]==k) && (pProc->GetLateralFromIndices()[q]==i)){sum+=_aCumulativeLatBal[jss]/area; }
        jss++;
      }
    }
  }
  
  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns cumulative gross flux between unit iFrom and iTo in HRU k
///
/// \param k [in] index of HRU
/// \param iFrom [in] index of storage compartment
/// \param iTo [in] index of storage compartment
/// \return cumulative gross flux between unit iFrom and iTo in HRU k
//  does not address fluxes due to lateral flux
double CModel::GetCumulFluxBetween(const int k,const int iFrom,const int iTo) const
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((k<0) || (k>=_nHydroUnits),"CModel::GetCumulativeFlux: bad HRU index",RUNTIME_ERR);
  ExitGracefullyIf((iFrom<0) || (iTo>=_nStateVars),"CModel::GetCumulativeFlux: bad state var index",RUNTIME_ERR);
#endif
  int js=0;
  double sum=0;
  int q;
  const int *iFromp; 
  const int *iTop;
  int nConn;
  for(int j = 0; j < _nProcesses; j++)
  {
    iFromp=_pProcesses[j]->GetFromIndices();
    iTop  =_pProcesses[j]->GetToIndices();
    nConn =_pProcesses[j]->GetNumConnections();
    for (q = 0; q < nConn; q++)//each process may have multiple connections
    {
      if( (iTop  [q]== iTo) && (iFromp[q]== iFrom)){ sum+=_aCumulativeBal[k][js]; }
      if( (iFromp[q]== iTo) && (iTop  [q]== iFrom)){ sum-=_aCumulativeBal[k][js]; }
      js++;
    }
  }
  
  return sum;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns options structure model
/// \return pointer to transport model
//
const optStruct  *CModel::GetOptStruct() const{ return _pOptStruct; }

//////////////////////////////////////////////////////////////////
/// \brief Returns transport model
/// \return pointer to transport model
//
CTransportModel  *CModel::GetTransportModel() const{return _pTransModel;}

/*****************************************************************
   Watershed Diagnostic Functions
    -aggregate data from subbasins and HRUs
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average total precipitation rate at all HRUs [mm/d]
///
/// \return Area-weighted average of total precipitation rate [mm/d] over all HRUs
//
double CModel::GetAveragePrecip() const
{
  double sum(0);
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=_pHydroUnits[k]->GetForcingFunctions()->precip*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current area-weighted average snowfall rate [mm/d] at all HRUs
///
/// \return Current area-weighted average snowfall rate [mm/d] at all HRUs
//
double CModel::GetAverageSnowfall() const
{
  double sum(0);
  const force_struct *f;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      f=_pHydroUnits[k]->GetForcingFunctions();
      sum+=(f->precip*f->snow_frac)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns current area-weighted average forcing functions at all HRUs
///
/// \return Current area-weighted average forcing functions at all HRUs
//
force_struct CModel::GetAverageForcings() const
{
  static force_struct Fave;
  const force_struct *pF_hru;
  double area_wt;
  ZeroOutForcings(Fave);

  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      pF_hru =_pHydroUnits[k]->GetForcingFunctions();
      area_wt=_pHydroUnits[k]->GetArea()/_WatershedArea;

      Fave.precip          +=area_wt*pF_hru->precip;
      Fave.precip_daily_ave+=area_wt*pF_hru->precip_daily_ave;
      Fave.precip_5day     +=area_wt*pF_hru->precip_5day;
      Fave.snow_frac       +=area_wt*pF_hru->snow_frac;

      Fave.temp_ave       +=area_wt*pF_hru->temp_ave;
      Fave.temp_daily_min +=area_wt*pF_hru->temp_daily_min;
      Fave.temp_daily_max +=area_wt*pF_hru->temp_daily_max;
      Fave.temp_daily_ave +=area_wt*pF_hru->temp_daily_ave;
      Fave.temp_month_max +=area_wt*pF_hru->temp_month_max;
      Fave.temp_month_min +=area_wt*pF_hru->temp_month_min;
      Fave.temp_month_ave +=area_wt*pF_hru->temp_month_ave;

      Fave.temp_ave_unc   +=area_wt*pF_hru->temp_ave_unc;
      Fave.temp_min_unc   +=area_wt*pF_hru->temp_min_unc;
      Fave.temp_max_unc   +=area_wt*pF_hru->temp_max_unc;

      Fave.air_dens       +=area_wt*pF_hru->air_dens;
      Fave.air_pres       +=area_wt*pF_hru->air_pres;
      Fave.rel_humidity   +=area_wt*pF_hru->rel_humidity;

      Fave.cloud_cover    +=area_wt*pF_hru->cloud_cover;
      Fave.ET_radia       +=area_wt*pF_hru->ET_radia;
      Fave.SW_radia       +=area_wt*pF_hru->SW_radia;
      Fave.SW_radia_unc   +=area_wt*pF_hru->SW_radia_unc;
      Fave.SW_radia_net   +=area_wt*pF_hru->SW_radia_net;
      Fave.LW_radia       +=area_wt*pF_hru->LW_radia;
      Fave.day_length     +=area_wt*pF_hru->day_length;
      Fave.day_angle      +=area_wt*pF_hru->day_angle;   //not really necc.

      Fave.wind_vel       +=area_wt*pF_hru->wind_vel;

      Fave.PET            +=area_wt*pF_hru->PET;
      Fave.OW_PET         +=area_wt*pF_hru->OW_PET;
      Fave.PET_month_ave  +=area_wt*pF_hru->PET_month_ave;

      Fave.potential_melt +=area_wt*pF_hru->potential_melt;

      Fave.recharge       +=area_wt*pF_hru->recharge;

      Fave.subdaily_corr  +=area_wt*pF_hru->subdaily_corr;
    }
  }
  return Fave;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns area weighted average of state variables across all modeled HRUs
///
/// \param i [in] State variable index
/// \return Area weighted average of state variables across all modeled HRUs
//
double CModel::GetAvgStateVar(const int i) const
{
  //Area-weighted average
#ifdef _STRICTCHECK_
  ExitGracefullyIf((i<0) || (i>=_nStateVars),"CModel GetAvgStateVar::improper index",BAD_DATA);
#endif
  double sum(0.0);
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum+=(_pHydroUnits[k]->GetStateVarValue(i)*_pHydroUnits[k]->GetArea());
    }
  }
  return sum/_WatershedArea;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified forcing function over watershed
///
/// \param &forcing_string [in] string identifier of forcing function to assess
/// \return Area-weighted average of specified forcing function
//
double CModel::GetAvgForcing (const string &forcing_string) const
{
  //Area-weighted average
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if (_pHydroUnits[k]->IsEnabled())
    {
      sum    +=_pHydroUnits[k]->GetForcing(forcing_string)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of specified cumulative flux over watershed
///
/// \param i [in] index of storage compartment
/// \param to [in] true if evaluating cumulative flux to storage compartment, false for 'from'
/// \return Area-weighted average of cumulative flux to storage compartment i
//
double CModel::GetAvgCumulFlux (const int i, const bool to) const
{
  //Area-weighted average
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum +=GetCumulativeFlux(k,i,to)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns area-weighted average of  cumulative flux between two compartments over watershed
///
/// \param iFrom [in] index of 'from' storage compartment
/// \param iTo [in] index of 'to' storage compartment
/// \return Area-weighted average of cumulative flux between two compartments over watershed
//
double CModel::GetAvgCumulFluxBet (const int iFrom, const int iTo) const
{
  //Area-weighted average
  double sum=0.0;
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      sum +=GetCumulFluxBetween(k,iFrom,iTo)*_pHydroUnits[k]->GetArea();
    }
  }
  return sum/_WatershedArea;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total channel storage [mm]
/// \return Total channel storage in all of watershed [mm]
//
double CModel::GetTotalChannelStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetChannelStorage();            //[m3]
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total reservoir storage [mm]
/// \return Total reservoir storage in all of watershed [mm]
//
double CModel::GetTotalReservoirStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetReservoirStorage();          //[m3]
  }
  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}
//////////////////////////////////////////////////////////////////
/// \brief Returns total rivulet storage distributed over watershed [mm]
/// \return Total rivulet storage distributed over watershed [mm]
//
double CModel::GetTotalRivuletStorage() const
{
  double sum(0);

  for (int p=0;p<_nSubBasins;p++)
  {
    sum+=_pSubBasins[p]->GetRivuletStorage();//[m3]
  }

  return sum/(_WatershedArea*M2_PER_KM2)*MM_PER_METER;
}

/*****************************************************************
   Model Creation Functions
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Adds additional HRU to model
///
/// \param *pHRU [in] (valid) pointer to HRU to be added
//
void CModel::AddHRU(CHydroUnit *pHRU)
{
  if (!DynArrayAppend((void**&)(_pHydroUnits),(void*)(pHRU),_nHydroUnits)){
    ExitGracefully("CModel::AddHRU: adding NULL HRU",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds HRU group
///
/// \param *pHRUGroup [in] (valid) pointer to HRU group to be added
//
void CModel::AddHRUGroup(CHRUGroup *pHRUGroup)
{
  for(int kk=0;kk<_nHRUGroups;kk++){
    if(pHRUGroup->GetName()==_pHRUGroups[kk]->GetName()){
      WriteWarning("CModel::AddHRUGroups: cannot add two HRU groups with the same name. Group "+pHRUGroup->GetName()+ " is duplicated in input.",true);
      break;
    }
  }
  if (!DynArrayAppend((void**&)(_pHRUGroups),(void*)(pHRUGroup),_nHRUGroups)){
    ExitGracefully("CModel::AddHRUGroup: adding NULL HRU Group",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds sub basin
///
/// \param *pSB [in] (valid) pointer to Sub basin to be added
//
void CModel::AddSubBasin(CSubBasin *pSB)
{
  if (!DynArrayAppend((void**&)(_pSubBasins),(void*)(pSB),_nSubBasins)){
    ExitGracefully("CModel::AddSubBasin: adding NULL HRU",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds gauge to model
///
/// \param *pGage [in] (valid) pointer to Gauge to be added to model
//
void CModel::AddGauge  (CGauge *pGage)
{
  if (!DynArrayAppend((void**&)(_pGauges),(void*)(pGage),_nGauges)){
    ExitGracefully("CModel::AddGauge: adding NULL Gauge",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds gridded forcing to model
///
/// \param *pGrid [in] (valid) pointer to ForcingGrid to be added to model
//
void CModel::AddForcingGrid  (CForcingGrid *pGrid)
{
  if (!DynArrayAppend((void**&)(_pForcingGrids),(void*)(pGrid),_nForcingGrids)){
    ExitGracefully("CModel::AddForcingGrid: adding NULL ForcingGrid",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds transient parameter to model
///
/// \param *pTP [in] (valid) pointer to transient parameter to be added to model
//
void CModel::AddTransientParameter(CTransientParam   *pTP)
{
  if (!DynArrayAppend((void**&)(_pTransParams),(void*)(pTP),_nTransParams)){
    ExitGracefully("CModel::AddTransientParameter: adding NULL transient parameter",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds class change to model
///
/// \param *pTP [in] (valid) pointer to transient parameter to be added to model
//
void CModel::AddPropertyClassChange(const string HRUgroup,
                                    const class_type tclass,
                                    const string new_class,
                                    const time_struct &tt)
{
  class_change *pCC=NULL;
  pCC=new class_change;
  pCC->HRU_groupID=DOESNT_EXIST;
  for (int kk = 0; kk < _nHRUGroups; kk++){
    if (!strcmp(_pHRUGroups[kk]->GetName().c_str(), HRUgroup.c_str()))
    {
      pCC->HRU_groupID=kk;
    }
  }
  if (pCC->HRU_groupID == DOESNT_EXIST){
    string warning = "CModel::AddPropertyClassChange: invalid HRU Group name: " + HRUgroup+ ". HRU group names should be defined in .rvi file using :DefineHRUGroups command. ";
    ExitGracefullyIf(pCC->HRU_groupID == DOESNT_EXIST,warning.c_str(),BAD_DATA_WARN); return;
  }
  pCC->newclass=new_class;
  if ((tclass == CLASS_LANDUSE) && (CLandUseClass::StringToLUClass(new_class) == NULL)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid land use class specified",BAD_DATA_WARN);return;
  }
  if ((tclass == CLASS_VEGETATION) && (CVegetationClass::StringToVegClass(new_class) == NULL)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid vegetation class specified",BAD_DATA_WARN);return;
  }
  if ((tclass == CLASS_HRUTYPE) && (StringToHRUType(new_class) == HRU_INVALID_TYPE)){
    ExitGracefully("CModel::AddPropertyClassChange: invalid HRU type specified",BAD_DATA_WARN);return;
  }

  pCC->tclass=tclass;
  if ((tclass != CLASS_VEGETATION) && (tclass != CLASS_LANDUSE) && (tclass!=CLASS_HRUTYPE)){
    ExitGracefully("CModel::AddPropertyClassChange: only vegetation, land use, and HRU type classes may be changed during the course of simulation",BAD_DATA_WARN);return;
  }

  //convert time to model time
  pCC->modeltime= TimeDifference(_pOptStruct->julian_start_day,_pOptStruct->julian_start_year ,tt.julian_day, tt.year);
  if ((pCC->modeltime<0) || (pCC->modeltime>_pOptStruct->duration)){
    string warn;
    warn="Property Class change dated "+tt.date_string+" does not occur during model simulation time";
    WriteWarning(warn,_pOptStruct->noisy);
  }

  //cout << "PROPERTY CLASS CHANGE " << pCC->HRU_groupID << " " << pCC->tclass << " "<<pCC->modeltime<<endl;
  if (!DynArrayAppend((void**&)(_pClassChanges),(void*)(pCC),_nClassChanges)){
    ExitGracefully("CModel::AddPropertyClassChange: adding NULL property class change",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds observed time series to model
///
/// \param *pTS [in] (valid) pointer to observed time series to be added to model
//
void CModel::AddObservedTimeSeries(CTimeSeriesABC *pTS)
{
  if (!DynArrayAppend((void**&)(_pObservedTS),(void*)(pTS),_nObservedTS)){
    ExitGracefully("CModel::AddObservedTimeSeries: adding NULL observation time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds observation weighting time series to model
///
/// \param *pTS [in] (valid) pointer to observation weights time series to be added to model
//
void CModel::AddObservedWeightsTS(CTimeSeriesABC   *pTS)
{
  if (!DynArrayAppend((void**&)(_pObsWeightTS),(void*)(pTS),_nObsWeightTS)){
    ExitGracefully("CModel::AddObservedWeightsTS: adding NULL observed weights time series",BAD_DATA);}
}
//////////////////////////////////////////////////////////////////
/// \brief Adds diagnostic to model
///
/// \param *pDiag [in] (valid) pointer to diagnostic to be added to model
//
void CModel::AddDiagnostic(CDiagnostic       *pDiag)
{
  if (!DynArrayAppend((void**&)(_pDiagnostics),(void*)(pDiag),_nDiagnostics)){
    ExitGracefully("CModel::AddDiagnostic: adding NULL diagnostic",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds model output time to model
///
/// \param timestamp [in] ISO-format timestamp string YYYY-MM-DD hh:mm:ss
/// \note must be called after model starttime & duration are set
//
void CModel::AddModelOutputTime   (const time_struct &tt_out, const optStruct &Options)
{

  //local time
  double t_loc = TimeDifference(Options.julian_start_day,Options.julian_start_year,tt_out.julian_day,tt_out.year);

  ExitGracefullyIf(t_loc<0,
                   "AddModelOutputTime: Cannot have model output time prior to start of simulation",BAD_DATA_WARN);
  if (t_loc>Options.duration){
    WriteWarning("AddModelOutputTime: model output time specified after end of simulation. It will be ignored",Options.noisy);
  }

  //add to end of array
  double *tmpOut=new double [_nOutputTimes+1];
  for (int i=0;i<_nOutputTimes;i++){
    tmpOut[i]=_aOutputTimes[i];
  }
  tmpOut[_nOutputTimes]=t_loc;
  delete [] _aOutputTimes;
  _aOutputTimes=tmpOut;
  _nOutputTimes++;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds a hydrological process to system
/// \details Adds a hydrological process that moves water or energy from state
/// variable (e.g., an array of storage units) i[] to state variables j[]
///
/// \param *pHydroProc [in] (valid) pointer to Hydrological process to be added
//
void CModel::AddProcess(CHydroProcessABC *pHydroProc)
{
  ExitGracefullyIf(pHydroProc==NULL                 ,"CModel AddProcess::NULL process"  ,BAD_DATA);
  for (int q=0;q<pHydroProc->GetNumConnections();q++)
  {
    int i=pHydroProc->GetFromIndices()[q];
    int j=pHydroProc->GetToIndices  ()[q];
    ExitGracefullyIf((i<0) && (i>=_nStateVars),"CModel AddProcess::improper storage index",BAD_DATA);
    ExitGracefullyIf((j<0) && (j>=_nStateVars),"CModel AddProcess::improper storage index",BAD_DATA);
  }
  if (!DynArrayAppend((void**&)(_pProcesses),(void*)(pHydroProc),_nProcesses)){
    ExitGracefully("CModel::AddProcess: adding NULL hydrological process",BAD_DATA);}
}

//////////////////////////////////////////////////////////////////
/// \brief Adds state variables during model construction
/// \details Called during model construction (only from parse of RVI file) to dynamically
/// generate a list of state variables in the model
///
/// \note Soil structure (or any SV with multiple layers) needs to be created in Model Constructor beforehand (is this true for other multilayer SVs??)
///
/// \param *aSV [in] array of state variables types to be added
/// \param *aLev [in] array of state variable levels to be added
/// \param nSV [in] Integer number of SVs to be added (size of aSV[] and aLev[])
//
void  CModel::AddStateVariables( const sv_type *aSV,
                                 const int     *aLev,
                                 const int      nSV)
{
  int i,ii;
  bool found;
  for (ii=0;ii<nSV;ii++)
  {
    found=false;
    for (i=0;i<_nStateVars;i++){
      if ((_aStateVarType [i]==aSV [ii]) &&
          (_aStateVarLayer[i]==aLev[ii])){found=true;}
    }
    //found=(_aStateVarIndices[(int)(aSV[ii])][aLev[ii]]!=DOESNT_EXIST);
    if (!found)
    {
      sv_type *tmpSV=new sv_type[_nStateVars+1];
      int     *tmpLy=NULL;
      tmpLy=new int    [_nStateVars+1];
      ExitGracefullyIf(tmpLy==NULL,"CModel::AddStateVariables",OUT_OF_MEMORY);
      for (i=0;i<_nStateVars;i++)//copy old arrays
      {
        tmpSV[i]=_aStateVarType  [i];
        tmpLy[i]=_aStateVarLayer [i];
      }
      //set values of new array items
      tmpSV[_nStateVars]=aSV[ii];
      tmpLy[_nStateVars]=aLev[ii];
      //add index to state variable lookup table
      ExitGracefullyIf((int)(aSV[ii])>MAX_STATE_VARS,
                       "CModel::AddStateVariables: bad type specified",RUNTIME_ERR);
      ExitGracefullyIf((aLev[ii]<-1) || (aLev[ii]>=MAX_SV_LAYERS),
                       "CModel::AddStateVariables: bad layer index specified",RUNTIME_ERR);
      if (aLev[ii]==DOESNT_EXIST){_aStateVarIndices[(int)(aSV[ii])][0       ]=_nStateVars;}
      else                       {_aStateVarIndices[(int)(aSV[ii])][aLev[ii]]=_nStateVars;}
      //delete & replace old arrays
      delete [] _aStateVarType;  _aStateVarType=NULL;
      delete [] _aStateVarLayer; _aStateVarLayer=NULL;
      _aStateVarType =tmpSV;
      _aStateVarLayer=tmpLy;
      _nStateVars++;
      ExitGracefullyIf(_nStateVars>MAX_STATE_VARS,
                       "CModel::AddStateVariables: exceeded maximum number of state variables in model",RUNTIME_ERR);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Adds aquifer state variables during model construction
/// \details Should be called once during parse of .rvi file when :AquiferLayers command is specified
///
/// \param nLayers [in] number of aquifer layers in model
//
void CModel::AddAquiferStateVars(const int nLayers)
{
  _nAquiferLayers=nLayers;
  sv_type *aSV =new sv_type [_nAquiferLayers];
  int     *aLev=new int       [_nAquiferLayers];
  for (int i=0;i<_nAquiferLayers;i++)
  {
    aSV [i]=GROUNDWATER;
    aLev[i]=i;
  }
  AddStateVariables(aSV,aLev,_nAquiferLayers);
  delete [] aSV;
  delete [] aLev;
}

//////////////////////////////////////////////////////////////////
/// \brief Adds custom output object
///
/// \param *pCO [in] (valid) pointer to Custom output object to be added
//
void CModel::AddCustomOutput(CCustomOutput *pCO)
{
  if (!DynArrayAppend((void**&)(_pCustomOutputs),(void*)(pCO),_nCustomOutputs)){
    ExitGracefully("CModel::AddCustomOutput: adding NULL custom output",BAD_DATA);}
}


/*****************************************************************
   Other Manipulator Functions
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Sets lake storage state variable index based on SV type and layer
///
/// \param *lak_sv [out] SV type
///     \param lev [in] Layer of interest
//
void CModel::SetLakeStorage   (const sv_type      lak_sv, const int lev)
{
  _lake_sv=GetStateVarIndex(lak_sv,lev);
  ExitGracefullyIf(_lake_sv==DOESNT_EXIST,
                   "CModel::SetLakeStorage: non-existent state variable",BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief Sets state variable defined by SV and lev to be an aggregated variable
///
/// \param SV [in] SV type to be aggregated
///     \param lev [in] Layer of SV to be aggregated
/// \param group_name [in] HRU group over which aggregation takes place
//
void CModel::SetAggregatedVariable(const sv_type SV, const int lev, const string group_name)
{
  int i=GetStateVarIndex(SV,lev);
  ExitGracefullyIf(i==DOESNT_EXIST,
                   "CModel::SetAggregatedVariable: non-existent state variable",BAD_DATA);
  for (int kk=0;kk<_nHRUGroups;kk++){
    if (!_pHRUGroups[kk]->GetName().compare(group_name)){_pHRUGroups[kk]->SetAsAggregator(i);}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Sets HRU group for which storage output files are generated
///
/// \param *pOut [in] (assumed valid) pointer to HRU group
//
void CModel::SetOutputGroup(const CHRUGroup *pOut){
  _pOutputGroup=pOut;
}

//////////////////////////////////////////////////////////////////
/// \brief sets number of layers used to simulate snowpack
/// \param nLayers # of snow layers
//
void  CModel::SetNumSnowLayers     (const int          nLayers)
{
  ExitGracefullyIf(nLayers<0,
                   "CModel::SetNumSnowLayers: cannot set negative number of snow layers",BAD_DATA);
  sv_type *aSV=new sv_type [nLayers];
  int     *aLev=new int [nLayers];
  for (int m=0;m<nLayers;m++){
    aSV[m]=SNOW;
    aLev[m]=m;
  }
  AddStateVariables(aSV,aLev,nLayers);
  delete [] aSV;
  delete [] aLev;
}

bool IsContinuousFlowObs(CTimeSeriesABC *pObs,long SBID);

//////////////////////////////////////////////////////////////////
/// \brief overrides streamflow with observed streamflow
/// \param SBID [in] valid subbasin identifier of basin with observations at outflow
//
void CModel::OverrideStreamflow   (const long SBID)
{
  for (int i=0;i<_nObservedTS; i++)
  {
    if (IsContinuousFlowObs(_pObservedTS[i],SBID))
    {
      //check for blanks in observation TS
      bool bad=false;
      for (int n=0;n<_pObservedTS[i]->GetNumValues();n++){
        if (_pObservedTS[i]->GetValue(n)==RAV_BLANK_DATA){bad=true;break;}
      }
      if (bad){
        WriteWarning("CModel::OverrideStreamflow::cannot override streamflow if there are blanks in observation data",_pOptStruct->noisy);
        return;
      }

      long downID=GetSubBasinByID(SBID)->GetDownstreamID();
      if (downID!=DOESNT_EXIST)
      {
        //Copy time series of observed flows to new time series
        string name="Inflow_Hydrograph_"+to_string(SBID);
        CTimeSeries *pObs=dynamic_cast<CTimeSeries *>(_pObservedTS[i]);
        CTimeSeries *pTS =new CTimeSeries(name,*pObs);//copy time series
        pTS->SetTag(to_string(SBID));

        //add as inflow hydrograph to downstream
        GetSubBasinByID(downID)->AddInflowHydrograph(pTS);
        GetSubBasinByID(SBID)->SetDownstreamID(DOESNT_EXIST);
        return;
      }
      else{
        WriteWarning("CModel::OverrideStreamflow: overriding streamflow at an outlet subbasin has no impact on model operation",_pOptStruct->noisy);
      }
    }
  }
}

/*****************************************************************
   Routines called repeatedly during model simulation
------------------------------------------------------------------
*****************************************************************/

//////////////////////////////////////////////////////////////////
/// \brief Increments water/mass/energy balance
/// \details add cumulative amount of water/mass/energy for each process connection
/// \note  called from main solver
///
/// \param q_star [in] Integer index of hydrologic process connection
/// \param k [in] Integer index of HRU
/// \param moved [in] Double amount balance is to be incremented (mm or mg/m2 or MJ/m2)
//
void CModel::IncrementBalance(const int    q_star,
                              const int    k,
                              const double moved){
#ifdef _STRICTCHECK_
  ExitGracefullyIf(q_star>_nTotalConnections,"CModel::IncrementBalance: bad index",RUNTIME_ERR);
#endif
  _aCumulativeBal[k][q_star]+=moved;
  _aFlowBal      [k][q_star]= moved;
}
//////////////////////////////////////////////////////////////////
/// \brief Increments lateral flow water/energy balance
/// \details add cumulative amount of water/energy for each lateral process connection
/// \note  called from main solver
///
/// \param jss [in] Integer index of hydrologic lateral process connection
/// \param moved [in] Double amount balance is to be incremented [mm-m2] or [MJ] or [mg]
//
void CModel::IncrementLatBalance( const int jss,
                                  const double moved)//[mm-m2] or [MJ]
{
  _aCumulativeLatBal[jss]+=moved;
  _aFlowLatBal      [jss]=moved; 
}
//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative mass & energy added to system (precipitation/ustream basin flows, etc.)
/// \details Increment cumulative precipitation based on average preciptation and corresponding timestep [mm]
///
/// \param &Options [in] Global model options information
/// \param &tt [in] current time
//
void CModel::IncrementCumulInput(const optStruct &Options, const time_struct &tt)
{
  double area;
  _CumulInput+=GetAveragePrecip()*Options.timestep;

  for (int p=0;p<_nSubBasins;p++){
    area = _WatershedArea*M2_PER_KM2;
    _CumulInput+=_pSubBasins[p]->GetIntegratedSpecInflow(tt.model_time,Options.timestep)/area*MM_PER_METER;//converted to [mm] over  basin
  }

  _pTransModel->IncrementCumulInput(Options,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief Increments cumulative outflow from system for mass balance diagnostics
/// \details Increment cumulative outflow according to timestep, flow, and area of basin
///
/// \param &Options [in] Global model options information
//
void CModel::IncrementCumOutflow(const optStruct &Options)
{
  double area=(_WatershedArea*M2_PER_KM2);
  for(int p=0;p<_nSubBasins;p++)
  {
    bool outflowdisabled;
    CSubBasin *pSBdown=NULL;
    pSBdown=GetSubBasinByID(_pSubBasins[p]->GetDownstreamID());
    outflowdisabled=((pSBdown==NULL) || (!pSBdown->IsEnabled()));

    if(_pSubBasins[p]->IsEnabled()){
      if((_aSubBasinOrder[p]==0) || (outflowdisabled))//outlet does not drain into another subbasin
      {
        _CumulOutput+=_pSubBasins[p]->GetIntegratedOutflow(Options.timestep)/area*MM_PER_METER;//converted to [mm] over entire watershed
      }
      _CumulOutput+=_pSubBasins[p]->GetReservoirLosses(Options.timestep)/area*MM_PER_METER;
    }
  }
  _pTransModel->IncrementCumulOutput(Options);
}

//////////////////////////////////////////////////////////////////
/// \brief Updates values of user-specified transient parameters, updates changes to land use class
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateTransientParams(const optStruct   &Options,
                                   const time_struct &tt)
{
  int nn=(int)((tt.model_time+REAL_SMALL)/Options.timestep);//current timestep index
  for (int j=0;j<_nTransParams;j++)
  {
    class_type ctype=_pTransParams[j]->GetParameterClassType();
    string     pname=_pTransParams[j]->GetParameterName();
    string     cname=_pTransParams[j]->GetParameterClass();
    double     value=_pTransParams[j]->GetTimeSeries()->GetSampledValue(nn);

    if (ctype==CLASS_SOIL)
    {
      CSoilClass::StringToSoilClass(cname)->SetSoilProperty(pname,value);
    }
    else if (ctype==CLASS_VEGETATION)
    {
      CVegetationClass::StringToVegClass(cname)->SetVegetationProperty(pname,value);
    }
    else if (ctype==CLASS_TERRAIN)
    {
      CTerrainClass::StringToTerrainClass(cname)->SetTerrainProperty(pname,value);
    }
    else if (ctype==CLASS_LANDUSE)
    {
      CLandUseClass::StringToLUClass(cname)->SetSurfaceProperty(pname,value);
    }
    else if (ctype==CLASS_GLOBAL)
    {
      CGlobalParams::SetGlobalProperty(pname,value);
    }
  }
  int k;
  for (int j = 0; j<_nClassChanges; j++)
  {
    if ((_pClassChanges[j]->modeltime > tt.model_time - TIME_CORRECTION) &&
        (_pClassChanges[j]->modeltime < tt.model_time + Options.timestep))
    {//change happens this time step

      int kk = _pClassChanges[j]->HRU_groupID;
      if (_pClassChanges[j]->tclass == CLASS_LANDUSE){
        //cout << "LAND USE CHANGE!"<<endl;
        for (int k_loc = 0; k_loc < _pHRUGroups[kk]->GetNumHRUs();k_loc++)
        {
          k=_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex();
          CLandUseClass *lult_class= CLandUseClass::StringToLUClass(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeLandUse(lult_class);
        }
      }
      else if (_pClassChanges[j]->tclass == CLASS_VEGETATION){
        //cout << "VEGETATION CHANGE! "<< _pClassChanges[j]->modeltime << " "<<tt.model_time<<endl;
        for (int k_loc = 0; k_loc < _pHRUGroups[kk]->GetNumHRUs();k_loc++)
        {
          k=_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex();
          CVegetationClass *veg_class= CVegetationClass::StringToVegClass(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeVegetation(veg_class);
        }
      }
      else if (_pClassChanges[j]->tclass == CLASS_HRUTYPE){
        //cout << "HRU TYPE CHANGE! "<< _pClassChanges[j]->modeltime << " "<<tt.model_time<<endl;
        for (int k_loc = 0; k_loc < _pHRUGroups[kk]->GetNumHRUs();k_loc++)
        {
          k=_pHRUGroups[kk]->GetHRU(k_loc)->GetGlobalIndex();
          HRU_type typ;
          typ=StringToHRUType(_pClassChanges[j]->newclass);
          _pHydroUnits[k]->ChangeHRUType(typ);
        }
      }
    }
  }

}
//////////////////////////////////////////////////////////////////
/// \brief Recalculates HRU derived parameters
/// \details Recalculate HRU derived parameters that are based upon time-of-year/day and SVs (storage, temp)
/// \remark tt is model time, not global time
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::RecalculateHRUDerivedParams(const optStruct    &Options,
                                         const time_struct  &tt)
{
  for (int k=0;k<_nHydroUnits;k++)
  {
    if(_pHydroUnits[k]->IsEnabled())
    {
      _pHydroUnits[k]->RecalculateDerivedParams(Options,tt);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Updates values stored in modeled time series of observation data
///
/// \param &Options [in] Global model options information
/// \param &tt [in] Current time structure
//
void CModel::UpdateDiagnostics(const optStruct   &Options,
                               const time_struct &tt)
{
  if (_nDiagnostics==0){return;}

  int n=(int)(floor((tt.model_time+TIME_CORRECTION)/Options.timestep));//current timestep index
  double value, obsTime;
  int layer_ind;
  CSubBasin *pBasin=NULL;  
  for (int i=0;i<_nObservedTS;i++)
  {
    string datatype=_pObservedTS[i]->GetName();
    sv_type svtyp=CStateVariable::StringToSVType(datatype,layer_ind,false);

    if (datatype=="HYDROGRAPH")
    {
      pBasin=GetSubBasinByID (s_to_l(_pObservedTS[i]->GetTag().c_str()));

      if ((Options.ave_hydrograph) && (tt.model_time!=0))
      {
        value=pBasin->GetIntegratedOutflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
      }
      else
      {
        value=pBasin->GetOutflowRate();
      }
    }
    else if (datatype == "RESERVOIR_STAGE")
    {
      pBasin=GetSubBasinByID (s_to_l(_pObservedTS[i]->GetTag().c_str()));
      value = pBasin->GetReservoir()->GetResStage();
    }
    else if (datatype == "RESERVOIR_INFLOW")
    {
      pBasin = GetSubBasinByID(s_to_l(_pObservedTS[i]->GetTag().c_str()));
      value = pBasin->GetIntegratedReservoirInflow(Options.timestep)/(Options.timestep*SEC_PER_DAY);
    }
		else if (datatype == "RESERVOIR_NETINFLOW")
		{
			pBasin = GetSubBasinByID(s_to_l(_pObservedTS[i]->GetTag().c_str()));
			CReservoir *pRes= pBasin->GetReservoir(); 
			CHydroUnit *pHRU= GetHRUByID(pRes->GetHRUIndex());
			double avg_area = pHRU->GetArea();

			if (pHRU == NULL) {avg_area = 0.0;}

			double tem_precip1 = pBasin->GetAvgForcing("PRECIP")*avg_area*M2_PER_KM2/MM_PER_METER / (Options.timestep*SEC_PER_DAY); 
			double losses = pBasin->GetReservoirEvapLosses(Options.timestep) / (Options.timestep*SEC_PER_DAY);
			value = pBasin->GetIntegratedReservoirInflow(Options.timestep) / (Options.timestep*SEC_PER_DAY) + tem_precip1 - losses;
		}
    else if (svtyp!=UNRECOGNIZED_SVTYPE)
    {
      CHydroUnit *pHRU=NULL;
      pHRU=GetHRUByID(s_to_i(_pObservedTS[i]->GetTag().c_str()));
      string error="CModel::UpdateDiagnostics: Invalid HRU ID specified in observed state variable time series "+_pObservedTS[i]->GetName();
      ExitGracefullyIf(pHRU==NULL,error.c_str(),BAD_DATA);
      int sv_index=this->GetStateVarIndex(svtyp,layer_ind);
      value=pHRU->GetStateVarValue(sv_index);
    }
    else
    {
      if (tt.model_time ==0){
        WriteWarning("CModel::UpdateDiagnostics: invalid tag used for specifying Observation type",BAD_DATA);
      }
      value=0;
    }
    _pModeledTS[i]->SetValue (n,value);


    obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]); // time of the next observation

		/* K. Lee - The +Options.timestep is breaking irregular observation diagnostics.
		Do not fully understand what Nick intended to fix by adding +Options.timestep, but taking it out fixes the issue 
		and does not break regular observations when tested with the tutorial files.
		Regular observations filled with blanks and irregular observations yield same result when taken out
	
	    while((tt.model_time+Options.timestep >= obsTime+_pObservedTS[i]->GetSampledInterval()) &&  //N .Sgro Fix
	          (_aObsIndex[i]<_pObservedTS[i]->GetNumSampledValues()))
	    */
		while ((tt.model_time >= obsTime + _pObservedTS[i]->GetSampledInterval()) &&  
			(_aObsIndex[i]<_pObservedTS[i]->GetNumSampledValues()))
		{
      value=RAV_BLANK_DATA;
      // only set values within diagnostic evaluation times. The rest stay as BLANK_DATA
      if ((obsTime >= Options.diag_start_time) && (obsTime <= Options.diag_end_time))
      {
        value= _pModeledTS[i]->GetModelledValue(obsTime,_pObservedTS[i]->GetType());
      }

      _pModeledTS[i]->SetSampledValue(_aObsIndex[i],value);
      _aObsIndex[i]++;
      obsTime =_pObservedTS[i]->GetSampledTime(_aObsIndex[i]);
    }
  }
}
//////////////////////////////////////////////////////////////////
/// \brief Apply hydrological process to model
/// \details Method returns rate of mass/energy transfers rates_of_change [mm/d, mg/m2/d, or MJ/m2/d] from a set
/// of state variables iFrom[] to a set of state variables iTo[] (e.g., expected water movement [mm/d]
// from storage unit iFrom[0] to storage unit iTo[0]) in the given HRU using state_var[] as the expected
/// value of all state variables over the timestep.
///
/// \param j [in] Integer process indentifier
/// \param *state_var [in] Array of state variables for HRU
/// \param *pHRU [in] Pointer to HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Time structure
/// \param *iFrom [in] Array (size: nConnections)  of Indices of state variable losing mass or energy
/// \param *iTo [in] Array (size: nConnections)  of indices of state variable gaining mass or energy
/// \param &nConnections [out] Number of connections between storage units/state vars
/// \param *rates_of_change [out] Double array (size: nConnections) of loss/gain rates of water [mm/d], mass [mg/m2/d], and/or energy [MJ/m2/d]
/// \return returns false if this process doesn't apply to this HRU, true otherwise
//
bool CModel::ApplyProcess ( const int          j,                    //process identifier
                            const double      *state_var,            //array of state variables for HRU
                            const CHydroUnit  *pHRU,                 //pointer to HRU
                            const optStruct   &Options,
                            const time_struct &tt,
                                  int         *iFrom,                //indices of state variable losing water or heat
                                  int         *iTo,                  //indices of state variable gaining water or heat
                                  int         &nConnections,         //number of connections between storage units/state vars
                                  double      *rates_of_change) const//loss/gain rates of water [mm/d] and energy [MJ/m2/d]
{
#ifdef _STRICTCHECK_
  ExitGracefullyIf((j<0) && (j>=_nProcesses),"CModel ApplyProcess::improper index",BAD_DATA);
#endif
  CHydroProcessABC *pProc=_pProcesses[j];

  nConnections=pProc->GetNumConnections();//total connections: nConnections+nCascades
  if (!_aShouldApplyProcess[j][pHRU->GetGlobalIndex()]){return false;}
  if (!pHRU->IsEnabled()){return false;}

  for (int q=0;q<nConnections;q++)
  {
    iFrom[q]=pProc->GetFromIndices()[q];
    iTo  [q]=pProc->GetToIndices  ()[q];
    rates_of_change[q]=0.0;
  }

  pProc->GetRatesOfChange(state_var,pHRU,Options,tt,rates_of_change);

  //special cascade handling (prior to applying constraints)
  //------------------------------------------------------------------------
  if (pProc->HasCascade())
  {
    int nCascades;
    static double max_state_var[MAX_STATE_VARS];
    nCascades=pProc->GetNumCascades();
    for (int i=0;i<_nStateVars;i++){//should only calculate for participating compartments
      max_state_var[i]=pHRU->GetStateVarMax(i,state_var,Options);
    }
    for (int q=0;q<nCascades;q++)
    {
      iFrom[nConnections-nCascades+q]=pProc->GetCascadeFromIndex();
      iTo  [nConnections-nCascades+q]=pProc->GetCascadeToIndices()[q];
    }
    pProc->Cascade(rates_of_change,state_var,&max_state_var[0],Options.timestep);
  }

  pProc->ApplyConstraints(state_var,pHRU,Options,tt,rates_of_change);
  return true;
}
//////////////////////////////////////////////////////////////////
/// \brief Apply lateral exchange hydrological process to model
/// \details Method returns rate of mass/energy transfers rates_of_change [mm/d, mg/m2/d, or MJ/m2/d] from a set
/// of state variables iFrom[] in HRUs kFrom[] to a set of state variables iTo[] in HRUs kTo[] (e.g., expected water movement [mm/d]
// from storage unit iFrom[0] in HRU kFrom[0] to storage unit iTo[0] in HRU kFrom[0]) in the given HRU using state_vars[][] as the expected
/// value of all state variables over the timestep.
///
/// \param j [in] Integer process indentifier
/// \param **state_vars [in] Array of state variables for all HRUs (size: [nHRUs][nStateVars])
/// \param *pHRU [in] Pointer to HRU
/// \param &Options [in] Global model options information
/// \param &tt [in] Time structure
/// \param *kFrom [in] Array (size: nLatConnections)  of Indices of HRU losing mass or energy
/// \param *kTo [in] Array (size: nLatConnections)  of indices of HRU gaining mass or energy
/// \param *iFrom [in] Array (size: nLatConnections)  of Indices of state variable losing mass or energy
/// \param *iTo [in] Array (size: nLatConnections)  of indices of state variable gaining mass or energy
/// \param &nConnections [out] Number of connections between storage units/state vars
/// \param *exchange_rates [out] Double array (size: nConnections) of loss/gain rates of water [mm-m2/d], mass [mg/d], and/or energy [MJ/d] 
/// \return returns false if this process doesn't apply to this HRU, true otherwise
//
bool CModel::ApplyLateralProcess( const int          j,
                                  const double* const* state_vars,
                                  const optStruct   &Options,
                                  const time_struct &tt,
                                        int         *kFrom,
                                        int         *kTo,
                                        int         *iFrom,
                                        int         *iTo,
                                        int         &nLatConnections,
                                        double      *exchange_rates) const
{
  CLateralExchangeProcessABC *pLatProc;

  ExitGracefullyIf((j<0) && (j>=_nProcesses),"CModel ApplyProcess::improper index",BAD_DATA);

  nLatConnections=_pProcesses[j]->GetNumLatConnections();

  if(nLatConnections==0){return false;}

  pLatProc=(CLateralExchangeProcessABC*)_pProcesses[j]; // Cast
  if (!_aShouldApplyProcess[j][pLatProc->GetFromHRUIndices()[0]]){return false;} //JRC: is the From/0 appropriate?
  
  
  for (int q=0;q<nLatConnections;q++)
  {
    iFrom[q]=pLatProc->GetLateralFromIndices()[q];
    iTo  [q]=pLatProc->GetLateralToIndices()[q];
    kFrom[q]=pLatProc->GetFromHRUIndices()[q];
    kTo  [q]=pLatProc->GetToHRUIndices()[q];
    exchange_rates[q]=0.0;
  }
  for(int q=0;q<nLatConnections;q++)
  {
    if(!_pHydroUnits[kFrom[q]]->IsEnabled()) { return false; }
    if(!_pHydroUnits[kTo  [q]]->IsEnabled()) { return false; } //ALL participating HRUs must be enabled to apply
  }
  pLatProc->GetLateralExchange(state_vars,_pHydroUnits,Options,tt,exchange_rates);

  return true;
}


