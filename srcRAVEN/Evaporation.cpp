/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team

  Routines for calculating PET:
  -Penman Monteith Equation
  -Penman-Combination
  -Priestly Taylor Equation
  -Hargreaves Equation
  -Hargreaves Equation
  =Shuttle-Wallace
  -Granger
  ----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "Properties.h"
#include "Radiation.h"
#include "HydroUnits.h"
#include "Model.h"

double EvapGranger(const force_struct *F,const CHydroUnit *pHRU);

//////////////////////////////////////////////////////////////////
/// \brief Calculates PET using Makink 1957 method \cite Lu2005JotAWRA
/// \ref From Makkink, 1957 as described in "A comparison of six potential
/// evapotranspiration methods for regional use in the southeastern
/// united states", American Water Resources Association, 2005. \cite Makkink1957JIWE
/// \note This is a utility function called by EstimatePET
/// \remark Added by Graham Stonebridge, Fall 2011
/// \param *F [in] Reference to model forcing functions
/// \return Calculated PET [mm/d]
//
double Makkink1957Evap(const force_struct *F)
{
  double PET;
  double gamma;     //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]

  sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
  de_dT    =GetSatVapSlope           (F->temp_ave,sat_vap);
  LH_vapor =GetLatentHeatVaporization(F->temp_ave);
  gamma    =GetPsychometricConstant  (F->air_pres,LH_vapor);

  PET=0.61*(de_dT/(de_dT+gamma))*F->SW_radia*23.8846/58.5-0.12;

  return max(0.0,PET);
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates PET using Turc 1961 method \cite Lu2005JotAWRA
/// \ref From Turc, L., 1961. Evaluation de besoins en eau d�irrigation, ET potentielle, Ann. Agron. 12:13-49. \cite turc1961AA
/// as defined in "A comparison of six potential 
/// evapotranspiration methods for regional use in the southeastern
/// united states", American Water Resources Association, 2005.
/// \note This is a utility function called by EstimatePET
/// \ref Added by Graham Stonebridge, Fall 2011
/// \param *F [in] Reference to model forcing functions
/// \return Calculated PET [mm/d]
//
double TurcEvap(const force_struct *F)
{
  double t,evp;
  t=max(0.0,F->temp_daily_ave);
  evp=0.013*t/(t+15)*((F->SW_radia)*23.8846+50);
  if (F->rel_humidity<0.5){evp*=(1+(50-F->rel_humidity*100)/70);}
  return max(0.0,evp);
}

//////////////////////////////////////////////////////////////////
/// \brief Returns evaporation rate [mm/d]
/// \details Returns evaporation rate using the Penman-monteith equation, which
/// is an energy-balance approach
/// \ref  Adapted from Dingman, Howell, T.A and Evett, S.R. (USDA-Agricultural research service) \cite Howell2004 \cite Dingman2004
/// http://www.cprl.ars.usda.gov/pdfs/PM COLO Bar 2004 corrected 9apr04.pdf
/// \note This is a utility function called by EstimatePET
///
/// \param *F [in] Forcing functions for a specific HRU over the current time step
/// \param &atmos_cond [in] Atmospheric conductance [mm/s]
/// \param &canopy_cond [in] Canopy conductance [mm/s]
/// \return Evaporation rate [mm/d]
//
double PenmanMonteithEvap(const force_struct     *F,
                          const double       &atmos_cond,   //[mm/s]
                          const double       &canopy_cond)  //[mm/s]
{
  double numer, denom;
  double gamma;     //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de* / dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]
  double vapor_def; //vapor deficit [kPa]

  sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
  de_dT    =GetSatVapSlope           (F->temp_ave,sat_vap);
  LH_vapor =GetLatentHeatVaporization(F->temp_ave);
  gamma    =GetPsychometricConstant  (F->air_pres,LH_vapor);
  vapor_def=sat_vap*(1.0-F->rel_humidity);

  //Calculate evaporation - Dingman eqn 7-56
  numer =de_dT*max(F->SW_radia_net + F->LW_radia,0.0); //[kPa/K*MJ/m2/d]
  numer+=(F->air_dens)*SPH_AIR*vapor_def*(atmos_cond*SEC_PER_DAY/MM_PER_METER);//[kPa/K*MJ/m2/d]
  denom = (de_dT+gamma*(1.0+(atmos_cond/canopy_cond)))*LH_vapor*DENSITY_WATER; //[kPa/K*MJ/m3]
  if (canopy_cond==0){return 0.0;}//zero conductance means no ET

  return numer/denom*MM_PER_METER;//[mm/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Returns potential evaporation rate [mm/d] \cite Dingman2004
/// \details returns the potential evaporation rate [mm/d] using the Penman-Combination
//  equation which is an energy and mass balance approach
/// \remark Adapted from Dingman pg 285-6
/// \note This is a utility function called by EstimatePET
/// \param *F [in] Forcing functions for a specific HRU over the current time step
/// \param &vert_trans Vertical transmissivity [m2/kg]
/// \return Potential evaporation rate [mm/d]
//
double PenmanCombinationEvap(const force_struct *F,
                             const double &vert_trans)    //[m*s^2/kg]
{
  double gamma;     //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]
  double vapor_def; //vapor deficit [kPa]

  double numer,denom;

  sat_vap               =GetSaturatedVaporPressure(F->temp_ave);
  de_dT                 =GetSatVapSlope           (F->temp_ave,sat_vap);
  LH_vapor              =GetLatentHeatVaporization(F->temp_ave);
  gamma                 =GetPsychometricConstant  (F->air_pres,LH_vapor);

  vapor_def   =sat_vap*(1.0-(F->rel_humidity));

  /// \ref Calculate evaporation - Dingman eqn 7-33
  numer =de_dT*max((F->SW_radia_net) + (F->LW_radia),0.0);
  numer+=gamma*vert_trans*DENSITY_WATER*LH_vapor*(F->wind_vel)*vapor_def;

  denom =DENSITY_WATER*LH_vapor*(de_dT + gamma);

  return max(0.0,(numer/denom))*MM_PER_METER;//[mm/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Returns potential evaporation rate [mm/d] \cite Stannard1993WRR
/// \details  returns the potential evaporation rate [mm/d] using the
/// Priestley-Taylor equation
/// \ref adapted from SWAT
/// \note This is a utility function called by EstimatePET
/// \param *F [in] Forcing functions for a specific HRU over the current time step
/// \return Potential evaporation rate [mm/d]
//
double PriestleyTaylorEvap(const force_struct *F)
{
  double gamma;     //psychometric "constant" [kPa/K]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double sat_vap;   //Saturation vapor pressure

  sat_vap =GetSaturatedVaporPressure(F->temp_ave);
  de_dT   =GetSatVapSlope                                               (F->temp_ave,sat_vap);
  LH_vapor=GetLatentHeatVaporization(F->temp_ave);
  gamma   =GetPsychometricConstant  (F->air_pres,LH_vapor);

  return 1.28 * (de_dT/(de_dT+gamma))*max(F->SW_radia_net+F->LW_radia,0.0)/LH_vapor/DENSITY_WATER*MM_PER_METER;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns potential evaporation rate [mm/d]
/// \details  returns the potential evaporation rate [mm/d] using the
/// Hargreaves equation
/// \ref adapted from WATFLOOD Manual \cite Kouwen2011
/// \note This is a utility function called by EstimatePET
/// \param *F [in] Forcing functions for a specific HRU over the current time step
/// \return Potential evaporation rate [mm/d]
//
double HargreavesEvap(const force_struct *F)//[C]
{
  const double HARGREAVES_CONST=0.0075;
  double Ra;
  double rel_hum,delT; //relative humidity[0..1], monthly temperature range [F]
  double Ct;

  Ra=F->ET_radia;//[MJ/m2/d]
  Ra/=(LH_VAPOR*DENSITY_WATER/MM_PER_METER);//[MJ/m2/day]->[mm/d]

  rel_hum=F->rel_humidity;

  Ct=0.125;
  if (rel_hum>=0.54){Ct=0.035*pow(100*(1.0-rel_hum),0.333);}

  delT=CelsiusToFarenheit(F->temp_month_max)-
       CelsiusToFarenheit(F->temp_month_min);

  // \todo [optimize] - move this check to initialization routine
  ExitGracefullyIf(F->temp_month_max==NOT_SPECIFIED,
                   "PET_HARGREAVES requires minimum and maximum monthly temperatures",BAD_DATA);

  return max(0.0,HARGREAVES_CONST*Ra*Ct*sqrt(delT)*CelsiusToFarenheit(F->temp_ave));
}

//////////////////////////////////////////////////////////////////
/// \brief Returns potential evaporation rate [mm/d]
/// \details  returns the potential evaporation rate [mm/d] using the Hargreaves 1985 equation
/// \remark adapted from History and Evaluation of Hargreaves Evapotranspiration
/// Equation, GH Hargreaves, RG Allen, J. of Irrigation and Drainage Eng. \cite Allen2003JoIaDE
/// Jan/Feb 2003
/// \note This is a utility function called by EstimatePET
/// \param *F [in] Forcing functions for a specific HRU over the current time step
/// \return Potential evaporation rate [mm/d]
//
double Hargreaves1985Evap(const force_struct *F)//[C]
{
  const double HARGREAVES_CONST=0.0023;
  double Ra;   //maximum incoming solar radiation
  double delT; //range in monthly temperatures

  Ra=F->ET_radia;//[MJ/m2/day]
  Ra/=(LH_VAPOR*DENSITY_WATER/MM_PER_METER);//[MJ/m2/day]->[mm/d]

  delT=F->temp_daily_max-F->temp_daily_min;
  delT=max(delT,0.0);

  return max(0.0,HARGREAVES_CONST*Ra*sqrt(delT)*(F->temp_ave+17.8));
}
/*****************************************************************
Jensen Haise Evaporation
******************************************************************
/// returns the potential evaporation rate (mm/d) using the Jensen &
/// Haise (1963) model as outlined in the PRMS Manual \cite jensen1963JoIDD
----------------------------------------------------------------*/
/*double JensenHaise1963Evap(const force_struct *F,
  const double       &sat_vap_max,// monthly ave. saturated vapor pressure in hottest summer month[KPa]
  const double       &sat_vap_min,// monthly ave. saturated vapor pressure in coldest winter month[KPa],//min summer monthly temp [C]
  const double       &elev,//elevation [masl]
  const double       &julian_day)
  {
  double Ta,ctx,Rin,ch,c1,cts;

  c1=68.0-(3.6*FEET_PER_METER*elev*0.001);
  ch=50/(sat_vap_max-sat_vap_min);
  cts=1.0/(c1+13.0*ch);

  Ta =CelsiusToFarenheit(F->temp_daily_ave);

  ctx=27.5-0.25*(sat_vap_max-sat_vap_min)*MB_PER_KPA-((FEET_PER_METER*elev)/1000);

  Rin=(F->LW_radia+F->SW_radia_net);
  Rin/=(LH_VAPOR*DENSITY_WATER/MM_PER_METER);//[MJ/m2/day]->[mm/d]

  return cts*(Ta-ctx)*Rin; */

//////////////////////////////////////////////////////////////////
/// \brief Returns PET [mm/d] using Hamon 1961 method
/// \ref Hamon (1961) model as outlined in the PRMS Manual \cite Hamon1961
/// \note This is a utility function called by EstimatePET
/// \param *F [in] Model forcing functions
/// \return Calculated PET [mm/d]
//
double Hamon1961Evap(const force_struct *F)
{
  double abs_hum,sat_vap;
  sat_vap=GetSaturatedVaporPressure(F->temp_daily_ave);//KPa

  abs_hum=216.7*(sat_vap*MB_PER_KPA)/(F->temp_daily_ave+ZERO_CELSIUS);//abs. humidity, g/m3 (may wish to make separate function of T)s

  return 0.0055*4.0*abs_hum*F->day_length*F->day_length*MM_PER_INCH;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates PET using known canopy properties and current forcing functions over time step
/// \param *F [in] Model forcing functions
/// \param *pHRU [in] Reference to the specified HRU
/// \param evap_type [in] Method of evaporation calculation selected
/// \param ref_elevation [in] reference elevation for forcing function estimation
/// \param &tt [in] Current model time
/// \return Calculated PET [mm/d]
//
double CModel::EstimatePET(const force_struct &F,
                           const CHydroUnit   *pHRU,
                           const double       &wind_measurement_ht,
                           const double       &ref_elevation,
                           const evap_method   evap_type ,
                           const optStruct    &Options,
                           const time_struct  &tt)
{
  double PET=0.0;

  switch(evap_type)
  {
  case(PET_CONSTANT):
  {
    PET =3.0; break;
  }
  case(PET_DATA):
  {
    PET=F.PET; break;//calculated direct from Gauge
  }
  case(PET_FROMMONTHLY):
  {
    double peRatio=1.0+HBV_PET_TEMP_CORR*(F.temp_ave_unc-F.temp_month_ave);
    peRatio=max(0.0,min(2.0,peRatio));
    PET=F.PET_month_ave*peRatio;
    break;
  }
  case(PET_MONTHLY_FACTOR):
  {
    double forest_corr,max_temp;
    double PET_corr=pHRU->GetSurfaceProps()->forest_PET_corr; //(A0PEFO in UBCWM)
    double Fc      =pHRU->GetSurfaceProps()->forest_coverage;

    //max_temp=F.temp_daily_max;//MQ UBCWM
    max_temp=F.temp_max_unc;
    forest_corr=((PET_corr)*Fc + 1.0*(1.0-Fc)); //XEVAP3

    //orographic corrections - necessary evil for UBCWM emulation; can't separate
    const double A0PELA=0.9;
    double refelev=ref_elevation;
    if(Options.keepUBCWMbugs){refelev=g_debug_vars[4];}
    double XEVAP2=A0PELA * 0.001 * (refelev - pHRU->GetElevation());
    PET=forest_corr*(F.PET_month_ave*max_temp+XEVAP2); //PET_month_ave actually stores Monthly evap factors [mm/d/K]

    PET=max(PET,0.0); 

    break;
  }
  case(PET_PENMAN_MONTEITH):
  {
    double can_cond;
    double zero_pl,rough;
    double vap_rough_ht,atmos_cond;
    double ref_ht;

    can_cond=pHRU->GetVegVarProps()->canopy_conductance;
    zero_pl =pHRU->GetVegVarProps()->zero_pln_disp;
    rough   =pHRU->GetVegVarProps()->roughness;
    ref_ht  =pHRU->GetVegVarProps()->reference_height; //default veght+2.0 m

    if(wind_measurement_ht>ref_ht){ref_ht=wind_measurement_ht;} //correction if real measurement height data is available

    vap_rough_ht=0.1*rough;//=1.0*rough for dingman

    atmos_cond=CalcAtmosphericConductance(F.wind_vel,ref_ht,zero_pl,rough,vap_rough_ht);

    PET=PenmanMonteithEvap(&F,atmos_cond,can_cond); break;
  }
  case(PET_PENMAN_COMBINATION):
  {
    double zero_pl,z0,vert_trans,ref_ht;

    zero_pl =pHRU->GetVegVarProps()->zero_pln_disp;
    z0   =pHRU->GetVegVarProps()->roughness;
    ref_ht  =pHRU->GetVegVarProps()->reference_height; //default veght+2.0 m

    if(wind_measurement_ht>ref_ht){ref_ht=wind_measurement_ht;} //correction if real measurement height data is available

    vert_trans =GetVerticalTransportEfficiency(F.air_pres,ref_ht,zero_pl,z0);

    PET   =PenmanCombinationEvap(&F,vert_trans); break;
  }
  case(PET_JENSEN_HAISE):
  {
    //double sat_vap_max
    ExitGracefully("PET_JENSEN_HAISE",STUB);
    PET=0.0;break;//JensenHaise1963Evap(&F,sat_vap_max,sat_vap_min,elev,julian_day);
  }
  case(PET_HAMON):
  {
    PET=Hamon1961Evap(&F); break;
  }
  case(PET_PRIESTLEY_TAYLOR):
  {
    PET=PriestleyTaylorEvap(&F); break;
  }
  case(PET_HARGREAVES):
  {
    PET=HargreavesEvap(&F); break;
  }
  case(PET_HARGREAVES_1985):
  {
    PET=Hargreaves1985Evap(&F); break;
  }
  case(PET_TURC_1961):
  {
    PET=TurcEvap(&F); break;
  }
  case(PET_MAKKINK_1957):
  {
    PET=Makkink1957Evap(&F);break;
  }
  case(PET_SHUTTLEWORTH_WALLACE):
  {
    //PET=ShuttleworthWallaceEvap(&F,matric_pot,&S,&G,&CV);  // \todo [funct] (need to import additional data)
    ExitGracefully("EstimatePET:Shuttleworth Wallace",STUB);
  }
  case(PET_PENMAN_SIMPLE33) :
  {
    double Rs =F.SW_radia;   //[MJ/m2/d]
    double R_et =F.ET_radia; //[MJ/m2/d]
    double Tave=F.temp_ave;  //[C]
    double rel_hum=F.rel_humidity; //[0..1]

    PET = 0.047*Rs*sqrt(Tave + 9.5) - 2.4*pow(Rs / R_et, 2.0) + 0.09*(Tave + 20)*(1-rel_hum);
    break;
  }
  case(PET_PENMAN_SIMPLE39) :
  {
    double Rs =F.SW_radia;   //[MJ/m2/d]
    double R_et =F.ET_radia; //[MJ/m2/d]
    double Tave=F.temp_ave;  //[C]
    double rel_hum=F.rel_humidity; //[0..1]

    PET = 0.038*Rs*sqrt(Tave + 9.5) - 2.4*pow(Rs / R_et, 2.0) + 0.075*(Tave + 20)*(1-rel_hum);
    break;
  }
  case(PET_MOHYSE) :
  {
    double lat_rad=pHRU->GetLatRad();
    double declin=CRadiation::SolarDeclination(F.day_angle);
    double cpet=CGlobalParams::GetParams()->MOHYSE_PET_coeff;

    PET = cpet/PI*acos(-tan(lat_rad)*tan(declin))*exp((17.3*F.temp_ave)/(238+F.temp_ave));
    break;
  }
  case(PET_OUDIN) :
  {
    PET=max(F.ET_radia/DENSITY_WATER/LH_VAPOR*(F.temp_daily_ave+5.0)/100.0,0.0);
    break;
  }
  case (PET_GRANGER):
    PET=EvapGranger(&F,pHRU);
    break;
  default:
  {
    ExitGracefully("CModel::UpdateHRUPET: Invalid Evaporation Type",BAD_DATA); break;
  }
  }

  if (PET<(-REAL_SMALL)){
    string warn="negative PET ("+to_string(PET)+" mm/d) calculated in CModel::UpdateHRUPET";
    WriteWarning(warn,false);
    PET=0.0;
  }

  double veg_corr=pHRU->GetVegetationProps()->PET_veg_corr;
  return PET*veg_corr;
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates ground air resistance [s/m]
/// \param &G [in] Reference to ground surface properties
/// \param &VV [in] Reference to vegetation properties
/// \param &wind_vel [in] Wind velocity [m/s]
/// \return Ground air resistance [s/m]
/// from Brook90
//
double GroundAirResistance(const surface_struct &G,   //ground surface
                           const veg_var_struct &VV,
                           const double         &wind_vel)
{
  double Rga;  //ground-air resistance (aka RAS), s/m
  double closed_roughness,closed_zerodisp;
  double ustar,KH;

  //function of height, ground roughness, refht, z0, canopy rough,
  closed_roughness=CVegetationClass::CalcClosedRoughness(VV.height);
  closed_zerodisp =CVegetationClass::CalcClosedZeroPlaneDisp(VV.height,closed_roughness);
  upperswap(closed_roughness,G.roughness);

  ustar=VON_KARMAN*wind_vel/(log((VV.reference_height-VV.zero_pln_disp)/VV.roughness));
  KH   =VON_KARMAN*ustar*(VV.height-VV.zero_pln_disp);

  Rga  =VV.height/(WIND_EXTINCT*KH)*exp(WIND_EXTINCT * (1.0-G.roughness/VV.height));
  Rga-=exp(-WIND_EXTINCT*(closed_roughness + closed_zerodisp)/VV.height);
  upperswap(Rga,1);
  return Rga;
}

//////////////////////////////////////////////////////////////////
/// \brief Returns PET [mm/d] and ground evaporation rate [mm/d] using Shuttleworth Wallace method
/// \ref From Shuttleworh and Wallace, 1985 \cite shuttleworth1985QJRMS
/// \ref as implemented in  Brook90 routines SWGRA and SWPE
///
/// \param &F [in] Reference to model forcing functions
/// \param &matric_pot [in] Matric potential
/// \param &S [in] Reference to soil properties class
/// \param &G [in] Reference to surface properties class
/// \param &VV [in] Reference to vegetation properties class
/// \return Calculated PET [mm/d]
//
double ShuttleworthWallaceEvap(const force_struct   *F,
                               const double         &matric_pot,
                               const soil_struct    &S,
                               const surface_struct &G,   //ground surface
                               const veg_var_struct &VV)        //canopy
{
  double gamma;     //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]
  double vapor_def; //vapor deficit [kPa]

  //following should be calculated as canopy/ground parameters
  double Rss=0.5;//ground-evap resistance, s/m
  Rss=CSoilClass::CalcSoilResistance(matric_pot,S);

  double Rsc;//canopy surface resistance, s/m
  Rsc=1.0/VV.canopy_conductance*MM_PER_METER;

  double Rga; //ground-air resistance (aka RAS), s/m
  //Rga=GroundAirResistance(G,VV,F->wind_vel);
  double closed_roughness,closed_zerodisp;
  double ustar,KH;
  double zero_pl=VV.zero_pln_disp;

  //function of height, ground roughness, refht, z0, canopy rough,
  closed_roughness=CVegetationClass::CalcClosedRoughness(VV.height);
  closed_zerodisp =CVegetationClass::CalcClosedZeroPlaneDisp(VV.height,closed_roughness);
  upperswap(closed_roughness,G.roughness);

  ustar=VON_KARMAN*F->wind_vel/(log((VV.reference_height-zero_pl)/VV.roughness));
  KH   =VON_KARMAN*ustar*(VV.height-zero_pl);

  Rga  =VV.height/(WIND_EXTINCT*KH)*exp(WIND_EXTINCT * (1.0-G.roughness/VV.height));
  Rga-=exp(-WIND_EXTINCT*(closed_roughness + closed_zerodisp)/VV.height);
  upperswap(Rga,1);//why?

  double Raa;//boundary layer resistance, s/m
  Raa =log((VV.reference_height- zero_pl)/(VV.height-zero_pl))/(VON_KARMAN*ustar);
  Raa+=(VV.height / (WIND_EXTINCT * KH)) * (-1.0 + exp(WIND_EXTINCT * (VV.height - closed_roughness - closed_zerodisp) / VV.height));

  double Rac;//leaf-air resistance, s/m
  double leaf_width=1.0; //leaf width, m
  double UH = (ustar / VON_KARMAN) * log((VV.height - zero_pl) / VV.roughness );
  double RB = (100.0 * WIND_EXTINCT) * sqrt(leaf_width / UH) / (1.0 - exp(-WIND_EXTINCT / 2.0));
  Rac = RB / (LEAF_PROJ_RAT * VV.LAI + PI * VV.SAI);

  //Replace above with
  //CalculateResistances(G,VV,F->wind_vel,Raa,Rac,Rga);
  /*CalculateResistances(const surface_struct &G,   //ground surface
    const veg_var_struct &VV,
    const double         &wind_vel,
    double &R_boundary,
    double &R_leaf_air,
    double &R_ground_air);*/

  sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
  de_dT    =GetSatVapSlope           (F->temp_ave,sat_vap);
  LH_vapor =GetLatentHeatVaporization(F->temp_ave);
  gamma    =GetPsychometricConstant  (F->air_pres,LH_vapor);
  vapor_def=sat_vap*(1.0-(F->rel_humidity));

  double AE; //Availabale energy at canopy top, [MJ/d/m2]
  double AE_grnd;//Available energy at ground, [MJ/d/m2]
  AE = (F->SW_radia_net + F->LW_radia);
  AE_grnd =AE*VV.skyview_fact;

  double Rsoil=(de_dT+gamma)*Rga+gamma*Rss; //ground,canopy, and leaf resistances [Kpa-s/K-m]
  double Rcan =(de_dT+gamma)*Rac+gamma*Rsc;
  double Ratm =(de_dT+gamma)*Raa;

  double Ccs=1.0/(1.0+Rsoil*Ratm/(Rcan *(Rsoil+Ratm)));//conductances [unitless]
  double Ccc=1.0/(1.0+Rcan *Ratm/(Rsoil*(Rcan +Ratm)));

  double vpd1=vapor_def-de_dT*Rga*(AE - AE_grnd)/HCP_AIR;//revised vapor deficits [kPa]
  double vpd2=vapor_def-de_dT*Rac*(     AE_grnd)/HCP_AIR;

  double PMS =((Raa+Rga)*de_dT*AE+HCP_AIR*vpd1) / ((de_dT+gamma) * (Raa + Rga) + gamma * Rss);
  double PMC =((Raa+Rac)*de_dT*AE+HCP_AIR*vpd2) / ((de_dT+gamma) * (Raa + Rac) + gamma * Rsc);
  double LE=Ccc*PMC+Ccs*PMS;//total latent heat flux density

  double vpd3 = vapor_def + Raa*(de_dT*AE-(de_dT+gamma)*LE)/HCP_AIR;

  double PET,EVAP;
  PET =(Rac*de_dT*(AE - AE_grnd)+HCP_AIR*vpd3)/((de_dT+gamma)*Rac+gamma*Rsc);
  PET =PET/(DENSITY_WATER*LH_vapor)*MM_PER_METER; //mm/d=[MJ/d/m2]*[m3/kg]*[kg/MJ]*[mm/m]=mm/d

  EVAP =(Rga*de_dT*(     AE_grnd)+HCP_AIR*vpd3)/((de_dT+gamma)*Rga+gamma*Rss);
  EVAP =EVAP/(DENSITY_WATER*LH_vapor)*MM_PER_METER;

  return PET;
}

//////////////////////////////////////////////////////////////////
/// \brief Adjust wind velocity
/// \details returns ratio of wind speed at reference height (above canopy) to
/// wind speed at weather station.
/// \ref Adapted from Brook90s adaptation of Brutsaert 1982 \cite Federer2010
/// \param *F [in] Current model forcing functions
/// \param *VV [in] Reference to vegetation properties of current HRU
/// \param ws_fetch [in] Weather station fetch [m]
/// \param ws_measurement_ht [in] Weather station measurement height for wind, above any zero plane [m]
/// \param ws_roughness [in] Weather station roughness parameter
/// \return Adjusted wind velocity ratio (v_{refht}/v_{met station})
// ROUTINE CURENTLY UNUSED
//
double AdjustWindVel(const force_struct *F,
                     const veg_var_struct *VV,
                     const double ws_fetch,//weather station fetch, m
                     const double ws_measurement_ht,//weather station measurement height for wind,above any zero plane, m
                     const double ws_roughness)//weather station roughness parameter, m
{
  double tmp;
  if (ws_measurement_ht>1.0)
  {
    //Brutsaert (1982) equation 7-39
    double HIBL=0.334*pow(ws_fetch,0.875)*pow(ws_roughness,0.125);//height of internal boundary layer, m
    //Brutsaert equations 7-41 and 4-3
    tmp =log(HIBL / ws_roughness);
    tmp*=log((VV->reference_height - VV->zero_pln_disp) / VV->roughness);
    tmp/=(log(HIBL / VV->roughness) * log(ws_measurement_ht / ws_roughness));
    return tmp;
  }
  else{
    return 1.0;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates ground evaporation using Shuttleworth Wallace method
/// \details Shuttleworth and Wallace (1985) ground evaporation when transpiration known (adapted from Brook90 SWGE)
/// \param actual_transpiration Actual transpiration [mm/d]
/// \param &F [in] Reference to model forcing functions
/// \param &matric_pot [in] Matric potential
/// \param &S [in] Reference to soil properties class
/// \param &G [in] Reference to surface properties class
/// \param &VV [in] Reference to vegetation properties class
/// \return Calculated evaporation rate [mm/d]
//
double GroundEvaporation(const double actual_transpiration,//[mm/d]
                         const force_struct   *F,
                         const double         &matric_pot,
                         const soil_struct    &S,
                         const surface_struct &G,   //ground surface
                         const veg_var_struct &VV)      //canopy
{
  double gamma;       //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]
  double vapor_def; //vapor deficit [kPa]

  double Rga=GroundAirResistance(G,VV,F->wind_vel);

  double Rss;//ground-evap resistance, s/m
  Rss=CSoilClass::CalcSoilResistance(matric_pot,S);

  double Raa;//boundary layer resistance, s/m
  Raa=0.1;   //CalcBoundaryLayerResistance();

  sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
  de_dT    =GetSatVapSlope           (F->temp_ave,sat_vap);
  LH_vapor =GetLatentHeatVaporization(F->temp_ave);
  gamma    =GetPsychometricConstant  (F->air_pres,LH_vapor);
  vapor_def=sat_vap*(1.0-(F->rel_humidity));

  double trans=actual_transpiration*DENSITY_WATER*LH_vapor; //[mm/d]->[MJ/m2/d]

  double AE;     //Available energy at canopy top, [MJ/d/m2]
  double AE_grnd;//Available energy at ground, [MJ/d/m2]
  AE = (F->SW_radia_net + F->LW_radia);
  AE_grnd =AE*VV.skyview_fact;

  double Rs = (de_dT + gamma) * Rga + gamma * Rss;
  double Ra = (de_dT + gamma) * Raa;
  double LE = (Rs* trans + HCP_AIR*vapor_def + de_dT*(Rga*AE_grnd + Raa*AE)) / (Rs + Ra);

  return (LE - trans)/DENSITY_WATER/LH_vapor;//[mm/d]

}

//////////////////////////////////////////////////////////////////
/// \brief Calculates snow evaporation using Shuttleworth Wallace method
/// \details Shuttleworth and Wallace (1985) potential snow evaporation (adapted from Brook90 SWGE)
/// \remark Not currently used
/// \param &F [in] Reference to model forcing functions
/// \param &snow_temp [in] Temperature of snow [C]
/// \param &G [in] Reference to surface properties class
/// \param &VV [in] Reference to vegetation properties class
/// \return Calculated evaporation rate of snow [mm/d]
//
double SnowEvaporation(const force_struct   *F,
                       const double         &snow_temp,
                       const surface_struct &G,   //ground surface
                       const veg_var_struct &VV)      //canopy
{
  double gamma;                   //psychometric "constant" [kPa/K]
  double LH_vapor;        //latent heat of vaporization [MJ/kg]
  double sat_vap;   //Saturation vapor pressure [kPa]

  double Rga=GroundAirResistance(G,VV,F->wind_vel);//s/m
  double Raa;//boundary layer resistance, s/m
  Raa=0.1;//CalcBoundaryLayerResistance();

  sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
  LH_vapor =GetLatentHeatVaporization(F->temp_ave);
  gamma    =GetPsychometricConstant  (F->air_pres,LH_vapor);

  double esnow=0.61;//kPa
  if (F->temp_ave<-0.1){
    esnow = GetSaturatedVaporPressure(min(F->temp_ave,snow_temp));
  }
  double tmp;
  tmp=0.3*(MM_PER_METER/LH_SUBLIM/DENSITY_WATER)*(HCP_AIR/gamma)*(esnow-sat_vap)/(Raa+Rga);
  //[-][kg/MJ]*[m3/kg]*[J/m3/K]*[K/kPa]*[kPa]*[m/s]-->[mm/s*J/MJ]
  return tmp*MJ_PER_J*SEC_PER_DAY;//=[m/d]
}

//////////////////////////////////////////////////////////////////
/// \brief Calculates evaporation using Blaney-Criddle method
/// \todo [add funct] Implement Blaney-Criddle evaporation
/// \details from Blaney & Criddle, 1962 \cite Blaney1962
/// \ref Added by GrahamStonebridge, Fall 2011
/// \param *F [in] Model forcing functions
/// \return Calculated PET [mm/d]
//
double Blaneycriddle(const force_struct *F)
{
  double rho;
  rho=1.0;
  ExitGracefully("Blaneycriddle",STUB);
  return rho*(0.46*F->temp_month_ave+8);
}
//////////////////////////////////////////////////////////////////
/// \brief calculates drying power modifier for wind thru canopy
/// \details from Granger \cite ?
/// \ref ported from classEvap in CRHM 
/// \param wind_vel [in] wind velocity in m/s
/// \param veg_ht [in] canopy height in meters
/// \return drying power modifying coefficient
//
double GetDryingPower(const double &wind_vel,const double &veg_ht) //u in m/s; canopy height in m
{
  // returns Drying power f(u) (mm/d/kPa)
  double Z0 = veg_ht*100.0/7.6;
  return (8.19 + 0.22*Z0) + (1.16 + 0.08*Z0)*wind_vel;
}
//////////////////////////////////////////////////////////////////
/// \brief Calculates evaporation using Granger method
/// \details from Granger \cite ?
/// \ref ported from classEvap in CRHM 
/// \param *F [in] Model forcing functions
/// \param *pHRU pointer to HRU
/// \return Calculated PET [mm/d]
//
double EvapGranger(const force_struct *F,const CHydroUnit *pHRU)
{
  //

  double gamma;     //psychometric "constant" [kPa/K]
  double LH_vapor;  //latent heat of vaporization [MJ/kg]
  double de_dT;     //Vapor pressure-temp slope=de*/dT [kPa/K]
  double sat_vap;   //Saturation vapor pressure [kPa]
  double PET=0.0;

  double F_Qg=0.1;//fraction to ground flux
  double Q = (F->SW_radia_net+F->LW_radia)*(1.0 - F_Qg); // CRHM supposes daily value (mm/d)

  if(Q > 0.0)
  {
    sat_vap  =GetSaturatedVaporPressure(F->temp_ave);
    de_dT    =GetSatVapSlope(F->temp_ave,sat_vap);
    LH_vapor =GetLatentHeatVaporization(F->temp_ave);
    gamma    =GetPsychometricConstant(F->air_pres,LH_vapor);

    double veg_ht=pHRU->GetVegVarProps()->height;
    double ea_mod=GetDryingPower(F->wind_vel,veg_ht)*(sat_vap*(1.0-F->rel_humidity));

    double D=0;
    if(ea_mod > 0.0) { D = min(ea_mod / (ea_mod + Q),1.0); }
    double G = 1.0 / (0.793 + 0.2*exp(4.902*D)) + 0.006 * D;

    PET=  (de_dT * Q + gamma * ea_mod)/(de_dT + gamma / G);
  }
  return max(PET,0.0);
}