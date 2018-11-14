/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------
  Properties.h
  ------------------------------------------------------------------
  Contains definitions of property structures and declarations of
  default methods for their creation
  should also get a routine for checking validity of params
  visible to all classes
  struct soil_struct,
  struct veg_struct, veg_var_struct,
  struct surface_struct
  default methods for their creation.
  Should also get a routine for checking validity of params
  visible to all classes:
  - struct soil_struct,
  - struct veg_struct, veg_var_struct,
  - struct surface_struct
  ----------------------------------------------------------------*/
#ifndef DEFAULT_PROPS_H
#define DEFAULT_PROPS_H

#include "RavenInclude.h"
struct force_struct;

///////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects a specific soil type.
/// \details Soil properties calculated from soiltype.
/// \note If additional properties are added, default values/formulae for values need to be
/// specified in CSoilClass::AutoCalculateSoilProps
//
struct soil_struct
{
  //specified (required) parameters
  double org_con;           ///< [0..1]    percentage organic content
  double sand_con;          ///< [0..1]    percentage sand content
  double clay_con;          ///< [0..1]    percentage clay content

  //specified or automatically calculated parameters
  //standard soil parameters
  double porosity;          ///< [0..1]    soil porosity
  double stone_frac;        ///< [0..1]    stone volume fraction
  double bulk_density;      ///< [kg/m3]   bulk dry density
  //double specific_storage;///< [-]       specific storage

  double cap_ratio;         ///< [-]       total water storage cap.=thick*cap_ratio ; cap_ratio=poro*(1-stone_f)

  //Thermal properties
  double heat_capacity;     ///< [J/m3/K]   saturated *volumetric* heat capacity
  double thermal_cond;      ///< [W/m/K]    saturated soil thermal conductivity

  //unsaturated flow parameters
  double hydraul_cond;      ///< [mm/d]  saturated hydraulic conductivity
  double clapp_b;           ///< [-]       clapp-hornberger exponent  \math \f$ {\psi}={\psi_a}*S^{-b} \f$
  double clapp_n;           ///< [-]       clapp-hornberger transition parameters for
  double clapp_m;           ///< [mm]      high saturation: \math \f$ {\psi}=m * (S - n) * (S - 1)\f$
  double sat_res;           ///< [0..1]    Hydroscopic minimum saturation / residual saturation
  double sat_wilt;          ///< [0..1]    wilting point saturation
  double field_capacity;    ///< [0..1]    saturation at field capacity
  double air_entry_pressure;///< [-mm]     air entry pressure (a positive value)
  double wilting_pressure;  ///< [-mm]     Wilting pressure
  double wetting_front_psi; ///< [-mm]     Wetting front matric potential (for GA model)
  double ksat_std_deviation;///< [-]       standard deviation in log(k_sat)

  //evaporation parameters
  double evap_res_fc;       ///< [d/mm]     soil evaporation resistance at field capacity
  double shuttleworth_b;    ///< [-]        exponent in relation of soil evap res to water potential \math \f$ Res = {evap_res_fc} * ({\Psi}/{\psi_{fc}})^b \f$
  double PET_correction;    ///< [-]        fraction of PET used for soil evap

  //albedo parameters
  double albedo_wet;        ///< [-]        saturated soil albedo
  double albedo_dry;        ///< [-]        bone dry soil albedo

  //special params (for specific models)
  double VIC_zmin;          ///< [mm]      Xinanjiang parameters for VIC infiltration model \ref (Woods et al, 1992) \cite Wood1992JoGR
  double VIC_zmax;          ///< [mm]
  double VIC_alpha;         ///< [-]
  double VIC_evap_gamma;    ///< [-]       power law exponent for VIC soil evaporation

  double b_exp;             ///< [-]       power law exponent for VIC_ARNO infiltration model \ref (Clark et al, 2008)

  double max_perc_rate;     ///< [mm/d]    VIC/ARNO/GAWSER percolation rate - user defined between 0.010 - 1000.0
  double perc_coeff;        ///< [1/d]     Linear percolation storage coefficient
  double perc_n;            ///< [-]       VIC/ARNO percolation exponent - user defined between 1.00 - 20.00
  double SAC_perc_alpha;    ///< [-]       Sacremento percolation multiplier - user defined between 1.00 - 250.00
  double SAC_perc_expon;    ///< [-]       Sacremento percolation exponent - user defined between 1.00 - 5.00
  double perc_aspen;        ///< [mm/d]    constant max percolation rate for PERC_ASPEN model (S. Grass, 2018)

  double max_baseflow_rate; ///< [mm/d]    max baseflow rate (e.g., VIC_ARNO)- user defined between 0.001 - 10000.00
  double baseflow_n;        ///< [-]       VIC/ARNO baseflow exponent - user defined between 1.0 - 10.0
  double baseflow_coeff;    ///< [1/d]     Linear Baseflow storage coefficient [Q_base=K*S]
  double baseflow_thresh;   ///< [0..1]    threshold saturation for onset of baseflow

  double max_cap_rise_rate; ///< [mm/d]    HBV max capillary rise rate

  double max_interflow_rate;///< [mm/d]    PRMS max_interflow rate
  double interflow_coeff;   ///< [1/d]     Linear Interflow storage coefficient

  double HBV_beta;          ///< [?]       soil infiltration param from HBV (move to Surface struct?)

  double UBC_evap_soil_def; ///< [mm]      soil deficit at which AET=0.1*PET (UBCWM P0EGEN)
  double UBC_infil_soil_def;///< [mm]      soil deficit at which effective impermeable fraction depletes to 0.1 (UBCWM P0AGEN)

  double GR4J_x2;           ///< [mm/d]    GR4J maximum groundwater exchange rate
  double GR4J_x3;           ///< [mm]      GR4J reference storage for baseflow/GW exchange

  double exchange_flow;     ///< [mm/d]    exchange flow rate with mixing zone (greater than 0)

  //Transport parameters (all per soil/constituent combination)
  double retardation   [MAX_CONSTITUENTS]; ///< [-] soil-specific retardation factor
  double mineraliz_rate[MAX_CONSTITUENTS]; ///< [1/d] soil-specific mineralization rate 
  double loss_rate     [MAX_CONSTITUENTS]; ///> [1/d] soil-specific unspecified linear loss rate
  double transf_coeff  [MAX_CONSTITUENTS][MAX_CONSTITUENTS]; ///< [1/d] linear transformation coefficients (one per species-species combination) 
  double stoichio_coeff[MAX_CONSTITUENTS][MAX_CONSTITUENTS]; ///< [1/d] linear transformation coefficients (one per species-species combination) 
  
};

///////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects a specific type of vegetation
/// \details Contains vegetation properties calculated from vegetation class
/// \note If additional properties are added, default values/formulae for values need to be
/// specified in CVegetationClass::AutoCalculateSoilProps
//
struct veg_struct
{
  string vegetation_name;   ///<           name of vegetation class

  //specified (required) properties
  double max_height;        ///< [m]       maximum vegetation height
  double max_leaf_cond;     ///< [mm/s]    maximum leaf conductance
  double max_LAI;           ///< [m2/m2]   maximum Leaf Area Index

  //estimable parameters - usually unchanged
  double svf_extinction;    ///< [-]       extinction coefficient used to calculate skyview factor; SVF=exp(-extinct*LAI)
  double rain_icept_fact;   ///< [-]       relates percentage of throughfall of rain to LAI+SAI (~0.06)
  double snow_icept_fact;   ///< [-]       relates percentage of throughfall of snow to LAI+SAI (~0.04)
  double rain_icept_pct;    ///< [0..1]    percentage of rain intercepted (maximum) (only on canopy portion)
  double snow_icept_pct;    ///< [0..1]    percentage of snow intercepted (maximum) (only on canopy portion)

  double SAI_ht_ratio;      ///< [m2/m3]   ratio of stem area index to height
  double trunk_fraction;    ///< [0..1]    fraction of canopy attributed to tree trunk
  double stemflow_frac;     ///< [0..1]    ~0.03
  //double max_bottom_elev; ///< [m]       bottom of canopy (e.g., CLM v3.5)

  double albedo;            ///< [-]       visible/near-infrared albedo of leaf
  double albedo_wet;        ///< [-]       albedo of wet leaf

  double relative_ht[12];   ///< [m/m]     relative vegetation height over 12 mos (ht=rel_ht*max_ht)
  double relative_LAI[12];  ///< [m/m]     relative vegetation LAI over 12 mos

  double max_capacity;      ///< [mm]      maximum canopy storage capacity
  double max_snow_capacity; ///< [mm]      maximum canopy snow storage capacity

  double veg_dens;          ///< [1/m2]    vegetation count per meter squared (range: 1.0-500.0; recommended 300.0 for crops and grass, 1.0 for forests and shrubs)
  double veg_diam;          ///< [m]       vegetation diameter (range 0.0-2.0; recommended 0.003 for crops and grass, 0.5-1.0 for forests and shrubs [m]
  double veg_mBeta;         ///< [-]       mBeta parameter

  double PET_veg_corr;      ///< [0..1]    vegetation-based PET correction (multiplicative)  

  //root properties
  double root_extinct;      ///< [-]       extinction coefficient for roots, exp(-ext*z)
  double max_root_length;   ///< [mm/m2]   root length per unit canopy area
  double min_resistivity;   ///< [d/mm]    1.0/max_conductivity
  double xylem_frac;        ///< [0..1]    fraction of plant resistance in xylem (used to calculate xylem resistance)
  double rootradius;        ///< [mm]      average root radius (used to calculate cowan alpha)
  double psi_critical;      ///< [-mm]     minimum plant leaf water potential
                            //                the critical potential at which stomates close

  //special process-specific params (may not be used)
  double drip_proportion;   ///< [1/d]     drip proportion for bucket drip model: \math \f$ dS/dt=-dp*S \f$
  double max_intercept_rate;///< [mm/d]    maximum rate of rainfall interception
  double CHU_maturity;      ///< [-]       crop heat unit maturity; level at which PET is maximized

  //Transport parameters 
  double uptake_moderator   [MAX_CONSTITUENTS]; ///< [-] vegetation-specific uptake factor
};

///////////////////////////////////////////////////////////////////
/// \brief Contains derived vegetation / canopy / root properties - unique for each HRU
//
struct veg_var_struct
{
  //derived properties (temporally variable)-------------------------------------------
  //Updated in CVegegtationClass::RecalculateCanopyParams
  double LAI;               //< [m2/m2] Leaf Area Index
  double SAI;               //< [m2/m2] Stem Area Index
  double height;            //< [m]     vegetation height
  double capacity;          //< [mm]    rain storage capacity
  double snow_capacity;     //< [mm SWE]snow storage capacity
  double leaf_cond;         //< [mm/s]  leaf conductance (same as stomatal conductance?)
  double canopy_conductance;//< [mm/s]  leaf conductance corrected for shelter & LAI
  double rain_icept_pct;    //< [0..1]  percentage of rain intercepted  (only on canopy portion)
  double snow_icept_pct;    //< [0..1]  percentage of snow intercepted  (only on canopy portion)

  double shelter_factor;    ///< [-]    accounts for sheltered leaves : about 0.5-1.0
  double skyview_fact;      ///< [0..1] skyview factor, pct of ground visible from sky

  double roughness;         ///< [m]    surface roughness parameter for momentum transfer
  double zero_pln_disp;     ///< [m]    zero-plane displacement, height where wind vel. goes to zero
  double reference_height;  ///< [m]    reference height for air properties above ground

  //derived root properties (temporally variable)-------------------------------------------
  double root_length;       ///< [mm/m2]root length per unit land area
  double resistivity;       ///< [d/mm] plant resistance to water flow

  //should be portioned out by soil layer (made array in HRU?)
  double rel_rootden;       ///< [m/m3] relative values of root length per unit volume
  double root_resistance;   ///< [d]    root resistance for layer
  double cowan_alpha;       ///< [MPa]  modified Cowan alpha for layer

};

//////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects the qualites of the surface
/// \note All are calculated from LULT and fixed with time
//
struct surface_struct
{
  string landuse_name;      ///<           name of land use linked to this data

  double impermeable_frac;  ///< [-]       fraction of surface that is impermeable

  double roughness;         ///< [m]       roughness parameter of ground surface
                            //             used for evapotranspiration calcs
  double forest_coverage;   ///< [0..1]    fraction of land covered by canopy
  double forest_sparseness; ///< [0..1]    sparseness of canopy in land covered by forest
  double UBC_icept_factor;  ///< [0..1]    effective forest cover for interception (C0TREE*C0CANY)
  double wind_exposure;     ///< [0..1]    basin wind coefficient (1 for unforested areas)

  //snow parameters
  double melt_factor;       ///< [mm/d/C]  maximum snow melt factor used in degree day and hybrid snowmelt models
  double min_melt_factor;   ///< [mm/d/C]  minimum snow melt factor used in degree day and hybrid snowmelt models
  double max_melt_factor;   ///< [mm/d/C]  maximum snow melt factor used in degree day and hybrid snowmelt models
  double DD_melt_temp;      ///< [C]       melting base temperature in degree day methods
  double DD_aggradation;    ///< [1/mm]    rate of increase of melt factor with cumulative snowmelt (Kcum in HMETS)
  double refreeze_factor;   ///< [mm/d/C]  maximum refreeze factor used in degree day and hybrid snowmelt models
  double DD_refreeze_temp;  ///< [C]       degree day refreeze temperature (refreeze base temperature)
  double refreeze_exp;      ///< [-]       empirical exponent for DD freezing equation 
  double HBV_melt_for_corr; ///< [-]       HBV snowmelt forest correction (MRF in HBV-EC)
  double HBV_melt_asp_corr; ///< [-]       HBV snowmelt aspect correction (AM in HBV-EC)
  double snow_patch_limit;  ///< [mm]      SWE limit below which snow does not completely cover ground.  Used as a threshold for snow-patching algorithms (default=0.0).
  double fetch;             ///< [m]       distance of unobstructed wind flow (g.t. 300m)
  double conv_melt_mult;    ///< [-]       convection melt multiplier
  double cond_melt_mult;    ///< [-]       condensation melt multiplier
  double rain_melt_mult;    ///< [-]       rain melt multiplier

  //Glacier parameters
  double glac_storage_coeff;  ///< [-]     maximum linear storage coefficient for glacial melt =K*G
  double HBV_melt_glacier_corr;///< [-]    degree day correction factor for glacial melt (MRG in HBV-EC)
  double HBV_glacier_Kmin;  ///< [-]       minimum linear storage coefficient for glacial melt =K*G (KG1 in HBV-EC)
  double HBV_glacier_Ag;    ///< [1/mm]    extinction coefficient for diminishing storage coefficient with snow depth atop glacier, AG (K=Kmin+(Kmax-Kmin)*exp(-ag*snow_depth))

  double CC_decay_coeff;    ///< [1/d]     linear decay coefficient for decreasing cold content

  //runoff parameters
  double SCS_CN;            ///< [-]       SCS curve number (for antecedent wetness condition II)
  double partition_coeff;   ///< [0..1]    simple partitioning coefficient: percentage of rainfall that runs off
                            //             used for rational method
  double SCS_Ia_fraction;   ///< [0..1]    fraction of rainfall/melt which is initially abstracted to depression storage (default=0.2)
  double HMETS_runoff_coeff;///< [0..1]    HMETS algorithm runoff coeff - maximum fraction of ponded water which can infiltrate

  //surface storage variables
  double max_sat_area_frac; ///< [-]       PRMS maximum saturated area (pct)- user defined between 0.050 - 0.950
  double b_exp;             ///< [-]       ARNO/VIC b exponent - user defined between 0.001 - 3.0

  double dep_max;           ///< [mm]      maximum amount of water that can be stored in depressions
  double abst_percent;      ///< [0..1]    percentage of rainfall/melt which is abstracted to depression storage
  double dep_max_flow;      ///< [mm/d]    outflow rate with full depression storage (dep_stor=dep_max)
  double dep_n;             ///< [-]       power law coefficient for depression outflow
  double dep_threshhold;    ///< [mm]      threshold storage at which flow commences
  double dep_k;             ///< [1/d]     linear overflow coefficient Q=k*(h-dep_threshhold)
  double dep_seep_k;        ///< [1/d]     linear seepage coefficient Q=k*h
  double dep_crestratio;    ///<-[m/m]     ratio of crest width to square root of HRU area (0-5) ~1.5

  double ow_PET_corr;       ///< [-]       fraction of PET to apply to open water evaporation
  double lake_PET_corr;     ///< [-]       fraction of PET to apply to lake evaporation
  double forest_PET_corr;   ///< [-]       fraction of PET to apply to forest evapotranspiration
  double AET_coeff;         ///< [-]       linear AET coefficient

  double lake_rel_coeff;    ///< [1/d]     linear lake storage coefficient 

  double GR4J_x4;           ///< [d]       GR4J time routing parameter
  double gamma_shape;       ///< [-]       gamma distribution shape parameter for convolution routing
  double gamma_scale;       ///< [-]       gamma distribution scaling parameter for convolution routing
  double gamma_shape2;      ///< [-]       gamma distribution shape parameter for convolution routing (v2)
  double gamma_scale2;      ///< [-]       gamma distribution scaling parameter for convolution routing (v2)
};

//////////////////////////////////////////////////////////////////
/// \brief Structure that, when instantiated, contains defining information that reflects the qualities of the terrain
/// \note Calculated from topography, all are fixed with time
//
struct terrain_struct
{
  double drainage_density;  ///< [km/km2] length of *all* rivers in a basin divided by the area of the drainage basin
  double hillslope_length;  ///< [m]      representative hillslope length within the watershed

  double lambda;            ///< [m]      TOPMODEL mean of the log-transformed topographic index (between 5.0 - 10.0)
};

////////////////////////////////////////////////////////////////////
/// \brief UBC lapse rate codes
///
//
struct UBC_lapse
{ //defined in UBCWM manual
  double A0STAB;  ///< [0..1] Precipitation gradient modification factor (default 0.01)
  double A0TLXM;  ///< [�C / 1000 m] Lapse rate for maximum temperatures when  the station elevation is less than 2000 m (default 10)
  double A0TLNM;  ///< [�C / 1000 m] Lapse rate for minimum temperatures when the station elevation is less than 2000 m (default 0.5)
  double A0TLXH;  ///< [�C / 1000 m] Lapse rate for maximum temperatures when the station elevation is greater than 2000 m (default 6.4)
  double A0TLNH;  ///< [�C / 1000 m] Lapse rate for minimum temperatures when the station elevation is greater than 2000 m (default 2.0)
  double P0TEDL;  ///< [�C / 1000 m] Lapse rate of maximum temperature range for elevations below 2000 m (default 6.0)
  double P0TEDU;  ///< [�C / 1000 m] Lapse rate of maximum temperature range for elevations above 2000 m (default 0.0)
  double E0LLOW;  ///< [m] Reference elevation for precipitation gradient P0GRADL
  double E0LMID;  ///< [m] Elevation above which precipitation gradient P0GRADM applies. Set at approx. 1/2 barrier height
  double E0LHI;   ///< [m] Elevation above which the precipitation gradient P0GRADU applies.  Set at approx. 2/3 barrier height
  double P0GRADL; ///< [%, 0..20] Precipitation gradient factor for elevations below E0LMID (default 5)
  double P0GRADM; ///< [%, 0..20] Precipitation gradient factor for elevations below E0LHI (default 0)
  double P0GRADU; ///< [%, 0..20] Precipitation gradient factor for elevations above E0LHI (default 0)
  double A0PELA;  ///< [mm / 1000 m] Lapse rate of potential evapotranspiration factor A0EDDF (default 0.9)
  double max_range_temp; ///< [�C] (A0TERM)
  double A0PPTP;  ///< [mm/d] Threshold precipitation for temperature lapse rate (default 5.0)
};

////////////////////////////////////////////////////////////////////
/// \brief UBC Watershed model snow parameters
///
//
struct UBC_snow_par{
  double ALBASE;          ///< P0ALBASE, albedo exponential decay threshold value (0.65)
  double ALBREC;          ///< P0ALBREC, recession constant (0.9/day)
  double MAX_CUM_MELT;    ///< P0ALBMLX, "max annual melt" (~4000mm)
  double ALBSNW;          ///< P0ALBSNW  total daily snowfall required to bring albedo to that of new snow [mm/d] ,  15
};

////////////////////////////////////////////////////////////////////
/// \brief Stores all global model parameters
///
//
struct global_struct
{
  double           adiabatic_lapse;     ///< [�C/km]
  double           wet_adiabatic_lapse; ///< [�C/km]
  double           precip_lapse;        ///< [mm/d/km] precipitation lapse rate for orographic correction

  double           rainsnow_temp;       ///< [�C] rain/snow halfway transition temperature  (-0.15 for dingman, 0.0 for HBV, ~1 for UBC)
  double           rainsnow_delta;      ///< [�C] range of rain-snow transition zone (around rainsnow_temp) (4.0 for HBV, 2.0 for UBC)
  double           snow_SWI;            ///< [0..1] ~0.05, 0.07 irreducible water saturation fraction of snow 
                                        //   =WHC in HBV-EC, WATCAP in UBCWM, SWI in GAWSER
  double           snow_SWI_min;        ///< [0..1] minimum irreducible sat fraction 
  double           snow_SWI_max;        ///< [0..1] maximum irreducible sat fraction
  double           SWI_reduct_coeff;    ///< [1/mm] SWI reduction factor with cumulative snowmelt 
  double           snow_temperature;    ///< [�C]   default snow temperature if not explicitly modelled
  double           snow_roughness;      ///< [mm]  roughness height of snow
  double           min_snow_albedo;     ///< [0..1] very old snow/glacier albedo (~0.3)
  double           max_snow_albedo;     ///< [0..1] albedo of fresh snow  (~0.95)

  double           avg_annual_snow;     ///< [mm] avg annual snow as SWE
  double           avg_annual_runoff;   ///< [mm] avg annual runoff from basin

  double           max_reach_seglength; ///< [km] maximum reach segment length
  double           max_SWE_surface;     ///< [mm] maximum SWE in surface snow layer

  UBC_lapse        UBC_lapse_params;    ///< parameters for orographic corrections
  UBC_snow_par     UBC_snow_params;     ///< UBC snow parameters
  double           UBC_s_corr[12];      ///< UBC solar radiation corrections
  double           UBC_n_corr[12];      ///<

  double           UBC_GW_split;        ///< UBC GW Split parameter [0..1]
  double           UBC_exposure_fact;   ///< UBC Sun exposure factor of forested areas (F0ERGY) [0..1]
  double           UBC_cloud_penet;     ///< UBC Fraction of solar radiation penetrating cloud cover [0..1]
  double           UBC_LW_forest_fact;  ///< UBC mulitplier of temperature to estimate LW radiation in forests (P0BLUE*P0LWVF) [mm/d/K]
  double           UBC_flash_ponding;   ///< UBC ponding threshold for flash factor (V0FLAS) [mm]

  double           airsnow_coeff;       ///< [1/d] air/snow heat transfer coefficient

  double           MOHYSE_PET_coeff;    ///< [mm/d] MOHYSE PET constant (PET @ equinox when temperature is 10 degC)

  double           TOC_multiplier;      ///< [mm] time of concentration multiplier
};
#endif
