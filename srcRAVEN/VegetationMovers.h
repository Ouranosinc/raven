/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvCanopyEvap
  CmvCanopySnowEvap
  CmvCanopyDrip
  ----------------------------------------------------------------*/
#ifndef VEG_MOVERS_H
#define VEG_MOVERS_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

////////////////////////////////////////////////////////////////////
/// \brief Models of canopy evaporation
//
enum canevap_type
{
  CANEVP_RUTTER,                      ///< rutter conceptual model
  CANEVP_MAXIMUM,                     ///< evaporates at PET rate
  CANEVP_ALL                            ///< HBV model
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for canopy evaporation
/// \details Calculates loss of water from canopy to atmosphere
//
class CmvCanopyEvap: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  canevap_type                   type; ///< Type of canopy evapormyation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvCanopyEvap(canevap_type eval_type);        //general constructor
  ~CmvCanopyEvap();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(canevap_type eval_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for canopy snowpack evaporation
/// \details Calculates loss of water from canopy snowpack to atmosphere
//
class CmvCanopySnowEvap: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  canevap_type                   type; ///< Type of canopy snowpack evaporation

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvCanopySnowEvap(canevap_type eval_type);    //general constructor
  ~CmvCanopySnowEvap();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(canevap_type eval_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

////////////////////////////////////////////////////////////////////
/// \brief Models of canopy drip to land surface
//
enum candrip_type
{
  CANDRIP_RUTTER,   ///< drip rate determined solely by canopy capacity overflow
  CANDRIP_SLOWDRAIN ///< drip rate linearly proportional to storage
};

////////////////////////////////////////////////////////////////////
/// \brief Data abstraction for loss of water from canopy to land surface
//
class CmvCanopyDrip: public CHydroProcessABC
{
private:/*------------------------------------------------------*/

  candrip_type  type; ///< Model of canopy drip used

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvCanopyDrip(candrip_type drip_type,
                int                                      to_index);
  ~CmvCanopyDrip();

  //inherited functions
  void Initialize();
  void GetRatesOfChange(const double              *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;
  void ApplyConstraints(const double      *state_vars,
                        const CHydroUnit  *pHRU,
                        const optStruct   &Options,
                        const time_struct &tt,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(candrip_type drip_type,
                                           sv_type *aSV, int *aLev, int &nSV);
};

#endif

