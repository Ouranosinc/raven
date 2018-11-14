/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2017 the Raven Development Team
  ----------------------------------------------------------------
  class definitions:
  CmvAbstraction
  ----------------------------------------------------------------*/

#ifndef ABSTRACTION_H
#define ABSTRACTION_H

#include "RavenInclude.h"
#include "HydroProcessABC.h"

///////////////////////////////////////////////////////////////////
/// \brief Method of calculating abstraction of rainfall/snowmelt
//
enum abstraction_type{
  ABST_PERCENTAGE, ///< Abstraction by percentage
  ABST_FILL,       ///< Fill abstraction
  ABST_SCS         ///< Abstraction using SCS method
};

///////////////////////////////////////////////////////////////////
/// \brief Calculates the abstraction of rainfall/snowmelt (ponded water) to depression storage
//
class CmvAbstraction: public CHydroProcessABC
{
private:/*------------------------------------------------------*/
  abstraction_type type; ///< Model of abstaction selected

public:/*-------------------------------------------------------*/
  //Constructors/destructors:
  CmvAbstraction(const abstraction_type absttype);
  ~CmvAbstraction();

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
                        const time_struct &t,
                        double      *rates) const;

  void        GetParticipatingParamList   (string  *aP, class_type *aPC, int &nP) const;
  static void GetParticipatingStateVarList(abstraction_type  snalbtype,
                                           sv_type *aSV,
                                           int     *aLev,
                                           int     &nSV);
};
#endif
