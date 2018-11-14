/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2018 the Raven Development Team
  ----------------------------------------------------------------*/
#ifndef UNITTESTING_H
#define UNITTESTING_H

#include "RavenInclude.h"

void DateTest();
void OpticalAirMassTest();
void ClearSkyTest();
void ShortwaveTest();
void ShortwaveGenerator();
void JulianConvertTest();
void SmartLookupUnitTest();
void SmartIntervalTest();
void GammaTest();

void RavenUnitTesting(const optStruct &Options);
#endif
