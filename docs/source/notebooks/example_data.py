#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 10:25:14 2019

@author: ets
"""

from pathlib import Path

TESTS_HOME = Path(__file__).parent.parent.parent.parent
TD = TESTS_HOME / 'tests' / 'testdata'
CFG_FILE = [str(TESTS_HOME / 'test.cfg'), ]

TESTDATA = dict()
TESTDATA['gr4j-cemaneige'] = \
    {'pr': TD / 'gr4j_cemaneige' / 'pr.nc',
     'tas': TD / 'gr4j_cemaneige' / 'tas.nc',
     'evap': TD / 'gr4j_cemaneige' / 'evap.nc'}

TESTDATA['raven-gr4j-cemaneige-nc-ts'] = TD / 'raven-gr4j-cemaneige' / 'Salmon-River-Near-Prince-George_meteo_daily.nc'
TESTDATA['raven-gr4j-cemaneige-nc-rv'] = tuple((TD / 'raven-gr4j-cemaneige').glob('raven-gr4j-salmon.rv?'))

TESTDATA['raven-mohyse-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-mohyse'] = TD / 'raven-mohyse'
TESTDATA['raven-mohyse-rv'] = tuple((TD / 'raven-mohyse').glob('raven-mohyse-salmon.rv?'))
TESTDATA['raven-mohyse-ts'] = tuple((TD / 'raven-mohyse').glob('Salmon-River-Near-Prince-George_*.rvt'))

TESTDATA['raven-hmets-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-hmets'] = TD / 'raven-hmets'
TESTDATA['raven-hmets-rv'] = tuple((TD / 'raven-hmets').glob('raven-hmets-salmon.rv?'))
TESTDATA['raven-hmets-ts'] = tuple((TD / 'raven-hmets').glob('Salmon-River-Near-Prince-George_*.rvt'))

TESTDATA['raven-hbv-ec-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']
TESTDATA['raven-hbv-ec'] = TD / 'raven-hbv-ec'
TESTDATA['raven-hbv-ec-rv'] = tuple((TD / 'raven-hbv-ec').glob('raven-hbv-ec-salmon.rv?'))
TESTDATA['raven-hbv-ec-ts'] = tuple((TD / 'raven-hbv-ec').glob('Salmon-River-Near-Prince-George_*.rvt'))

TESTDATA['ostrich-gr4j-cemaneige'] = TD / 'ostrich-gr4j-cemaneige'
TESTDATA['ostrich-gr4j-cemaneige-rv'] = \
    tuple(TESTDATA['ostrich-gr4j-cemaneige'].glob("*.rv?")) + tuple(TESTDATA['ostrich-gr4j-cemaneige'].glob('*.t??'))
TESTDATA['ostrich-gr4j-cemaneige-nc-ts'] = TESTDATA['raven-gr4j-cemaneige-nc-ts']

TESTDATA['ostrich-mohyse'] = TD / 'ostrich-mohyse'
TESTDATA['ostrich-mohyse-rv'] = \
    tuple(TESTDATA['ostrich-mohyse'].glob("*.rv?")) + tuple(TESTDATA['ostrich-mohyse'].glob('*.t??'))
TESTDATA['ostrich-mohyse-nc-ts'] = TESTDATA['raven-mohyse-nc-ts']

TESTDATA['ostrich-hmets'] = TD / 'ostrich-hmets'
TESTDATA['ostrich-hmets-rv'] = \
    tuple(TESTDATA['ostrich-hmets'].glob("*.rv?")) + tuple(TESTDATA['ostrich-hmets'].glob('*.t??'))
TESTDATA['ostrich-hmets-nc-ts'] = TESTDATA['raven-hmets-nc-ts']

TESTDATA['ostrich-hbv-ec'] = TD / 'ostrich-hbv-ec'
TESTDATA['ostrich-hbv-ec-rv'] = \
    tuple(TESTDATA['ostrich-hbv-ec'].glob("*.rv?")) + tuple(TESTDATA['ostrich-hbv-ec'].glob('*.t??'))
TESTDATA['ostrich-hbv-ec-nc-ts'] = TESTDATA['raven-hbv-ec-nc-ts']

TESTDATA['donnees_quebec_mrc_poly'] = TD / 'donneesqc_mrc_poly' / 'donnees_quebec_mrc_polygones.gml'
TESTDATA['watershed_vector'] = TD / 'watershed_vector' / 'LSJ_LL.zip'
TESTDATA['mrc_subset'] = TD / 'donneesqc_mrc_poly' / 'mrc_subset.gml'
TESTDATA['melcc_water'] = TD / 'melcc_water_management' / 'zone_gestion_leau_saintlaurent.gpkg'

# TODO: Replace the following files with subsets and set originals as production data
TESTDATA['earthenv_dem_90m'] = TD / 'earthenv_dem_90m' / 'earthenv_dem90_southernQuebec.tiff'
TESTDATA['hydrobasins_lake_na_lev12'] = TD / 'usgs_hydrobasins' / 'hybas_lake_na_lev12_v1c.zip'
TESTDATA['simfile_single'] = TD / 'hydro_simulations' / 'raven-gr4j-cemaneige-sim_hmets-0_Hydrographs.nc'
TESTDATA['input2d']=TD / 'input2d' / 'input2d.nc'