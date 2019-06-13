#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:24:45 2019

@author: ets
"""

import time

import pandas as pd
from birdy import WPSClient

url = "http://localhost:9099/wps"
wps = WPSClient(url, progress=True)

nbBasin = 5797

for i in range(1, nbBasin + 1):

    # Read the current dataset
    tmp = pd.read_csv('gauged_catchment_properties.csv')

    # Read the land use data
    Forest = tmp['forest']
    Shrubs = tmp['shrubs']
    Wetland = tmp['wetland']
    Crops = tmp['crops']
    Urban = tmp['urban']
    Grass = tmp['grass']
    SnowIce = tmp['snowIce']
    Water = tmp['water']
    RunSuccessLandUse = tmp['RunSuccessLandUse']

    # Read the terrain data
    Elevation = tmp['elevation']
    Slope = tmp['slope']
    Aspect = tmp['aspect']
    RunSuccessTerrain = tmp['RunSuccessTerrain']

    # Read the shape data
    Area = tmp['area']
    Latitude = tmp['latitude']
    Longitude = tmp['longitude']
    Gravelius = tmp['gravelius']
    Perimeter = tmp['perimeter']
    RunSuccessShape = tmp['RunSuccessShape']
    StationID = tmp['StationID']
    ID = tmp['ID']

    # Track the catchment number from the list
    ID[i - 1] = i

    # LAND USE ANALYSIS
    if RunSuccessLandUse[i - 1] == 0:
        for attempt in range(1):
            try:

                print('Doing Basin: ' + str(i))
                shape2 = '/home/ets/SHP_ZIP/Basin_' + str(i) + '.zip'

                resp = wps.nalcms_zonal_stats(shape=shape2, select_all_touching=True, band=1, simple_categories=True)

                while resp.isNotComplete():
                    time.sleep(1)

                [statistics] = resp.get(asobj=True)
                lu = statistics['features'][0]['properties']['land-use']
                total = sum(lu.values())
                landUse = {k: v / total for (k, v) in lu.items()}

                Forest[i - 1] = landUse['Forest']
                Grass[i - 1] = landUse['Grass']
                Wetland[i - 1] = landUse['Wetland']
                Water[i - 1] = landUse['Water']
                Urban[i - 1] = landUse['Urban']
                Shrubs[i - 1] = landUse['Shrubs']
                Crops[i - 1] = landUse['Crops']
                SnowIce[i - 1] = landUse['SnowIce']
                RunSuccessLandUse[i - 1] = 1

                # Save file
                zippedList = list(
                    zip(ID, StationID, Area, Latitude, Longitude, Gravelius, Perimeter, RunSuccessShape, Elevation,
                        Slope, Aspect, RunSuccessTerrain, Forest, Grass, Wetland, Water, Urban, Shrubs, Crops, SnowIce,
                        RunSuccessLandUse))
                frameList = pd.DataFrame(data=zippedList,
                                         columns=['ID', 'StationID', 'area', 'latitude', 'longitude', 'gravelius',
                                                  'perimeter', 'RunSuccessShape', 'elevation', 'slope', 'aspect',
                                                  'RunSuccessTerrain', 'forest', 'grass', 'wetland', 'water', 'urban',
                                                  'shrubs', 'crops', 'snowIce', 'RunSuccessLandUse'])
                frameList.to_csv('gauged_catchment_properties.csv')
            except Exception as e:
                print('{}: Failed, restarting iter {}'.format(e, attempt))
            else:
                break

    # TERRAIN ANALYSIS
    if RunSuccessTerrain[i - 1] == 0:
        for attempt in range(1):
            try:

                print('Doing Basin: ' + str(i))
                shape2 = '/home/ets/SHP_ZIP/Basin_' + str(i) + '.zip'

                resp = wps.terrain_analysis(shape=shape2, select_all_touching=True, projected_crs=3978)
                while resp.isNotComplete():
                    time.sleep(1)

                [properties] = resp.get(asobj=True)
                Elevation[i - 1] = properties[0]['elevation']
                Slope[i - 1] = properties[0]['slope']
                Aspect[i - 1] = properties[0]['aspect']
                RunSuccessTerrain[i - 1] = 1

                # Save File
                zippedList = list(
                    zip(ID, StationID, Area, Latitude, Longitude, Gravelius, Perimeter, RunSuccessShape, Elevation,
                        Slope, Aspect, RunSuccessTerrain, Forest, Grass, Wetland, Water, Urban, Shrubs, Crops, SnowIce,
                        RunSuccessLandUse))
                frameList = pd.DataFrame(data=zippedList,
                                         columns=['ID', 'StationID', 'area', 'latitude', 'longitude', 'gravelius',
                                                  'perimeter', 'RunSuccessShape', 'elevation', 'slope', 'aspect',
                                                  'RunSuccessTerrain', 'forest', 'grass', 'wetland', 'water', 'urban',
                                                  'shrubs', 'crops', 'snowIce', 'RunSuccessLandUse'])
                frameList.to_csv('gauged_catchment_properties.csv')

            except Exception as e:
                print('{}: Failed, restarting iter {}'.format(e, attempt))
            else:
                break

    # SHAPE ANALYSIS
    if RunSuccessShape[i - 1] == 0:
        for attempt in range(1):
            try:

                print('Doing Basin: ' + str(i))
                shape2 = '/home/ets/SHP_ZIP/Basin_' + str(i) + '.zip'

                resp = wps.shape_properties(shape=shape2, projected_crs=32198)

                while resp.isNotComplete():
                    time.sleep(1)

                [properties] = resp.get(asobj=True)
                Area[i - 1] = properties[0]['area'] / 1000000.0
                Longitude[i - 1] = properties[0]['centroid'][0]
                Latitude[i - 1] = properties[0]['centroid'][1]
                Gravelius[i - 1] = properties[0]['gravelius']
                Perimeter[i - 1] = properties[0]['perimeter']
                StationID[i - 1] = properties[0]['ID']
                RunSuccessShape[i - 1] = 1

                # Save File
                zippedList = list(
                    zip(ID, StationID, Area, Latitude, Longitude, Gravelius, Perimeter, RunSuccessShape, Elevation,
                        Slope, Aspect, RunSuccessTerrain, Forest, Grass, Wetland, Water, Urban, Shrubs, Crops, SnowIce,
                        RunSuccessLandUse))
                frameList = pd.DataFrame(data=zippedList,
                                         columns=['ID', 'StationID', 'area', 'latitude', 'longitude', 'gravelius',
                                                  'perimeter', 'RunSuccessShape', 'elevation', 'slope', 'aspect',
                                                  'RunSuccessTerrain', 'forest', 'grass', 'wetland', 'water', 'urban',
                                                  'shrubs', 'crops', 'snowIce', 'RunSuccessLandUse'])
                frameList.to_csv('gauged_catchment_properties.csv')

            except Exception as e:
                print('{}: Failed, restarting iter {}'.format(e, attempt))
            else:
                break
