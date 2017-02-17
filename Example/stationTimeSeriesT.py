#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (C) Copyright Bin Cao & Stephan Gruber
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# For variable codes and units of ERA-Interim, see: 
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
#==============================================================================
#A example of extracting time series from gridded data based on a list of 
#dictionaries describing stations. Time serie(s) are written into csv files,
#variables are rounded to a precision of 3 decimal places.
#==============================================================================

from os import path
from datetime import datetime

##### HOW TO RUN THIS #########################################################
#
# (1) Adapt the script below (settings, datarange) and run it
#
###############################################################################

# ================================ SETTING-UP =================================
dir_scr  = 'C:/OneDrive/GitHub/REDCAPP'# REDCAPP directory
dir_data = 'C:/OneDrive/GitHub/REDCAPP/Data'# Data directory
file_out  = ['C:/OneDrive/GitHub/REDCAPP/Result/upperAirT.csv',
             'C:/OneDrive/GitHub/REDCAPP/Result/coarseLandSurEffect.csv']# Fileout directory

execfile(path.join(dir_scr, 'downscaling.py'))

geop = path.join(dir_data,  'GEOP_testArea.nc')
sa   = path.join(dir_data,  'ecmwf_erai_sa_m_151201_151231.nc')
pl   = path.join(dir_data,  'ecmwf_erai_pl_m_151201_151231.nc')

# ================================== RUN ======================================
downscaling = downscaling(geop, sa, pl)

variable = 'Temperature'

#Format of time range desired (yy, mon, day, hh, min)
daterange = {'beg' : datetime(2015, 12, 01, 00, 00),
             'end' : datetime(2015, 12, 31, 18, 00)}

#Format of stations desired
#['name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber]             
stations=[{'name':'COV','lat': 46.41801198, 'lon': 9.821232448, 'ele': 3350.5},
          {'name':'SAM','lat': 46.52639523, 'lon': 9.878944266, 'ele': 1756.2}]


downscaling.extractStationDataCSV(daterange, variable, stations, file_out)