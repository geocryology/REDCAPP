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

##settings
dir_scr  = '/Users/bincao/Google Drive/Script/era-Ta'#script directory
dir_topo = '/Users/bincao/Documents/dem/GDEM'#DEM directory
dir_era  = '/Users/bincao/Documents/era'#era directory

execfile(path.join(dir_scr, 'downscaling.py'))

geop = path.join(dir_era, 'alps_geop.nc')#surface level
sa   = path.join(dir_era, 'alps_sa_79_15.nc')
pl   = path.join(dir_era, 'alps_pl_79_15.nc')

#run
downscaling = downscaling(geop, sa, pl)

variable = 'Temperature'

#Format of time range desired (yy, mon, day, hh, min)
daterange = {'beg' : datetime(2015, 01, 01, 00, 00),
             'end' : datetime(2015, 01, 10, 18, 00)}

file_out  = ['/Users/bincao/Desktop/pl_obs.csv',
             '/Users/bincao/Desktop/pl_dt.csv']
#Format of stations desired
#('name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber)             
stations = [{'name':'TAE','lat':47.47986561, 'lon':8.904870734, 'ele':585.80},
            {'name':'AAR','lat':47.38799123, 'lon':8.043881188, 'ele':454.87}]


downscaling.extractStationDataCSV(daterange, variable, stations, file_out)