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
#A example of extracting mean of given date range from gridded data 
#based on a fine-scale DEM. Mean values are written into a netcdf file,
#variables are rounded to a precision of 3 decimal places.
#==============================================================================

from os import path
from datetime import datetime

##### HOW TO RUN THIS ##########################################################
#
# (1) Adapt the script below (settings, datarange) and run it
#
################################################################################

##settings
dir_scr  = '/Users/bincao/Google Drive/Script/era-Ta'#toposcale directory
dir_topo = '/Users/bincao/Documents/dem/GDEM'#DEM directory
dir_era  = '/Users/bincao/Documents/era'


execfile(path.join(dir_scr, 'downscaling.py'))

dem  = path.join(dir_topo, 'example_alps.nc')
geop = path.join(dir_era, 'alps_geop.nc')
sa   = path.join(dir_era, 'alps_sa_79_15.nc')
pl   = path.join(dir_era, 'alps_pl_79_15.nc')

#run
downscaling = downscaling(geop, sa, pl, dem)

variable = 'Temperature'
#Format of time range desired (yy, mon, day, hh, min)
daterange = {'beg' : datetime(2015, 01, 01, 00, 00),
             'end' : datetime(2015, 01, 01, 18, 00)}

file_out  = '/Users/bincao/Desktop/fineT.nc'

downscaling.extractSpatialDataNCF(daterange, variable, file_out)

plt.imshow(temp)
plt.colorbar()