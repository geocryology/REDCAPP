#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (C) Copyright Bin Cao
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

# ================================ SETTING-UP =================================
dir_scr  = 'C:/OneDrive/GitHub/REDCAPP'#REDCAPP directory
dir_data = 'C:/OneDrive/GitHub/REDCAPP/Data'#DEM directory
file_out  = 'C:/OneDrive/GitHub/REDCAPP/Result/fineT.nc'#output directory

execfile(path.join(dir_scr, 'downscaling.py'))

dem  = path.join(dir_data, 'DEM_testArea.nc')
geop = path.join(dir_data, 'GEOP_testArea.nc')
sa   = path.join(dir_data, 'ecmwf_erai_sa_m_151201_151231.nc')
pl   = path.join(dir_data, 'ecmwf_erai_pl_m_151201_151231.nc')

# ================================== RUN ======================================
downscaling = downscaling(geop, sa, pl, dem)

variable = 'Temperature'
#Format of time range desired (yy, mon, day, hh, min)
daterange = {'beg' : datetime(2015, 12, 24, 00, 00),
             'end' : datetime(2015, 12, 24, 18, 00)}

downscaling.extractSpatialDataNCF(daterange, variable, file_out)