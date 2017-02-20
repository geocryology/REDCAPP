# -*- coding: utf-8 -*-
# 
# REanalysis Downscaling Cold Air Pooling Parameterization (REDCAPP)
#
# === COPYRIGHT AND LICENCE ====================================================
#
# Copyright 2017 Bin Cao & Stephan Gruber
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==============================================================================
#
# Example file for producing spatial grid of mean air temperature
#
# (1) Get REDCAPP at https://github.com/geocryology/REDCAPP
#
# (2) Register to ECMWF (free) 
# https://apps.ecmwf.int/registration/
#
# (3) Follow the instructions for "Installing your API key" on
# https://software.ecmwf.int/wiki/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch
#
# (4) Adapt the script below (settings and location) and run the script
#
# (5) Use other example scripts to learn how to use data.
#
# ==============================================================================



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
from redcapp import downscaling

##### HOW TO RUN THIS ##########################################################
#
# (1) Adapt the script below (settings, datarange) and run it
#
################################################################################

# ================================ SETTING-UP =================================
dir_data = '/Users/stgruber/src/REDCAPP/data' # directory with DEM and ERA data
file_out = 'spatialT.nc'#output directory


# ================================== RUN ======================================
dem  = path.join(dir_data, 'DEM_testArea.nc')
geop = path.join(dir_data, 'GEOP_testArea.nc')
sa   = path.join(dir_data, 'ecmwf_erai_sa_m_151201_151231.nc')
pl   = path.join(dir_data, 'ecmwf_erai_pl_m_151201_151231.nc')
file_out = path.join(dir_data, file_out)

downscaling = downscaling(geop, sa, pl, dem)

variable = 'Temperature'

#       format of time range (yyyy, mm, dd, hh, mm)
daterange = {'beg' : datetime(2015, 12, 24, 00, 00),
             'end' : datetime(2015, 12, 24, 18, 00)}

downscaling.extractSpatialDataNCF(daterange, variable, file_out)