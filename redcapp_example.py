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
# Example file containing: (a) data download,
#                          (b) obtaining spatialized mean air temperature
#                          (c) obtaining air temperature time series at point 
#
# (1) Get REDCAPP at https://github.com/geocryology/REDCAPP
#
# (2) Make sure that the directory containing thie file (redcapp_example.py) is 
#     contained in your PYTHONPATH.
#
# (3) Register to ECMWF (free) https://apps.ecmwf.int/registration/
#
# (4) Follow the instructions for "Installing your API key" on
# https://software.ecmwf.int/wiki/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch
#
# (5) Adapt the script below (settings and location) 
#
# (6) Run the script
#
# (7) Explore the results. Use a netcdf viewer to plot maps and time series.
#     Panoply (https://www.giss.nasa.gov/tools/panoply) is a good one.
#
# (8) Customise the code and use it for your project.
#
# (9) Please let us know how things work. We hope this is useful for you.
#
# ==============================================================================

from redcapp import redcapp_get, eraData
from datetime import datetime
from os import path


# directory containing all raw data and output data
dir_data = '/Users/stgruber/Desktop/data'

# Digital ELevation Model in ASCIIGRID format and lat/lon WGS84 grid
dem_file = 'lalala'

# output file names
dem_out  = 'DEM_fine-scale.nc' # high-resolution topography file 
spat_out = 'spatialT.nc'       # spatialised mean air temperature
#TODO: make this only one file per variable, only give name and the n 
tser_out  = ['C:/OneDrive/GitHub/REDCAPP/Result/upperAirT.csv',
             'C:/OneDrive/GitHub/REDCAPP/Result/coarseLandSurEffect.csv']# Fileout directory



#location: alps
date  = {'beg' : datetime(2015,1,1),
         'end' : datetime(2015,1,10)}
         
area  = {'north' : 46.65,
         'south' : 46.35,
         'west'  : 9.70, #positive is westwards of Greenwich
         'east'  : 9.95} #positive is westwards of Greenwich
         
elevation = {'min' : 0, 
             'max' : 4500}

#Format of stations for which to provide time series
#['name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber]             
stations=[{'name':'COV','lat': 46.41801198, 'lon': 9.821232448, 'ele': 3350.5},
          {'name':'SAM','lat': 46.52639523, 'lon': 9.878944266, 'ele': 1756.2}]

# make file names
dem_ncdf  = path.join(dir_data, 'DEM_testArea.nc')
geop      = path.join(dir_data, 'ecmwf_erai_to.nc')
spat_out  = path.join(dir_data, spat_out)

# === DOWNLOAD =================================================================
                                
rg = redcapp_get(date, area, elevation, dir_data, 5) 
rg.retrieve()

eraDownload = eraData()
eraDownload.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_sa_*'),1)
eraDownload.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_pl_*'),1)


# ==== IMPORT DEM ==============================================================

# TODO: Make function that allows user to import ASCII files. Provide example 
# data as ASCII file. In the import DEM function, have a flag:
# "geomorphometry = true" so that mrvbf, hypso, and range are already calculated
# and added to the dem netcdf.



# ==== SPATIALIZED MEAN TEMPERATURE ============================================

# TODO: make this a function that finds teh right files in the directory based on
# 'ecmwf_erai_sa_*' and 'ecmwf_erai_pl_*'. Otherwise the user has to adjust 
#  the file names manually
sa   = path.join(dir_data, 'ecmwf_erai_sa_m_151201_151231.nc')
pl   = path.join(dir_data, 'ecmwf_erai_pl_m_151201_151231.nc')

#TODO: DOES THIS ACTUALLY CONTAIN THE LSCF?

downscaling = downscaling(geop, sa, pl, dem_ncdf)
variable = 'Temperature'
downscaling.extractSpatialDataNCF(date, variable, spat_out)


# ==== AIR TEMPERATURE TIME SERIES =============================================

downscaling.extractStationDataCSV(date, variable, stations, tser_out)

