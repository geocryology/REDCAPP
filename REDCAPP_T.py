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
# =============================================================================

import matplotlib.pyplot as plt

from os import path
from datetime import datetime

# ============================== SETTING-UP ===================================
dir_scr  = 'C:/OneDrive/GitHub/REDCAPP'# toposcale directory
dir_data  = 'C:/OneDrive/GitHub/REDCAPP/data'# Data directory
file_out = 'C:/Users/BinCao/Desktop/topo_testArea.nc'# output directory


execfile(path.join(dir_scr, 'topography.py'))# topography script
execfile(path.join(dir_scr, 'downscaling.py'))# temperature script


#Input datasets
dem  = path.join(dir_data, 'DEM_testArea.nc')# fine-scale dem
geop = path.join(dir_data, 'GEOP_testArea.nc')# coarse-scale dem
sa   = path.join(dir_data, 'ecmwf_erai_sa_m_151201_151231.nc')# 2-m temperature
pl   = path.join(dir_data, 'ecmwf_erai_pl_m_151201_151231.nc')# pressure level temperature
demResoultion = 3./3600#in degree

# ================================== RUN ======================================
###upper-air temperature & land-surface effects
Downscaling = downscaling(geop, sa, pl, dem)

variable = 'Temperature'
#Format of time range desired (yy, mon, day, hh, min)
daterange = {'beg' : datetime(2015, 12, 24, 00, 00),
             'end' : datetime(2015, 12, 24, 18, 00)}

out_xyz_dem, lats, lons, shape = Downscaling.demGrid()
out_xyz_sur = Downscaling.surGrid(lats, lons, None)
pl,dt = Downscaling.spatialMean(variable, daterange, out_xyz_sur, 
                                 out_xyz_dem, shape)

##Land surface correction factors
#topographic factors
topo = topography(dem,demResoultion)
topo.describe()#dem description

mrvbf= topo.nmrvbf(out_xy = None, initTf = 50.0)
hypso= topo.coarseHypso(out_xy = None, bound = 30)
eleRange= topo.eleRange(out_xy = None, bound = 30)

topoEx = topoExport(mrvbf, hypso, eleRange, demFile = dem)
mrvbf, hypso, eleRange, lons, lats = topoEx.edgeClip()

#dem-derived lscf
lscf = landSurCorrectionFac(alpha= 0.64,betta= 1.34, sigma= 487)
LSCF = lscf.correctionFactor(hypso, mrvbf/8.0, eleRange)

# =================================== PLOTS ===================================

