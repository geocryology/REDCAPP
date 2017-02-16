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
#============================ INTRODUCTION ====================================
#A example of extracting topographic factors, e.g. mrvbf, hypsometirc, from
#fine scale of DEM.
#==============================================================================

########################## HOW TO RUN THIS ####################################
#
# (1) Adapt the script below (settings) and run it
# (2) Make sure the input DEM is larger than required area
#
###############################################################################

import numpy as np

from os import path

# =============================== SETTING-UP =================================
dir_scr  = 'C:/Users/CaoBin/Documents/REDCAPP' #toposcale directory
dir_topo = 'C:/Users/CaoBin/Documents/REDCAPP/data' #DEM directory
dir_era  = 'C:/Users/CaoBin/Documents/REDCAPP/data' #ERA directory
file_out = 'C:/Users/CaoBin/Desktop/topo_station.csv' #output directory

execfile(path.join(dir_scr, 'topography.py'))

dem  = path.join(dir_topo, 'DEM_testArea.nc')
demResoultion = 3./3600 #indegree

stations=[{'name':'COV','lat':46.41801198,'lon':9.821232448, 'ele':3350.5},
          {'name':'SAM','lat':46.54873596,'lon':9.853765880, 'ele':2490.0},
          {'name':'SAM','lat':46.52639523,'lon':9.878944266, 'ele':1756.2}]

# ==================================== RUN ===================================
#topographic factors simultions
topo = topography(dem, demResoultion)
topo.describe()

out_xyz = np.asarray([[s['lat'],s['lon'],s['ele']] for s in stations])
mrvbf = topo.nmrvbf(out_xy = out_xyz[:,:2], initTf = 50)
hypso = topo.siteHypso(out_xy = out_xyz, bound = 30)
eleRange = topo.eleRange(out_xy = out_xyz[:,:2], bound = 30)
            
#export
topoEx = topoExport(mrvbf, hypso, eleRange, stations = stations)
topoEx.stationTopo(file_out)






