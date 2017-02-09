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
#A example of extracting topographic factors, e.g. mrvbf, hypsometirc, from
#fine scale of DEM.
#==============================================================================

##### HOW TO RUN THIS #########################################################
#
# (1) Adapt the script below (settings) and run it
# (2) Make sure the input DEM is larger than required area
#
###############################################################################


#setting
dir_scr  = 'C:/OneDrive/Script/era-Ta'#toposcale directory
dir_topo = 'D:/Data/dem/GDEM'#DEM directory
dir_era  = 'D:/Data/era'
file_out = 'C:/Users/CaoBin/Desktop/topo_alps.nc'

execfile(path.join(dir_scr, 'topography.py'))

dem  = path.join(dir_topo, 'alps_arc3.nc')
demResoultion = 3./3600#indegree
topo = topography(dem,demResoultion)
topo.describe()


#run
topo.spatialTopo(file_out)


mrvbf3 = topo.mrvbf(out_xy = None, initTf = 50.0)




out_xy = np.array([[46.41801198, 9.821232448, 3350.520995],
                   [46.52639523, 9.878944266, 1756.194106]])

topo.mrvbf(out_xy = out_xy[:,:2], initTf = 50.0)