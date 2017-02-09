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

from os import path

#setting
dir_scr  = '/Users/bincao/Google Drive/Script/era-Ta'#toposcale directory
dir_topo = '/Users/bincao/Documents/dem/GDEM'#DEM directory
dir_era  = '/Users/bincao/Documents/Data/era'
file_out = '/Users/bincao/Desktop/topo.csv'

execfile(path.join(dir_scr, 'topography.py'))

dem  = path.join(dir_topo, 'example_alps.nc')
demResoultion = 3./3600#indegree
topo = topography(dem, demResoultion)
topo.describe()

stations = [{'name':'COV','lat':46.41801198,'lon':9.821232448, 'ele':3350.520995},
            {'name':'SAM','lat':46.52639523,'lon':9.878944266, 'ele':1756.194106}]

            
#run
topo.stationTopo(stations, file_out)





out_xy = np.asarray([[s['lat'],s['lon'],s['ele']] for s in stations])

mrvbf = topo.mrvbf(out_xy = out_xy[:,:2], initTf = 50)
hypso = topo.siteHypso(out_xy = out_xy, bound = 30)
eleRange = topo.eleRange(out_xy = out_xy[:,:2], bound = 30)