#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (C) Copyright Bin Cao
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# For variable codes and units of ERA-Interim, see: 
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
#==============================================================================

import numpy as np
import netCDF4 as nc
import csv

from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter,generic_filter,convolve,minimum_filter,maximum_filter
from math import radians


class topography(object):
    """
    Return object for topography that has methods for deriving topographic 
    factors based on DEM.
    
    (a) MRVBF:
        Multiresolution index of valley bottom flatness (MRVBF) is derived by 
    Gallant and Dowling (2003), detailed could be found:
        http://onlinelibrary.wiley.com/doi/10.1029/2002WR001426/abstract
    The first slope threshold is set to 50% so that the mrvbf is smoother than
    original one and could represent the cold air pooling better.
    (b) Hypsometric Position:
        Hypsometric position is calculated as the radio of the number of cells
    with higher elevation than given site to the total number cells in the 
    surrounding region and ranges from 1 (deepest valley) to 0 (highest peak).
        
    Args:
         demFile: input dem file
         demResoultion: the resolution of input dem in degree
         
    Example:
        dem  = 'example_alps.nc'
        demResolution = 3./3600
        topo = topography(dem, demResoultion)
        topo.spatialTopo('/Users/bincao/Desktop/topo.nc')
    """

    def __init__(self, demFile, demResolution):
        self.dem = nc.Dataset(demFile, 'r')
        self.R = 6371000 # the mean radius (in meter) of Earth
        self.resolution = demResolution #units = degree
        self.lon = self.dem['lon'][:]
        self.lat = self.dem['lat'][:]
        self.shape = [len(self.lat), len(self.lon)]
        self.size = len(self.lon)*len(self.lat)
    
    def describe(self):
        """Returns summary information of DEM file
        Example:
            topo.describe()
        """
        minlon = np.min(self.lon)
        maxlon = np.max(self.lon)
        minlat = np.min(self.lat)
        maxlat = np.max(self.lat)

        print 'Ranges:  West: %s,  East: %s, South: %s, North: %s '  %(minlon, maxlon, minlat, maxlat)
        print 'Resolution: ' + str(self.resolution)
        print 'Dem shape [Lat, Lon]: ' + str(self.shape)
        print 'Dem size:' + str(self.size)

    def pixelLength(self, lat, L=1):
        """Return the cellsize for L step in meter. The function
        Args:
            lat: latitude in list
            L: interger, step number
            
        Returns:
            yL: length of cellsize in lontitude direction
            [x*cellsize for x in xL]: length of cellsize in latitude direction
            
        Example:
            yi,xi = topo.pixelLength(topo.lat, L=3)
        """

        y1 = self.resolution *np.pi*self.R/180#lon
        x1 = [np.cos(radians(elem)) *y1 for elem in lat]#lat at base scale
        if L <= 2:
            cellsize = 1
        else:
            cellsize = 3**(L-2)
        yL = y1 * cellsize
        xL = x1[(cellsize-1)/2:len(lat):cellsize]

        return [yL, [x*cellsize for x in xL]]#x1[np.arange((cellsize-1)/2,len(lat),cellsize)]

    def scale(self, x, t, p):
        """
        Scale the input value onto [0,1]

        Args:
            x: Input value (> 0)
            t: Threshold parameter
            p: Shape parameter, larger values give more abrupt transitions
            between 1 (for x << t) and 0 (for x >> t).
            
        Returns:
            scaled value
            """

        return 1/(1 + (x/t)**p)


    def smoothDEM(self, ele):
        """
        Returns a smoothed DEM in 2D Gaussian kernel array format.
        The smoothing is performed with the Arc/Info focal mean function
        using an 11 × 11 Gaussian smoothing kernel with an radius of three
        cell (to corresopnd with the factor of 3 resolution change).
        Detailed introduciton could be found from equation (11) in MRVBF
        paper.
        
        Args:
            ele: Input fine-scale of elevation
            
        Returns:
            Smoothed DEM
            
        Example:
            sdemL = topo.smoothDEM(topo.dem['elevation'][:])
        """
        
        return gaussian_filter(ele,np.sqrt(4.5), truncate = 11)
        
        

    def refine(self, L, coarseValue, out_xy = None, limitSize = 80000000):
        """
        Return refined coarseValue by linear interpolation.
        
        Args:
            L: Interger, step number
            coarseValue: Values of coarse scale need be refined
            out_xy: Sites need to be refined
            limitSize: Maximum cell size input once. This could avoid RAM crash
            
            
        Returns:
            Refined value
            
        Example:
            
        """
        
        scale = 3**(L-2)
        latIndex = range((scale-1)/2,len(self.lat),scale)
        lonIndex = range((scale-1)/2,len(self.lon),scale)
        f = RegularGridInterpolator((self.lat[latIndex][::-1],self.lon[lonIndex]),
                coarseValue[::-1,:], method = 'linear',bounds_error = False)

        if not (out_xy is None):
            return f(out_xy)

        elif self.size <= limitSize:
            lon, lat = np.meshgrid(self.lon, self.lat)
            gridBase = np.array([lat.reshape(lat.size), lon.reshape(lon.size)]).T
            fineValue = f(gridBase)
            return fineValue.reshape((len(self.lat),len(self.lon)))
        else:
            fineValue = np.zeros((len(self.lat),len(self.lon)))
            chunkN = self.size/limitSize
            latLegth = len(self.lat)/chunkN

            for pieceN in range(0, chunkN-1):
                startLat = latLegth*pieceN
                endLat = latLegth*(pieceN+1)
                lon, lat = np.meshgrid(self.lon, self.lat[startLat:endLat])
                gridBase = np.array([lat.reshape(lat.size), 
                                     lon.reshape(lon.size)]).T
                fineValue[startLat:endLat,:]=f(gridBase).reshape(latLegth,len(self.lon))

            lon, lat = np.meshgrid(self.lon, self.lat[latLegth*(chunkN-1):])
            gridBase = np.array([lat.reshape(lat.size), 
                                 lon.reshape(lon.size)]).T
            fineValue[latLegth*(chunkN-1):,:] = f(gridBase).reshape(len(self.lat[latLegth*(chunkN-1):]),len(self.lon))

            return fineValue


    def flatness(self, ele, Tf, out_xy = None,  L = 1, Pf = 4):
        """
        Return flatness for step L, flatness is defined as the inversion of 
        slope and called by scale() function.
        Slope is modeled by standard finite difference techniques, and 
        described as a percentage or 100 times the tangent of slope angle.
        
        The detailed describtion of calculating slope could be found:
            http://help.arcgis.com/En/Arcgisdesktop/10.0/Help/index.html#/How_Slope_works/009z000000vz000000/

        Args:
            ele: array_like, input elevation
            Tf: threshold, transform slope to flatness by scale() function
            out_xy: array_like, the interested sites. If not given (None), 
                flatnessof all the DEM cells would be simulated. Default is 
                None.
            L: interger, step number
            Pf: Shape parameter used to scale the value to range of 0-1.

        Results:
            Return flatness value at L step
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            #flatness for the finest step.
            F1 = topo.flatness(ele, out_xy = None, Tf = 50)
        """
        
        #degree to meter
        yi,xi = self.pixelLength(self.lat, L)
        xi = np.array(xi)
        ##slope
        kernelLon = np.array([[1,0,-1],[2,0,-2],[1,0,-1]])
        kernelLat = np.array([[1,2,1],[0,0,0],[-1,-2,-1]])
        dLon = convolve(ele, kernelLon)*100/(8*xi[:,None])
        dLat = convolve(ele, kernelLat)*100/(8*yi)
        slope = np.sqrt(dLon**2+dLat**2)
        del dLon, dLat

        #slope to flatness
        flatness = self.scale(slope, Tf, Pf)
        del slope
        #interpolate to the given sites
        if L <= 2:
            if not (out_xy is None):
                f = RegularGridInterpolator((self.lat[::-1], self.lon),
                    flatness[::-1,:], method = 'linear',bounds_error = False)
                flatness = f(out_xy)
            return flatness

        else:
            flatness = self.refine(L=L, coarseValue=flatness, out_xy=out_xy)
            return flatness


    def __percentile(self, x):
        """Return the rank of the element in surrouding cells.
        The function is called by lowness function"""
        return np.sum(x <= x[(x.size-1)/2])


    def lowness(self, ele, out_xy=None, L=1, Tl=0.4, Pl=3, lowRadius=13):
        """
        Returns lowness for given sites, which is measured as the radio of 
        number of points of lower elevation to the total number of points in 
        the surrounding region.

        Args:
            ele: array_like, input elevation
            out_xy: array_like, the interested sites
            L: interger, step number
            lowRadius: interger, the DEM cell size, half the number of cells
            used for the remainder of the steps
            L: interger, step number
            Tl: threshold, used to transform elevation percentile to lowness
            by scale() function
            Pl: shape parameter,used to transform elevation percentile to
            lowness by scale() function,
                  larger values give more abrupt transitions
                  
        Returns:
            Return lowness value at L step.
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            #lowness for the finest step.
            L1 = topo.lowness(ele, out_xy=None, lowRadius=7)
        """
        
        size = float(lowRadius)**2
        pctl = generic_filter(ele, self.__percentile,
                                size = lowRadius)
        pctl /= size
        lowness = self.scale(pctl, Tl , Pl)
        if L <= 2:
            if not (out_xy is None):
                f = RegularGridInterpolator((self.lat[::-1], self.lon[:]),
                    lowness[::-1,:], method = 'linear',bounds_error = False)
                lowness = f(out_xy)
            return lowness
        else:
            lowness = self.refine(L = L, coarseValue=lowness, out_xy=out_xy)
            return lowness


    def finestScale(self, out_xy = None, initTf = 50.0):
        """Return MRVBF2 and CF2, which are combined values determined
        from the first and second-scale steps.
        
        Args:
            out_xy: Array-like, sites need to be simulated. If not 
            given (None), all the DEM cells will be simulated.
            
        Returns:
            MRVBF2: Combined MRVBF erived from finest-scale and second
                steps
            CF2: Combined flatness derived from finest-scale and second
                steps
                
        Example:
            MRVBF2, CF2 = topo.finestScale()
        """
        
        #first step
        #ele = self.dem['elevation'][:]
        F1 = self.flatness(self.dem['elevation'][:], out_xy=out_xy, Tf=initTf)
        L1 = self.lowness(self.dem['elevation'][:], out_xy=out_xy, lowRadius=7)
        PVF1 = F1*L1
        VF1 = 1 - self.scale(PVF1, 0.3,4)
        #second step
        F2 = self.flatness(self.dem['elevation'][:], out_xy=out_xy,Tf=initTf/2)
        L2 = self.lowness(self.dem['elevation'][:], out_xy=out_xy, lowRadius=13)
        PVF2 = F2*L2
        VF2 = 1 - self.scale(PVF2, 0.3,4)
        #mrvbf2
        w2 = 1 - self.scale(VF2, 0.4, 6.68)
        MRVBF2 = w2*(1+VF2) + (1-w2)*VF1
        CF2 = F1*F2

        return MRVBF2, CF2

  
    def nmrvbf(self,out_xy = None, initTf = 50.0):
        """Returns smoothed mrvbf.
        
        Args:
        out_xy: Array-like, sites need to be simulated. If not 
            given (None), all the DEM cells will be simulated.
        initTf:
        
        Returns:
            mrvbf: smoothed and scaled mrvbf
        """
        mrvbf, cf = self.finestScale(initTf = initTf, out_xy = out_xy)
        sdemL = self.smoothDEM(self.dem['elevation'][:])#smoothed base resolution dem
        meanKernel = np.full((3,3), 1.0/(3*3))
        for L in range(3,9):
            #L = 4
            shape = sdemL.shape
            Tf = initTf / (2**(L-1))
            #aggreate sdem for L step
            latIndex = range(1, shape[0], 3)
            lonIndex = range(1, shape[1], 3)
            sdemL = convolve(sdemL, meanKernel)[latIndex, :]
            sdemL = sdemL[:, lonIndex]

            #MRVBF for L step
            FL = self.flatness(sdemL, Tf, out_xy = out_xy,L = L)
            LL = self.lowness(sdemL, out_xy = out_xy, L= L, lowRadius = 13)
            cf = cf*FL
            PVFL = cf*LL
            VFL = 1 - self.scale(PVFL, 0.3,4)
            wL = 1 - self.scale(VFL, 0.4, np.log10((L-0.5)/0.1)/np.log10(1.5))
            mrvbf = wL*(L-1+VFL) + (1-wL)*mrvbf

        return mrvbf
        

    def aggregation(self, dem, scaleFactor = 3):
        """"Return a aggregated DEM, in which the cell size is is increased by
        a factor of scaleFactor.
        
        Args:
            dem: Input orginal fine-scale dem
            scaleFactor: A factor used to aggregate dem.
            
        Returns:
            aggDEM: Aggregated coarse-scale dem.
            latIndex: Index of aggregated dem in respect of origincal latitude
            lonIndex: Index of aggregated dem in respect of origincal lontitude
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            aggDem, latIndex, lonIndex = topo.aggregation(ele, 5)
        """
        
        shape = dem.shape
        meanKernel = np.full((scaleFactor,scaleFactor), 
                             1.0/(scaleFactor*scaleFactor))
        latIndex = range((scaleFactor-1)/2, shape[0], scaleFactor)
        lonIndex = range((scaleFactor-1)/2, shape[1], scaleFactor)
        aggDem = convolve(dem, meanKernel)[latIndex, :]
        aggDem = aggDem[:, lonIndex]

        return [aggDem, latIndex, lonIndex]


    def aroundArea(self, centerSite, size = 30):
        """Return the surrouding values of givens site.
        Args:
            centerSite: The center cell [lat, lon, ele] of the surrounding 
            area. All the cells with the distance <= size/2 will be clipped.
            size: Diameter of surrounding area.
        
        Returns:
            Array like clipped value with given size.
            
        Example:
            centerSite = [46.749374, 9.5506067, 1556.0]
            eleSubset = topo.aroundArea(centerSite)
        """
        
        #area size
        yi = self.pixelLength(self.lat)[0]
        lowRadius = int(size*1000/(2*yi))#radius
        #nearest cell index
        latPosition = np.abs(self.dem.variables['lat'][:]-centerSite[0]).argmin()
        lonPosition = np.abs(self.dem.variables['lon'][:]-centerSite[1]).argmin()
        #area cell index
        latli = latPosition - lowRadius
        latui = latPosition + lowRadius
        lonli = lonPosition - lowRadius
        lonui = lonPosition + lowRadius
        #nearest area with bound size
        eleSubset = self.dem.variables['elevation'][latli:latui, lonli:lonui]
        return eleSubset

    def siteHypso(self, out_xy, bound = 30):
        """Returns the hypsometric position in the surrouding area.
        
        Args:
            out_xyz: Site to simulate hypsometric position.
            bound: Diameter of surrounding area in km.
            
        Returns:
            lowness: Hyposmetric position in the surrouding area. lowness
            ranges from 1 (deepest valley) to 0 (highest peak).
            
        Example:
            out_xy = np.array([[46.749374, 9.5506067, 1556.0],
                               [46.665974, 9.6340027, 939.0]])
            lowness = topo.siteHypso(out_xy)
        """
        
        lowness = []
        for i, site in enumerate(out_xy):
            #print(site)
            #site = [ 46.4820675,    6.98736116]
            ind_ele = self.aroundArea(site, size = bound)#aroud ele
            pctl = np.sum(ind_ele >= out_xy[i,2])/float(ind_ele.size)
            lowness.append(pctl)
        return lowness


    def hypso(self, bound = 30, out_xy = None):
        """Return hypsometric position based on fine scale of dem.

        Args:
            bound: Diameter of surrounding size
            
        Returns:
            lowness: Hypsometric position of all dem cells derived from
            fine scale of dem.
            
        Examples:
            hypso = topo.hypso(bound = 30, out_xy = None)
        """
        
        yi = self.pixelLength(self.lat)[0]
        lowRadius = int((bound*1000/yi))
        pctl = generic_filter(self.dem['elevation'][:], 
                              self.__percentile, size = lowRadius)
        pctl = pctl/(float(lowRadius)**2)
        lowness = 1-pctl
        del pctl

        if not (out_xy is None):
            lowInterp = RegularGridInterpolator((self.lat[::-1],self.lon),
                        lowness[::-1], method = 'linear', bounds_error = False)
            lowness = lowInterp(out_xy)
        return lowness


    def coarseHypso(self, bound = 30, out_xy = None):
        """Return hypsometric position based on coarse scale of dem.

        Args:
            bound: Diameter of surrounding size
            
        Returns:
            lowness: Hypsometric position of all dem cells derived from
            fine scale of dem.
            
        Examples:
            hypso = topo.hypso(bound = 30, out_xy = None)
        """

        yi = self.pixelLength(self.lat)[0]
        scaleFactor = int(np.ceil(500//yi) // 2 * 2 + 1)
        aggDem = self.aggregation(self.dem['elevation'][:], scaleFactor)
        #lowness
        lowRadius = 61#bound*1000/500
        pctl = generic_filter(aggDem[0], 
                              self.__percentile, 
                              size = lowRadius)/(float(lowRadius)**2)
        lowness = 1 - pctl

        #refine
        lowInterp = RegularGridInterpolator((self.lat[aggDem[1]][::-1],
                                             self.lon[aggDem[2]]),
            lowness[::-1], method = 'linear', bounds_error = False)


        if not (out_xy is None):
            lowness = lowInterp(out_xy)

        else:
            lon, lat = np.meshgrid(self.lon, self.lat)
            gridBase = np.array([lat.reshape(lat.size), 
                                 lon.reshape(lon.size)]).T
            lowness = lowInterp(gridBase).reshape((len(self.lat), len(self.lon)))

        return lowness#,pctl

    
    def eleRange(self, bound = 30, out_xy = None):
        """Return elevation range for each dem cell within a area.
        
        Args:
            bound: Diameter of surrounding size
            
        Returns:
            rangeE: Array like, range of elevation range.
        
        
        Example:
            eleRange = topo.eleRange(out_xy = None, bound = 30)
        """
        lowRadius = bound*1000/(self.pixelLength(self.lat)[0])
        minEle = minimum_filter(self.dem['elevation'][:], size = lowRadius)
        maxEle = maximum_filter(self.dem['elevation'][:], size = lowRadius)
        rangeE = maxEle - minEle
        del minEle, maxEle
        if not (out_xy is None):
            rangeInterp = RegularGridInterpolator((self.lat[::-1],self.lon),
            rangeE[::-1], method = 'linear', bounds_error = False)
            rangeE = rangeInterp(out_xy)
        return rangeE
        

class topoExport(object):
    """"Returns a file of netCDF4 or csv contains all topographic factors 
    needed by REDCAPP. The edges of input mrvbf without data will be dropped.
    
    Args:
        demFile
        mrvbf: Array-like multiresolution index of valley bottoms flatness with NA at 
               edges. Array-like
        hypso: Array-like hypsometric position. 
        eleRange: Array-like elevation range
    
    Returns:
        A file of netCDF4 (csv) for spatial (station) topographoic information.
        
    Examples:
        #spatial topographic factors
        topoEx = topoExport(mrvbf, hypso, eleRange, dem = demFile)
        topoEx.spatialTopo()
        
        #station topographic factors
        topoEx = topoExport(mrvbf, hypso, eleRange, stations = stations)
        topoEx.stationTopo()       
    """
    
    def __init__(self, mrvbf, hypso, eleRange, demFile = None, stations = None):
        if not (demFile is None):
            self.dem = nc.Dataset(demFile, 'r')
        if not (stations is None):
            self.stations = stations
        self.mrvbf = mrvbf
        self.hypso = hypso
        self.eleRange = eleRange
        
        
    def edgeClip(self):
        """Drops the egde of mrvbf without data. The function is called
        by spatialTopo().
        """
        #the corner with values
        shape = self.mrvbf.shape
        center = [i/2 for i in shape]
        left = np.min(np.where(np.isfinite(self.mrvbf[center[0],:])))#left
        right = np.max(np.where(np.isfinite(self.mrvbf[center[0],:])))+1#right
        
        upper = np.min(np.where(np.isfinite(self.mrvbf[:,center[1]])))#upper
        low = np.max(np.where(np.isfinite(self.mrvbf[:,center[1]])))+1#low
        
        mrvbf = self.mrvbf[upper:low, left:right]
        hypso = self.hypso[upper:low, left:right]
        eleRange = self.eleRange[upper:low, left:right]
        lons = self.dem.variables['lon'][left:right]
        lats = self.dem.variables['lat'][upper:low]
        
        return mrvbf,hypso,eleRange,lons,lats
    
    def spatialTopo(self, file_out):
        """Export a netCDF4 file contains all topographic information.
        
        Args:
            file_out = 'C:/Users/CaoBin/Desktop/topo_testArea.nc'
            
        Example:
            topoEx.spatialTopo(file_out)
        """
        
        mrvbf,hypso,eleRange,lons,lats = self.edgeClip()
        mrvbf/=8.0
        
        #create nc file
        nc_root = nc.Dataset(file_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', mrvbf.shape[0])
        nc_root.createDimension('lon', mrvbf.shape[1])
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes  = nc_root.createVariable('lat', 'f4', ('lat'))
        Hypso  = nc_root.createVariable('hypso','f4',('lat','lon'),zlib=True)
        Mrvbf  = nc_root.createVariable('mrvbf','f4',('lat','lon'),zlib=True)
        RangeE = nc_root.createVariable('range','f4',('lat','lon'),zlib=True)
        
        longitudes.setncatts({'long_name': u"longitude"})
        latitudes.setncatts({'long_name': u"latitude"})
        Mrvbf.setncatts({'long_name': 
                         "normalized multiresolution index of valley bottom flatness"})
        Hypso.setncatts({'long_name': u"hyposmetric position"})
        RangeE.setncatts({'long_name': u"elevation range in prescirbef neighbourhood"})
        
        #assign variables
        longitudes[:] = lons
        latitudes[:] = lats
        Hypso[:] = hypso
        Mrvbf[:] = mrvbf
        RangeE[:] = eleRange
        
        #attribute
        nc_root.description = "DEM-derived topographic factors"
        longitudes.units = 'degree_east (demical)'
        latitudes.units = 'degree_north (demical)'
        
        nc_root.close()
        
    def stationTopo(self, file_out, initTf=50, bound=30):
        """Export a csv file contains all topographic information for given
        stations.
        
        Args:
            file_out = 'C:/Users/CaoBin/Desktop/topo_stations.nc'
            
        Example:
            topoEx = topoExport(mrvbf, hypso, eleRange, stations = stations)
            topoEx.stationTopo(file_out)
        """
        
        names = [s['name'] for s in self.stations]#station names
        #topographic values [hypso, mrvbf, elevationRange]
        values = np.array([self.hypso, self.mrvbf/8.0, self.eleRange]).T
        
        #write CSV
        with open(file_out, 'wb') as output_file:
           writer = csv.writer(output_file)
           writer.writerow(['station','hypso','mrvbf','eleR'])
           for n in range(len(names)):
               row = ['%.2f' % elem for elem in values[n]]
               row.insert(0,names[n])
               writer.writerow(row)
        output_file.close()
        
    
class landSurCorrectionFac(object):
    """
    Returns land surface correction factor based on input topographyic
    information.
    
    Args:
        alpha: Adjust constant value
        tau: A factor relating to fractional influence of surface effects on
            air temperature.
        kappa: A factor relating to cold air pooling on air temperature.
            
    Returns:
        lscf: Land surafce correction factor
        
    Example:
        lscf = lscf(alpha= 0.64, betta= 1.34, sigma= 487)
        LSCF = lscf.correctionFactor()
    
    """
    
    def __init__(self, alpha= 0.64, betta= 1.34, sigma= 487):
        self.sigma = sigma
        self.betta = betta
        self.alpha = alpha
    
    def scale(self, eleRange):
        return(np.exp(-eleRange/self.sigma))
        
    def correctionFactor(self, hypso, mrvbf, eleRange):
        s = self.scale(eleRange)
        p = mrvbf*(1-s)
        f = hypso* (1- s) + s
        lscf = self.betta*f + self.alpha*p
        return(lscf)
