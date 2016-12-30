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

import netCDF4 as nc
import numpy as np
import csv


from scipy.interpolate import RegularGridInterpolator
from bisect import bisect_left



class downscaling(object):
    """Rerurn object for downscaling that has methods for interpolationg
    upper-air temperature and surface influences at surface level
    based on disaggregating coarse-grid reanalysis and dem.
    
    Args:
        dem: A required fine-scale dem in netcdf format
    
    Example:
        dem  = 'example_alps.nc'
        geop = 'alps_geop.nc'
        sa   = 'alps_sa_79_15.nc'
        pl   = 'alps_pl_79_15.nc'
        
        downscaling = downscaling(dem, geop, sa, pl)
        
    """
    
    def __init__(self, geop, sa, pl, dem = None):
        self.g    = 9.80665 #Gravitational acceleration [m/s2]
        self.geop = nc.Dataset(geop)
        self.sa   = nc.Dataset(sa)
        self.pl   = nc.Dataset(pl)
        if not (dem is None):
            self.dem  = nc.Dataset(dem)
        
        
    def demGrid(self, stations = None):     
        """Return metadata of given stations or dem. 
        Format of stations desired [lat, lon, geop].
        
        Args: 
            stations: A list of dictionaries describing stations. If not given,
                metadata derived from given dem.
            
        Returns:
            out_xyz_dem: Metadata [lat, lon, geop] of input dem or sites
            lons: Longitude of input sites
            lats: Latitude of input sites
            shape: Shape of input dem or sites  
        """
        
        if not (stations is None):
            lats = [s['lat'] for s in stations]
            lons = [s['lon'] for s in stations]
            names = [s['name'] for s in stations]
            siteLocation = np.asarray([[s['lat'],s['lon'],
                                        s['ele']*self.g] for s in stations])
            shape = siteLocation.shape
            return siteLocation, lats, lons, shape, names        

        #out_xyz based on dem
        lons = self.dem.variables['lon'][:]
        lats = self.dem.variables['lat'][:]
        geop = self.dem.variables['elevation'][:]*self.g
        shape = geop.shape
        
        
        lons, lats = np.meshgrid(lons, lats)
        lons = lons.reshape(lons.size)
        lats = lats.reshape(lats.size)
        geop = geop.reshape(geop.size)
        
        out_xyz_dem = np.array([lats, lons, geop]).T

        return out_xyz_dem, lats, lons, shape
        
        
    def geoGrid(self):     
        """Returns original (NO interpolation) coarse gird metadata
        in format of [lat, lon, geop], based on geopotential file
        """
        longitude = self.geop['longitude'][:]
        latitude = self.geop['latitude'][:]
        lons, lats = np.meshgrid(longitude, latitude)
        lons = lons.reshape(lons.size)
        lats = lats.reshape(lats.size)
        geop = self.geop['z'][0,:,:]#geopotential
        geop = geop.reshape(geop.size)

        out_xyz_ori = np.array([lats, lons, geop]).T
        
        return out_xyz_ori
        
    
    def surGrid(self, lats, lons, stations):
        
        """Return interpolated surface geopotetial.
        
        Args: 
            lats: Latitude of intersted sites
            lons: Longitude of intersted sites
            stations: A list of dictionaries describing stations. If not given,
                metadata derived from given dem.
        
        Returns:
            out_xyz_sur: Fine-scale surface level metadata [lat, lon, geop],
            where latitude and lontitude are ontained from dem,
            while geop is interpolated from coarse geopotential file
            
        Example:
            downscaling = downscaling(dem, geop, sa, pl)
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_ori = downscaling.geoGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, out_xyz_dem[:,:2])
        """
        
        longitude = self.geop['longitude'][:]
        latitude = self.geop['latitude'][::-1]
        in_v = self.geop['z'][0,::-1,:]#geopotential
        fz = RegularGridInterpolator((latitude,longitude), in_v, 'linear')
        out_xy = np.array([lats, lons]).T

        if not (stations is None):
            lats = [s['lat'] for s in stations]
            lons = [s['lon'] for s in stations]
            out_xy = np.asarray([[s['lat'],s['lon']] for s in stations])
        
        z_interp = fz(out_xy)
        out_xyz_sur = np.array([lats, lons, z_interp]).T
        
        return out_xyz_sur

 

    def surTa(self, ind_time, out_xyz_sur):
        """Return interpolated 2-metre temperature.
        
        Args:
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            out_xyz_sur: 
        
        Returns:
            t_sa: interpolation fine-scale surface air temperatue based on 
            grid 2-meter temperature
            
        Example:
            dem  = 'example_alps.nc'
            geop = 'alps_geop.nc'
            sa   = 'alps_sa_79_15.nc'
            pl   = 'alps_pl_79_15.nc'
        
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)


            surTa = downscaling.surTa(0, out_xyz_sur)
        """
            
        in_v = self.sa['2 metre temperature'][ind_time,0,:,:]#geopotential
        in_v -= 273.15
        lat = self.sa.variables['lat'][:]
        lon = self.sa.variables['lon'][:]

        f_sa = RegularGridInterpolator((lat,lon), in_v, 'linear')
        t_sa = f_sa(out_xyz_sur[:,:2]) 

        return t_sa


    def gridValue(self, variable, ind_time):
        """Return original grid temperatures and geopotential of differnet
        pressure levels. The function are called by inLevelInterp() to
        get the input ERA-Interim values.
        
        Args: 
            variable: Given interpolated climate variable
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            
        Returns:
            gridT: Grid temperatures of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridLon: Grid longitude of pressure level variables
            gridLat: Grid latitude of pressure level variables
        
        Example:
            gridT,gridZ,gridLat,gridLon=downscaling.gridValue('Temperature',0)
            
        """
        
        gridT = self.pl.variables[variable][ind_time,:,:,:]
        gridZ = self.pl.variables['Geopotential'][ind_time,:,:,:]
        #x and y
        
        gridLat = self.pl['lat'][:]
        gridLon = self.pl['lon'][:]
        

        return gridT,gridZ,gridLat,gridLon


    def inLevelInterp(self,gridT, gridZ, gridLat, gridLon, out_xyz):
        """This is a 2D interpolatation, and returns interpolated temperatures
        of different pressure levels.
        
        Args:
            gridT: Grid temperatures of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridLat: Grid longitude of pressure level variables
            gridLon: Grid latitude of pressure level variables
            out_xyz: Given sites, which will be interpolated.
            
        Returns:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
        
        Examples:
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            #interpolate 2-meter temperature
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
        """
        
        shape = gridT.shape
        #create array to hold interpolation resultes
        t_interp = np.zeros([shape[0], len(out_xyz)])
        z_interp = np.zeros([shape[0], len(out_xyz)])

        #temperatue and elevation interpolation 2d
        for i in range(shape[0]):
            ft = RegularGridInterpolator((gridLat,gridLon), 
                                          gridT[i,:,:], 'linear')
            fz = RegularGridInterpolator((gridLat,gridLon), 
                                          gridZ[i,:,:], 'linear')
            t_interp[i,:] = ft(out_xyz[:,:2])#temperature
            z_interp[i,:] = fz(out_xyz[:,:2])#elevation

        t_interp -= 273.15

        return t_interp[::-1,:], z_interp[::-1,:]
        
    

        
    def fast1d(self, t_interp, z_interp, out_xyz):
        """This is a 1D interpoation. The function return interpolated 
        upper air temperature at the given sites by 
        interpolation between different pressure levels.
        
        Args:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            out_xyz: Given sites with elevation, which will be interpolated.
            
        Returns:
            dG:upper-air temperature at given sites
                
        Example: 
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
            
            #upper air temperature at the coarse and fine scale of elevation
            pl_sa = fast1d(t_interp, z_interp, out_xyz_sur)
            pl_obs = fast1d(t_interp, z_interp, out_xyz_dem)
        """
    
        ele = out_xyz[:,2]
        size = np.arange(out_xyz.shape[0])
        n = [bisect_left(z_interp[:,i], ele[i]) for i in size]
        n = [x+1 if x == 0 else x for x in n]
        
        lowN = [l-1 for l in n]
        
        upperT = t_interp[n,size]
        upperZ = z_interp[n,size]
        dG  = upperT-t_interp[lowN,size]#<0
        dG /= upperZ-z_interp[lowN,size]#<0
        dG *= out_xyz[:,2] - upperZ#>0
        dG += upperT
             
        return dG
    
          
        
    def interpAll(self, variable, ind_time, out_xyz_sur, out_xyz_obs):
        """Returns all needed interpolated temperatures at given time.
        
        Args:
            variable: Interpolated climated variable
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            out_xyz_sur: Interploated sites[lat, lon, geop], in which the 
            geopotential are interpolated from ERA-Interim geopotential file.
            out_xyz_obs: Interploated sites[lat, lon, geop], in which the 
            geopotential are gotten from DEM.
            
        Returns:
            pl_obs: Interpolated upper-air temperature at given sites
                with dem level
            pl_sur: Interpolated upper-air temperature at given sites
                with surface level
            
            
        Example:
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)
            pl_obs,pl_sur,t_sa = downscaling.interpAll(variable, ind_time, 
                                                       out_xyz_sur, 
                                                       out_xyz_dem)
            
            
        """
        gridT,gridZ,gridLat,gridLon = self.gridValue(variable, ind_time)
        t_interp,z_interp = self.inLevelInterp(gridT,gridZ,gridLat,gridLon,
                                               out_xyz_obs)
        pl_sur = self.fast1d(t_interp, z_interp, out_xyz_sur)
        pl_obs = self.fast1d(t_interp, z_interp, out_xyz_obs)
        t_sa = self.surTa(ind_time, out_xyz_sur)
        dt = pl_sur-t_sa
        
        return pl_obs, dt

    def spatialMean(self, variable, daterange, out_xyz_sur, out_xyz_obs, shape):
        """Return the MEAN upper-air temperature and 
        land surface influence during given date range and at given  area
        
        Args:
            variable: Climated variable to interpolate
            daterange: Date range to interpolate
            out_xyz_sur: Surface level sites to interpolate 
            out_xyz_obs: Topography sites to interpolate
                
        Returns:
            pl: Interpolated MEAN free-atmosphere during given date range and
            at given area.
            dt: Interpolated MEAN surface level land surface influences during 
            given date range and at given area.
            
        Example:
            
            dem  = 'example_alps.nc'
            geop = 'alps_geop.nc'
            sa   = 'alps_sa_79_15.nc'
        
            downscaling = downscaling(dem, geop, sa, pl)
            
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)
            pl,dt = downscaling.spatialMean(variable, daterange, 
                                            out_xyz_sur, out_xyz_dem, shape)
            
        """
        
        #obtain time range
        date_vec = nc.num2date(self.pl.variables['time'][:],
                                units = "seconds since 1970-1-1",
                                calendar='standard')


        #index of time steps to interpolate
        mask  = date_vec >= daterange.get('beg')
        mask *= date_vec <= daterange.get('end')
        ind_time_vec = np.arange(len(mask))[mask]
        out_time = date_vec[mask]

        sum_pl_obs = 0
        sum_dt = 0

        print("\nConducting downscaling now, have a cup of coffee please\n")

        for ind_out, ind_time in enumerate(ind_time_vec):
            print(out_time[ind_out])
            pl_obs,dt = self.interpAll(variable,ind_time, out_xyz_sur, 
                                                out_xyz_obs)
            sum_pl_obs += pl_obs
            sum_dt += dt

        dt = sum_dt/out_time.size
        pl = sum_pl_obs/out_time.size
        
        dt = dt.reshape(shape)
        pl = pl.reshape(shape)

        return pl,dt
       
    def stationTimeSeries(self, variable, daterange, out_xyz_sur, out_xyz_obs):
        """Return upper-air temperature and land surface influence
        at given time steps
        
        Args:
            variable: Climated variable to interpolate
            daterange: Date range to interpolate
            out_xyz_sur: Surface level sites to interpolate 
            out_xyz_obs: Topography sites to interpolate
            
        Returns:
            out_valu: 
            out_time: 
            
        Example:
            downscaling = downscaling(dem, geop, sa, pl)

            site = np.array([[47.38799123, 8.043881188, 454.87048*9.80665]])
            out_xyz_dem, lats, lons, shape = downscaling.demGrid(site)
            out_xyz_sur = downscaling.surGrid(lats, lons, site)
            
            out_valu,out_time=downscaling.stationTimeSeries(variable, 
                                                            daterange, 
                                                            out_xyz_sur,
                                                            out_xyz_dem)
        """
        
        #obtain time range
        date_vec = nc.num2date(self.pl.variables['time'][:],
                                units = "seconds since 1970-1-1",
                                calendar='standard')

        #index of time steps to interpolate
        mask  = date_vec >= daterange.get('beg')
        mask *= date_vec <= daterange.get('end')
        ind_time_vec = np.arange(len(mask))[mask]
        out_time = date_vec[mask]

        out_valu = np.zeros((2, ind_time_vec.size, 
                             out_xyz_obs.shape[0]))#pl_obs, dt

        print("\nConducting downscaling now, have a cup of coffee please\n")

        for ind_out, ind_time in enumerate(ind_time_vec):
            print(out_time[ind_out])
            out_valu[:,ind_out,:] = self.interpAll(variable, ind_time,
                                                   out_xyz_sur, out_xyz_obs)

        return out_valu, out_time

    
    def extractStationDataCSV(self, daterange, variable, stations, file_out):
        """
        Extracts time series from gridded data based on a list of dictionaries
        describing stations. Time serie(s) are written into csv files,
        variables are rounded to a precision of 3 decimal places.
        
        Args:
            daterange: Date range to interpolate. Format of time range 
                desired (yy, mon, day, hh, min)
            variable: Climate variable to interpolate
            stations: Sites to interpolate. Format of stations desired
                ('name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber)
        
        Returns: 
            Returns a csv file contains time series downscaled surface air 
            temperature and land surface influences (surface level) at given 
            sites.

        Example:
            
            daterange = {'beg' : datetime(2000, 01, 11, 00, 00),
                         'end' : datetime(2000, 01, 11, 06, 00)}
                         
            vairable = 'Temperature'
            
            stations = [{'name':'TAE', 'lat':47.47986561, 'lon':8.904870734, 'ele':85.80},
                        {'name':'AAR', 'lat':47.38799123, 'lon':8.043881188, 'ele':454.87}]
            
            file_out = ['/Users/bincao/OneDrive/pl_obs.csv',
                        '/Users/bincao/OneDrive/dT.csv']
            
            downscaling.extractStationDataCSV(daterange,variable,stations,file_out) 
        """
        
        #interpolate and extract values
        out_xyz_dem, lats, lons, shape, names = self.demGrid(stations)
        out_xyz_sur = self.surGrid(lats, lons, stations)
        values, time = self.stationTimeSeries(variable, daterange,
                                              out_xyz_sur, out_xyz_dem)

        #write CSV
        names.insert(0, 'Time_UTC')
        for i in range(len(file_out)):
              with open(file_out[i], 'wb') as output_file:
                 writer = csv.writer(output_file)
                 writer.writerow(names)
                 valu = values[i].tolist()
                 for n in range(len(time)):
                     row = ['%.3f' % elem for elem in valu[n]]
                     row.insert(0,time[n])
                     writer.writerow(row)
    
                     
    
    def extractSpatialDataNCF(self, daterange, variable, file_out):
        """
        Extracts mean of given date range from gridded data based on a 
        fine-scale DEM. Mean values are written into a netcdf file,
        variables are rounded to a precision of 3 decimal places.
        
        Args:
            daterange: Date range to interpolate. Format of time range 
                       desired (yy, mon, day, hh, min)
            variable: Climate variable to interpolate
            filout: Name of output netcdf file
        
        Returns:
            A netcdf file contains coordinate information, downscaled surface
            air temperature and land surface influence (surface level).
            

        Example:
            
            daterange = {'beg' : datetime(2000, 01, 11, 00, 00),
                         'end' : datetime(2000, 01, 11, 06, 00)}
                         
            vairable = 'Temperature'
            
            stations = [{'name':'TAE', 'lat':47.47986561, 'lon':8.904870734, 'ele':85.80},
                        {'name':'AAR', 'lat':47.38799123, 'lon':8.043881188, 'ele':454.87}]
            
            file_out  = '/Users/bincao/Desktop/subgrid.nc'
            
            downscaling.extractSpatialDataNCF(daterange, variable, file_out) 
        """
        out_xyz_dem, lats, lons, shape = self.demGrid()
        out_xyz_sur = self.surGrid(lats, lons, None)
        pl,dt = self.spatialMean(variable, daterange, out_xyz_sur, 
                                 out_xyz_dem, shape)
        
        #create nc file
        nc_root = nc.Dataset(file_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', shape[0])
        nc_root.createDimension('lon', shape[1])
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes = nc_root.createVariable('lat', 'f4', ('lat'))
        Ta = nc_root.createVariable('surface air temperature', 
                                     'f4', ('lat', 'lon'), zlib = True)
        dT = nc_root.createVariable('land surface influence at surface level', 
                                     'f4', ('lat', 'lon'), zlib = True)
        
        #assign variables
        lons = self.dem.variables['lon'][:]
        lats = self.dem.variables['lat'][:]
        longitudes[:] = lons
        latitudes[:] = lats
        Ta[:] = pl
        dT[:] = dt
        
        #attribute
        nc_root.description = "Downscaled surface air temperature and land "\
                              "surface influence (surface level) with a"\
                              "spatial resolution of input dem"
        longitudes.units = 'degree_east (demical)'
        latitudes.units = 'degree_north (demical)'
        Ta.units = 'celsius'
        dT.units = 'celsius'
        
        nc_root.close()
