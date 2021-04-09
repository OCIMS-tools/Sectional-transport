  GNU nano 4.9.2                                                 2_get_mercator_sections.py                                                            
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:35:56 2020

@author: cristinarusso
"""

import xarray as xr
import numpy as np
import pandas as pd


#%%
def geo_idx(dd, dd_array):
    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx


from math import radians, cos, sin, asin, sqrt

AVG_EARTH_RADIUS = 6371  # in km


def distance(lon1, lat1, lon2, lat2):
    """
    Calculates distance between two GPS positions (in decimal degrees)
    lon1,lat1 = start
    lon2,lat2 = end
    Returns distance in Km

    """

    from numpy import sin, cos, sqrt, arctan2, radians, add, subtract, multiply
    R = 6373.0  # Approximate earth radius in Km

    dlon = radians(lon2) - radians(lon1)
    dlat = radians(lat2) - radians(lat1)

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * arctan2(sqrt(a), sqrt(1 - a))

    d = R * c
    return d


def pad(data):
    bad_indexes = np.isnan(data)
    good_indexes = np.logical_not(bad_indexes)
    good_data = data[good_indexes]
    interpolated = np.interp(bad_indexes.nonzero()[0],
                             good_indexes.nonzero()[0], good_data)
    data[bad_indexes] = interpolated
    return data


def uv_rotate(u, v, angle):
    """
    USAGE:
        u_rotate = u*np.cos(angle) - v*np.sin(angle)
        v_rotate = u*np.sin(angle) + v*np.cos(angle)

        Rotate UV components by angle
    Input: 
        U, V, angle in degrees
    Outputs: 
        u_rotate, v_rotate
    """
    angle_rad = angle * np.pi / 180
    u_rot = u * np.cos(angle_rad) - v * np.sin(angle_rad)
    v_rot = u * np.sin(angle_rad) + v * np.cos(angle_rad)
    return u_rot, v_rot

#%%
def glorys_section(data,
                   section_coordinates,
                   time=True,
                   coastal_coordinate=None):

    section_lon, section_lat = section_coordinates

    if coastal_coordinate is None:

        section_dist = []
        for i in range(len(section_lon)):
            section_dist.append(
                distance(section_lon[0], section_lat[0], section_lon[i],
                         section_lat[i]))

        n = len(section_dist)

        if time == True:

            data_1 = np.zeros((len(data[:, 0, 0, 0]), len(data[0, :, 0, 0]),
                               n))  # Time, Depth, Stations

            for tt in range(len(data[:, 0, 0, 0])):
                for i in range(1, n):
                    data_1[tt, :, i] = data[tt, :,
                                            geo_idx(section_lat[i], lat),
                                            geo_idx(section_lon[i], lon)]
        else:
            data_1 = np.zeros((data[0, :, 0, 0], n))
            for i in range(1, n):
                data_1[:, i] = data[30, :,
                                    geo_idx(section_lon[i], lon),
                                    geo_idx(section_lat[i], lat)]

    else:
        coast_lon, coast_lat = coastal_coordinate

        section_dist = []
        for i in range(len(section_lon)):
            section_dist.append(
                distance(coast_lon, coast_lat, section_lon[i], section_lat[i]))

        n = len(section_dist)

        if time == True:

            data_1 = np.zeros((len(data[:, 0, 0, 0]), len(data[0, :, 0, 0]),
                               n))  # Time, Depth, Stations

            for tt in range(len(data[:, 0, 0, 0])):
                for i in range(1, n):
                    data_1[tt, :, i] = data[tt, :,
                                            geo_idx(section_lat[i], lat),
                                            geo_idx(section_lon[i], lon)]
        else:
            data_1 = np.zeros((data[0, :, 0, 0], n))
            for i in range(1, n):
                data_1[:, i] = data[30, :,
                                    geo_idx(section_lon[i], lon),
                                    geo_idx(section_lat[i], lat)]

    return data_1, section_dist
#%%
uvdata = xr.open_dataset('/home/ocims_platform/transport_auto/altimetry_data/mercator_ocean_daily.nc')
uvdata

#%%
udata = uvdata.variables['uo'][:]
vdata = uvdata.variables['vo'][:]

lat = uvdata.variables['latitude']
lon = uvdata.variables['longitude']
depth = uvdata.variables['depth']

np.shape(udata)

#%% ASCA - track 96
my_lon = np.linspace(27.4799, 29.088648,601)
my_lat = np.linspace(-33.3438, -36.03, 601)

transect = (my_lon,my_lat)
u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

sec_dist = []
for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
time = uvdata.time.values

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 'Lat':my_lat,'Lon':my_lon}

sec_data = xr.Dataset(sec_data_dict)


fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_96_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% track 172

my_lon = np.linspace(29.46, 31.31438,601)
my_lat = np.linspace(-31.81, -35.0381, 601)

transect = (my_lon,my_lat)

u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

np.shape(u_section_data)

#%%
sec_dist = []

for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
    
np.shape(sec_dist)
time = uvdata.time.values

#%%

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 
                 'Lon':my_lon, 'Lat':my_lat}

sec_data = xr.Dataset(sec_data_dict)
sec_data

fname ='/home/ocims_platform/transport_auto/transport_data/mercator_section_172_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% track 198

my_lon = np.linspace(22.25, 25.0547,601)
my_lat = np.linspace(-34.11, -38.5432, 601)

transect = (my_lon,my_lat) #(asca.Longitude.values, asca.Latitude.values)

#%%

u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

np.shape(u_section_data)

#%%
sec_dist = []

for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
    
np.shape(sec_dist)
#data = xr.open_mfdataset('/media/jethan/Minotaur/GLORYS_DATA/GLORYS_ACT/UV_data/Regridded/*.nc')
time = uvdata.time.values

#%%

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 
                 'Lon':my_lon, 'Lat':my_lat}

sec_data = xr.Dataset(sec_data_dict)
sec_data

fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_198_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% track 248

my_lon = np.linspace(31.18, 32.9768,601)
my_lat = np.linspace(-29.73 , -33.034078, 601)


transect = (my_lon,my_lat) #(asca.Longitude.values, asca.Latitude.values)

#%%

u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

np.shape(u_section_data)

#%%
sec_dist = []

for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
    
np.shape(sec_dist)
#%%
#data = xr.open_mfdataset('/media/jethan/Minotaur/GLORYS_DATA/GLORYS_ACT/UV_data/Regridded/*.nc')
time = uvdata.time.values

#%%

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 
                 'Lon':my_lon, 'Lat':my_lat}

sec_data = xr.Dataset(sec_data_dict)
sec_data

fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_248_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% track 20

my_lon = np.linspace(25.04, 27.204846, 601)
my_lat = np.linspace(-34.02 , -37.521523, 601)

transect = (my_lon,my_lat) #(asca.Longitude.values, asca.Latitude.values)

#%%

u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

np.shape(u_section_data)

#%%
sec_dist = []

for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
    
np.shape(sec_dist)
#%%
#data = xr.open_mfdataset('/media/jethan/Minotaur/GLORYS_DATA/GLORYS_ACT/UV_data/Regridded/*.nc')
time = uvdata.time.values

#%%

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 
                 'Lon':my_lon, 'Lat':my_lat}

sec_data = xr.Dataset(sec_data_dict)
sec_data

fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_20_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% track benguela_s

my_lon = np.linspace(15.15, 18.33, 601)
my_lat = np.linspace(-34 , -34, 601)

transect = (my_lon,my_lat)
u_section_data, section_dist = glorys_section(udata, transect)
v_section_data, section_dist = glorys_section(vdata, transect)

np.shape(u_section_data)

#%%
sec_dist = []

for i in range(len(section_dist)):
    sec_dist.append(float(section_dist[i]))
    
np.shape(sec_dist)

#%%
#data = xr.open_mfdataset('/media/jethan/Minotaur/GLORYS_DATA/GLORYS_ACT/UV_data/Regridded/*.nc')
time = uvdata.time.values

#%%

sec_data_dict = {'u_vel':(['time','depth','section_distance'],u_section_data), 
            'v_vel':(['time','depth','section_distance'],v_section_data), 
            'section_distance':sec_dist, 'time':time, 'depth':depth, 
                 'Lon':my_lon, 'Lat':my_lat}

sec_data = xr.Dataset(sec_data_dict)
sec_data

fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_benguela_s_data.nc'
sec_data.to_netcdf(fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')


