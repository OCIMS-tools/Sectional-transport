# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 09:15:41 2020

@author: cristinarusso
"""

import xarray as xr
import numpy as np
import math
import sectional_transport

#%%
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

# track 20 transport
fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_20_data.nc'
path = '/home/ocims_platform/transport_auto/transport_data/'


# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]


# Calculate along and cross track velocities from model data
u_data = uvdata.u_vel #u velocity component
v_data = uvdata.v_vel #v velocity component

dy = uvdata.Lat[-1].values - uvdata.Lat[0]
dx = uvdata.Lon[-1].values - uvdata.Lon[0]

along_track, cross_track = uv_rotate(u_data, v_data,
                                     -math.degrees(math.atan2(dy, dx)))  

# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(cross_track[0,:,1])-1,len(cross_track[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = cross_track.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':>
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_20_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')
#%%
#%%


fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_96_data.nc'
# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]


# Calculate along and cross track velocities from model data
u_data = uvdata.u_vel #u velocity component
v_data = uvdata.v_vel #v velocity component

dy = uvdata.Lat[-1].values - uvdata.Lat[0]
dx = uvdata.Lon[-1].values - uvdata.Lon[0]

along_track, cross_track = uv_rotate(u_data, v_data,
                                     -math.degrees(math.atan2(dy, dx)))  

# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(cross_track[0,:,1])-1,len(cross_track[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = cross_track.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':uvdata.time, 'depth':uvdata.depth[1:], 
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_96_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')


#%% Track 172 transport
fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_172_data.nc'


# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]


# Calculate along and cross track velocities from model data
u_data = uvdata.u_vel #u velocity component
v_data = uvdata.v_vel #v velocity component

dy = uvdata.Lat[-1].values - uvdata.Lat[0]
dx = uvdata.Lon[-1].values - uvdata.Lon[0]

along_track, cross_track = uv_rotate(u_data, v_data,
                                     -math.degrees(math.atan2(dy, dx)))  

# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(cross_track[0,:,1])-1,len(cross_track[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = cross_track.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':uvdata.time, 'depth':uvdata.depth[1:], 
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_172_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% Track 172 transport
fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_198_data.nc'


# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]

# Calculate along and cross track velocities from model data
u_data = uvdata.u_vel #u velocity component
v_data = uvdata.v_vel #v velocity component

dy = uvdata.Lat[-1].values - uvdata.Lat[0]
dx = uvdata.Lon[-1].values - uvdata.Lon[0]

along_track, cross_track = uv_rotate(u_data, v_data,
                                     -math.degrees(math.atan2(dy, dx)))  
                                     
# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(cross_track[0,:,1])-1,len(cross_track[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = cross_track.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':uvdata.time, 'depth':uvdata.depth[1:], 
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_198_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% Track 172 transport
fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_248_data.nc'

# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]


# Calculate along and cross track velocities from model data
u_data = uvdata.u_vel #u velocity component
v_data = uvdata.v_vel #v velocity component

dy = uvdata.Lat[-1].values - uvdata.Lat[0]
dx = uvdata.Lon[-1].values - uvdata.Lon[0]

along_track, cross_track = uv_rotate(u_data, v_data,
                                     -math.degrees(math.atan2(dy, dx)))  
                                     
# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(cross_track[0,:,1])-1,len(cross_track[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = cross_track.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':uvdata.time, 'depth':uvdata.depth[1:], 
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_248_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')

#%% Track benguela south transport
fname = '/home/ocims_platform/transport_auto/transport_data/mercator_section_benguela_s_data.nc'


# Importing model data
uvdata = xr.open_dataset(fname)

dist = uvdata.section_distance.values * 1000
m_depth = uvdata.depth.values
Model_depth = [round(float(i), 2) for i in m_depth]


# Calculate along and cross track velocities from model data
v_data = uvdata.v_vel #v velocity component


# Calculate Model and ACT sectional transport
hycom_sectional_data = np.zeros([len(v_data[0,:,1])-1,len(v_data[0,1,:])-1,1]) #created empty matrix which will be filled

cross_track = v_data.transpose('depth','section_distance','time')

Lon = uvdata.Lon
Lat = uvdata.Lat

for i in range(len(cross_track[1,1,:])):

        trans1, trans_tot = sectional_transport.transport(cross_track[:,:,i],uvdata.depth,uvdata.section_distance*1000)
        hycom_sectional_data[:,:,i] = trans1
        print(i)

# Creating netCDF for transport data
sec_data_dict = {'model_transport':(['depth','section_distance','time'],hycom_sectional_data),'section_distance':uvdata.section_distance[:-1], 'time':uvdata.time, 'depth':uvdata.depth[1:], 
                 'Lon':Lon, 'Lat':Lat}

sec_data = xr.Dataset(sec_data_dict)

fname = 'mercator_section_benguela_s_transport.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC')






