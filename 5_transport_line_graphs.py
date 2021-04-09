# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:38:48 2020

@author: cristinarusso
"""

#%% Importing modules
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

#%% Importing Data
path = '/home/ocims_platform/transport_auto/transport_data/'
fname = 'total_transports_all_sections_2.nc'
path_out = '/home/ocims_platform/transport_auto/transport_figures/'

ds=xr.open_dataset(path+fname)

#%% Plotting Line Graphs - 20

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across satellite track 20', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_20[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_20,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_20[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'20_total_transport.png',bbox_inches='tight')

#%% Plotting Line Graphs - 96

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across satellite track 96', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_96[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_96,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_96[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'96_total_transport.png',bbox_inches='tight')

#%% Plotting Line Graphs - 172

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across satellite track 172', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_172[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_172,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_172[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'172_total_transport.png',bbox_inches='tight')

#%% Plotting Line Graphs - 198

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across satellite track 198', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_198[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_198,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_198[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'198_total_transport.png',bbox_inches='tight')

#%% Plotting Line Graphs - 248

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across satellite track 248', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_248[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_248,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_248[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'248_total_transport.png',bbox_inches='tight')

#%% Plotting Line Graphs - Benguela

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Total Transport across 34S', fontweight='bold')
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_beng[-1],color= 'coral',marker='o',markersize=15)
h2=plt.plot(ds.time - np.timedelta64(12,'h'),ds.total_transport_beng,'k-')
plt.ylabel('Total Transport (Sv)',fontsize=14,fontweight='bold')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'));
plt.rcParams.update({'font.size': 14})
plt.plot(ds.time[-1] - np.timedelta64(12,'h'),ds.total_transport_beng[-1],'ko',markersize=6)
plt.grid(True)

fig.savefig(path_out+'beng_total_transport.png',bbox_inches='tight')


