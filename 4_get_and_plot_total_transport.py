# Importing modules
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import shutil

#%% Track 20
# Paths
path = '/home/ocims_platform/transport_auto/transport_data/'
fname = 'mercator_section_20_transport.nc'
path_out = '/home/ocims_platform/transport_auto/transport_figures/'
# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_20 = np.sum(hycom_transport[:,0:np.where(section_distance>=300)[0][0]])


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport across satellite track 20', fontweight='bold')
h2=plt.pcolor(section_distance[0:np.where(section_distance>=300)[0][0]],
           -depth[:],
           hycom_transport[:,0:np.where(section_distance>=300)[0][0]],
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'20_transport.png',bbox_inches='tight')

#%% Track 96
# Paths
fname = 'mercator_section_96_transport.nc'

# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_96 = np.sum(hycom_transport[:,0:np.where(section_distance>=300)[0][0]])


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport across satellite track 96', fontweight='bold')
h2=plt.pcolor(section_distance[0:np.where(section_distance>=300)[0][0]],
           -depth[:],
           hycom_transport[:,0:np.where(section_distance>=300)[0][0]],
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'96_transport.png',bbox_inches='tight')

#%% Track 172
# Paths
fname = 'mercator_section_172_transport.nc'

# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_172 = np.sum(hycom_transport[:,0:np.where(section_distance>=300)[0][0]])


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport across satellite track 172', fontweight='bold')
h2=plt.pcolor(section_distance[0:np.where(section_distance>=300)[0][0]],
           -depth[:],
           hycom_transport[:,0:np.where(section_distance>=300)[0][0]],
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'172_transport.png',bbox_inches='tight')

#%% Track 198
# Paths
fname = 'mercator_section_198_transport.nc'
# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_198 = np.sum(hycom_transport[:,0:np.where(section_distance>=300)[0][0]])


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)

fig = plt.figure(figsize =(15, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport across satellite track 198', fontweight='bold')
h2=plt.pcolor(section_distance[0:np.where(section_distance>=300)[0][0]],
           -depth[:],
           hycom_transport[:,0:np.where(section_distance>=300)[0][0]],
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'198_transport.png',bbox_inches='tight')

#%% Track 248

# Paths
fname = 'mercator_section_248_transport.nc'

# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_248 = np.sum(hycom_transport[:,0:np.where(section_distance>=300)[0][0]])


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)
fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport across satellite track 248', fontweight='bold')
h2=plt.pcolor(section_distance[0:np.where(section_distance>=300)[0][0]],
           -depth[:],
           hycom_transport[:,0:np.where(section_distance>=300)[0][0]],
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'248_transport.png',bbox_inches='tight')

#%% Benguela S

# Paths
fname = 'mercator_section_benguela_s_transport.nc'

# Importing Transport Dataset
ds = xr.open_dataset(path+fname)

# Importing Variables
hycom_transport = ds.model_transport[:,:,0]
depth = ds.depth
section_distance = ds.section_distance.values

# Computing total mean transport 
total_trans_beng = np.sum(hycom_transport)


# Plotting Transport Figure 

#midnorm = mpl.colors.DivergingNorm(vcenter=0)

fig = plt.figure(figsize =(10, 5),facecolor='none', edgecolor='none')
ax  = fig.add_axes([0.1, 0.135, 0.8, 0.8])

ax.set_title('Transport Southern Benguela (34S)', fontweight='bold')
h2=plt.pcolor(section_distance,
           -depth[:],
           hycom_transport,
             cmap='seismic_r',
             vmin=-0.05,vmax=0.05)


ax.set_xlabel('Section Distance (km)', fontweight = 'bold')
cbar_ax = fig.add_axes([0.915, 0.3, 0.02,0.4])
fig.colorbar(h2,cax = cbar_ax)
fig.text(0.915,0.715,'Transport \n(Sv)',fontweight = 'bold')
fig.text(0.02,0.575,'Depth (m)',fontweight = 'bold',rotation = 90)
plt.rcParams.update({'font.size': 16})
fig.savefig(path_out+'beng_transport.png',bbox_inches='tight')

#%% Total transport netCDFs

# Import original netcdf file
fname_old = 'total_transports_all_sections_2.nc'
df=xr.open_dataset(path+fname_old)

# Sorting out the time dimension
# adding new date to variable with old dates
time_new =ds.time
time_old = df.time
time_together=np.concatenate((time_old,time_new),axis=0)
time_together = np.squeeze(time_together)

# Importing old transports
tot_20 = df.total_transport_20
tot_96 = df.total_transport_96
tot_172 = df.total_transport_172
tot_198= df.total_transport_198
tot_248= df.total_transport_248
tot_beng= df.total_transport_beng

# Expanding dimension of old variables for concatenation
#tot_20 =np.expand_dims(tot_20,axis=0)
#tot_96 = np.expand_dims(tot_96,axis=0)
#tot_172 = np.expand_dims(tot_172,axis=0)
#tot_198 = np.expand_dims(tot_198,axis=0)
#tot_248 = np.expand_dims(tot_248,axis=0)
#tot_beng = np.expand_dims(tot_beng,axis=0)

# Expanding dimension of new variables for concatenation
total_trans_20 =np.expand_dims(total_trans_20,axis=0)
total_trans_96 = np.expand_dims(total_trans_96,axis=0)
total_trans_172 = np.expand_dims(total_trans_172,axis=0)
total_trans_198 = np.expand_dims(total_trans_198,axis=0)
total_trans_248 = np.expand_dims(total_trans_248,axis=0)
total_trans_beng = np.expand_dims(total_trans_beng,axis=0)

# Concatenating old variables with new variables
tot_20 =np.concatenate((tot_20, total_trans_20),axis=0)
tot_96 =np.concatenate((tot_96, total_trans_96),axis=0)
tot_172 =np.concatenate((tot_172, total_trans_172),axis=0)
tot_198 =np.concatenate((tot_198, total_trans_198),axis=0)
tot_248 =np.concatenate((tot_248, total_trans_248),axis=0)
tot_beng =np.concatenate((tot_beng, total_trans_beng),axis=0)

# Writing netCDF with transport timeseries
tot_trans_dict = {'total_transport_20':(['time'],tot_20),'total_transport_96':(['time'],tot_96),'total_transport_172':(['time'],tot_172),'total_transp>
sec_data = xr.Dataset(tot_trans_dict)
fname = 'total_transports_all_sections_1.nc'
sec_data.to_netcdf(path + fname,  unlimited_dims={'time':True}, format='NETCDF4_CLASSIC',mode='w')

# Making copy of netcdf and saving into same directory
shutil.copy(path+fname, path+'total_transports_all_sections_2.nc')

# Making copy of netcdf and saving into same directory
shutil.copy(path+fname, path+'total_transports_all_sections_3.nc')


