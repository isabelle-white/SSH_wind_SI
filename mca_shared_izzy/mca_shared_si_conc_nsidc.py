#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
import xeofs as xe
from geometry_izzyv1 import grad_sphere
from regression_izzyv1 import linregress_3D
import cartopy.crs as ccrs
crs = ccrs.PlateCarree()
import warnings
import os
warnings.filterwarnings('ignore')



# In[3]:


# constant
rho_a = 1.225 # air density kg m-3

Cd = 1.5e-3 # drag coeff


# In[4]:


gridtype = 'remapcon'


# In[5]:


# dot dataset
path = '/Users/iw2g24/PycharmProjects/SSH_project/'
print(path)
ds = xr.open_dataset(path + 'Data/dot_all_30bmedian_goco05c_sig3_1.nc')
print(ds)
# ERA5 remap dataset
ds_era5 = xr.open_dataset(path + 'Data/ERA5_u10v10_Antartica_2000-2023_' + gridtype + '.nc')


# In[6]:


print(ds_era5)


# In[7]:


#sea ice drift data from nsidc
directory = "/Users/iw2g24/PycharmProjects/SSH_project/mca_shared_izzy/si_conc_nsidc_G02202_V4/"
ds_example = xr.open_dataset(directory + 'seaice_conc_monthly_sh_197902_n07_v04r00.nc')
print(ds_example)
print(ds_example['xgrid'])

ds_example2 = xr.open_dataset(directory + 'seaice_conc_monthly_sh_198001_n07_v04r00.nc')
print(ds_example2)
print(ds_example2['xgrid'])

variable_name = "cdr_seaice_conc_monthly"

# Extract the data array for plotting
data_to_plot = ds_example2[variable_name].isel(tdim=0)
print("plotted data size",data_to_plot.shape)

# Plot using xarray's built-in method
data_to_plot.plot(cmap="coolwarm")  # Replace 'viridis' with your preferred color map
#plt.title(f"{variable_name} at time {ds_example['time']}")
plt.show()

print(ds_example2["xgrid"].shape)
print(ds_example2["ygrid"].shape)


# In[8]:


print(ds_example2["cdr_seaice_conc_monthly"].attrs)


# In[9]:


print(ds_example2.variables)


# In[10]:


ds_check = xr.open_dataset('nsidc_conc2.nc')
print(ds_check)


# In[11]:


ds_nsidc_conc = xr.open_dataset('nsidc_conc2.nc')
print(ds_nsidc_conc)


# In[12]:


data_to_plot = ds_nsidc_conc["cdr_seaice_conc_monthly"]
print("xgrid shape:", ds_nsidc_conc["xgrid"].shape)
print("ygrid shape:", ds_nsidc_conc["ygrid"].shape)
print("Data to plot shape:", data_to_plot.isel(time_new = 0, tdim=0).shape)


# In[13]:


# Plot sea ice conc using x and y coordinates
data_to_plot = ds_nsidc_conc["cdr_seaice_conc_monthly"]

plt.figure(figsize=(7,7))
plt.pcolormesh(
    ds_nsidc_conc["xgrid"],
    ds_nsidc_conc["ygrid"],
    data_to_plot.isel(time_new=0, tdim=0),
    cmap="coolwarm",
    shading="auto"
)
plt.colorbar(label="Sea Ice Conc")
plt.title("Sea Ice Drift with xgrid/ygrid")
plt.xlabel("xgrid")
plt.ylabel("ygrid")
plt.show()


# In[14]:


print("Data shape:", ds_nsidc_conc["cdr_seaice_conc_monthly"].isel(tdim=0).shape)
print("Longitude shape:", ds_nsidc_conc["lon"].shape)
print("Latitude shape:", ds_nsidc_conc["lat"].shape)


# In[15]:


print(ds_nsidc_conc["lon"])


# In[ ]:





# In[16]:


print(ds_nsidc_conc.values)
print(ds_nsidc_conc['lon'].shape)
print(ds_nsidc_conc['lat'].shape)


# In[17]:


print(ds_nsidc_conc['time_new'].dtype)


# In[18]:


time_nsidc = ds_nsidc_conc['time_new'].values
print("nsidc conc time", time_nsidc[284],time_nsidc[479])


# In[19]:


#reorder ds_nsidc_conc  to have time_new, lat, lon

# Open the dataset
ds_nsidc_conc = xr.open_dataset("nsidc_conc2.nc")

# Step 1: Extract lat, lon, and sea ice concentration
lat_nsidc = ds_nsidc_conc['lat'].values  # Shape (y, x)
lon_nsidc = ds_nsidc_conc['lon'].values  # Shape (y, x)
time_new = ds_nsidc_conc['time_new'].values  # Shape (time_new)

# Step 2: Flatten the lat and lon arrays into 1D arrays
lat_flat = lat_nsidc.flatten()
lon_flat = lon_nsidc.flatten()

# Step 3: Reshape the variable (e.g., 'sea_ice_conc') to match the flattened grid
var = ds_check['cdr_seaice_conc_monthly'].values  # Shape (time_new, tdim, y, x)

# Remove 'tdim' dimension by squeezing it
var_squeezed = var[:, 0, :, :]  # Shape becomes (time_new, y, x)

# Flatten 'y' and 'x' into a single dimension
var_flat = var_squeezed.reshape((var_squeezed.shape[0], -1))  # Shape (time_new, y*x)

# Step 4: Create a new dataset with dimensions (time_new, lat, lon)
ds_reordered = xr.Dataset(
    data_vars={
        'si_conc': (['time_new', 'points'], var_flat)  # Data flattened
    },
    coords={
        'time_new': time_new,      # Time dimension
        'lat': ('points', lat_flat),  # 1D lat coordinate
        'lon': ('points', lon_flat)   # 1D lon coordinate
    }
)

# Step 5: Verify the new dataset
print(ds_reordered)

# Step 6: Save the new dataset
output_file = "reordered_sea_ice_concentration.nc"
ds_reordered.to_netcdf(output_file)
print(f"Reordered dataset saved to: {output_file}")


# In[20]:


print(f"Number of points: {lat_flat.shape[0]} (should be equal to y * x)")
print(f"Sea Ice Concentration reshaped to: {var_flat.shape}")
print(332*316)


# In[21]:


# change coordinates from (lon,lat,time) to (time,lat,lon) as ERA5
lon_dot = ds['longitude'].values; lat_dot = ds['latitude'].values
dot_all = ds['dot'].values
dot = np.zeros((dot_all.shape[2], dot_all.shape[1], dot_all.shape[0]))
for t in range(0,dot_all.shape[2]):
    for i in range(0,dot_all.shape[1]):
        for j in range(0,dot_all.shape[0]):
            dot[t,i,j] = dot_all[j,i,t].copy()
seamask = dot[0].copy()/dot[0]
seamask[seamask == 0] = np.nan
del dot_all


# In[22]:


# make lon lat 2d
llon = np.zeros((dot.shape[1],dot.shape[2]))
for i in range(0,dot.shape[1]):
    llon[i,:] = lon_dot
llat = np.zeros((dot.shape[1],dot.shape[2]))
for i in range(0,dot.shape[2]):
    llat[:,i] = lat_dot


# In[23]:


print(dot.shape)
print(llon.shape)
print(llat.shape)
print(dot[0].shape)


# In[89]:


#plot dot with different colour schemes
# fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(20, 10))
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 3))

pc1 = axes.pcolormesh(llon,llat,dot[0], cmap='seismic')
axes.tick_params(axis='both', which='major', length=10, width=2, labelsize=14)
axes.set_title('SSA [m]', fontsize=20)
axes.set_xlabel('Longitude', fontsize=15)
axes.set_ylabel('Latitude', fontsize=15)
cbar1 = fig.colorbar(pc1, ax=axes)
cbar1.set_label('Dot [m]', fontsize=15)
cbar1.ax.tick_params(labelsize=14)

# pc1 = axes[0].pcolormesh(llon,llat,dot[0], cmap='seismic')
# axes[0].tick_params(axis='both', which='major', length=10, width=2, labelsize=14)
# axes[0].set_title('SSA [m]', fontsize=20)
# axes[0].set_xlabel('Longitude', fontsize=15)
# axes[0].set_ylabel('Latitude', fontsize=15)
# cbar1 = fig.colorbar(pc1, ax=axes[0])
# cbar1.set_label('Dot [m]', fontsize=15)
# cbar1.ax.tick_params(labelsize=14)

# pc2 = axes[1].pcolormesh(llon,llat,dot[0], cmap='coolwarm')
# axes[1].tick_params(axis='both', which='major', length=10, width=2, labelsize=14)
# axes[1].set_title('SSA [m]', fontsize=20)
# axes[1].set_xlabel('Longitude', fontsize=15)
# axes[1].set_ylabel('Latitude', fontsize=15)
# cbar2 = fig.colorbar(pc2, ax=axes[1])
# cbar2.set_label('Dot [m]', fontsize=15)
# cbar2.ax.tick_params(labelsize=14)
#
# pc3 = axes[2].pcolormesh(llon,llat,dot[0], cmap='viridis')
# axes[2].tick_params(axis='both', which='major', length=10, width=2, labelsize=14)
# axes[2].set_title('SSA [m]', fontsize=20)
# axes[2].set_xlabel('Longitude', fontsize=15)
# axes[2].set_ylabel('Latitude', fontsize=15)
# cbar3 = fig.colorbar(pc3, ax=axes[2])
# cbar3.set_label('Dot [m]', fontsize=15)
# cbar3.ax.tick_params(labelsize=14)

plt.tight_layout()
plt.show()


# In[25]:


# time
time_era5 = ds_era5['time'].values
time_dot = ds['time'].values
time_nsidc = ds_nsidc_conc['time_new'].values


# In[26]:


time_dot_array = np.array(time_dot, dtype='datetime64[ns]')
time_nsidc_array = np.array(time_nsidc, dtype='datetime64[ns]')

indices_dot_nsidc = np.where(np.isin(time_nsidc_array, time_dot_array))
print(indices_dot_nsidc)


# In[27]:


#print(time_era5)
print("dot time", time_dot[0],time_dot[-1])
print("era5 time", time_era5[30],time_era5[225])
print("nsidc conc time", time_nsidc[284],time_nsidc[479])


# In[28]:


tstart_era5 = 30
tend_era5 = 226 # +1 respect cell above

tstart_nsidc = 284
tend_nsidc = 480 #+1 respect cell above


# In[90]:


# # select sea ice concentration
# seaice_concentration = ds_nsidc_conc['cdr_seaice_conc_monthly'][tstart_nsidc:tend_nsidc].values
# print(seaice_concentration.shape)
# print(seaice_concentration.max())


# In[30]:


# DE-TREND DOT


# In[31]:


# linear regression in time on each grid point
n,slope,intercept,p_val,r_square,rmse = linregress_3D(dot)
# remove trend for dot
xt = np.zeros(dot.shape)
for t in range(0,dot.shape[0]):
    xt[t] = t
dot_detrended = dot - (slope*xt + intercept)


# In[32]:


plt.figure(figsize=(12,5))
plt.pcolormesh(llon, llat, slope, vmin=-6e-4, vmax=6e-4, cmap='seismic')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Slope of Dot')
plt.colorbar()


# In[33]:


dot_final = dot_detrended.copy()


# In[34]:


seaice_concentration = ds_nsidc_conc['cdr_seaice_conc_monthly'][tstart_nsidc:tend_nsidc].values
print(seaice_concentration.shape)

# Assuming 'time_new' is the time coordinate in the dataset
time_values = ds_nsidc_conc['time_new'][tstart_nsidc:tend_nsidc].values
print(time_values.shape)  # To verify the shape


# In[35]:


print(seaice_concentration.shape)
print((seaice_concentration).shape[0])
print((seaice_concentration[t,0]).shape)


# In[36]:


#remap data onto dot grid
# Target grid
lon_target = ds['longitude'].values  # (360,)
lat_target = ds['latitude'].values  # (64,)

# Create an empty list to store regridded data
regridded_data_list = []
time_steps_list = []


# Loop through each time step of the original sea ice data
for t in range(seaice_concentration.shape[0]):  # Loop over the time dimension
    # Extracting the sea ice concentration at time t, shape (332, 316)
    seaice = seaice_concentration[t, 0]

    # Coordinates for the original grid
    lon_seaice = lon_nsidc  # (332, 316)
    lat_seaice = lat_nsidc  # (332, 316)

    # Flatten original grid and sea ice data
    points = np.array([lon_seaice.flatten(), lat_seaice.flatten()]).T  # [lon, lat]
    values = seaice.flatten()

    # print(points.shape)
    # print(values.shape)


    # Generate target grid (meshgrid for lon and lat)
    lon_mesh, lat_mesh = np.meshgrid(lon_target, lat_target)  # 2D target grid

    # Interpolate onto the target grid
    seaice_regridded = griddata(points, values, (lon_mesh, lat_mesh), method='linear')

    # Add regridded data to list (append 2D array for this time step)
    regridded_data_list.append(seaice_regridded)

    # Store the corresponding time value
    time_steps_list.append(time_values[t])  # Add the corresponding time step

# Stack the regridded data along the time dimension
regridded_seaice_conc = np.stack(regridded_data_list, axis=0)  # Shape: (time, lat, lon)

# Now you can use regridded_data, and time_steps_list to refer to the time dimension
print(f"Regridded data shape: {regridded_seaice_conc.shape}")
#print(f"Time steps: {time_steps_list}")


# In[1]:


print(len(points))
print(len(values))


# In[37]:


# Create a new xarray DataArray with regridded data
regridded_seaice_da = xr.DataArray(
    regridded_seaice_conc,
    dims=("time", "lat", "lon"),
    coords={
        "time": time_steps_list,  # Assuming time steps are just integers
        "lat": lat_target,  # New latitude grid (64,)
        "lon": lon_target,  # New longitude grid (360,)
    },
    name="seaice_concentration",  # Variable name
)

# Create a Dataset with the regridded DataArray
regridded_seaice_ds = xr.Dataset(
    {"seaice_concentration": regridded_seaice_da}
)

# Print dataset to verify
print(regridded_seaice_ds)


# In[38]:


print(regridded_seaice_ds)

# Extract the time, latitude, and longitude coordinates
time_values = regridded_seaice_ds['time'].values
lat_values = regridded_seaice_ds['lat'].values
lon_values = regridded_seaice_ds['lon'].values

# Create a meshgrid for lat/lon coordinates for plotting
lon_mesh, lat_mesh = np.meshgrid(lon_values, lat_values)

# Plot the sea ice concentration at time step 0 (you can change the time index)
plt.figure(figsize=(12, 6))
plt.pcolormesh(llon, llat, regridded_seaice_conc[0], cmap='coolwarm', shading='auto')
plt.colorbar(label='Sea Ice Concentration')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(f'Sea Ice Concentration at Time {time_values[0]}')  # Use time label for title
plt.show()


# In[39]:


# Mask values > 1
regridded_seaice_filtered = np.where(regridded_seaice_conc >1,  np.nan, regridded_seaice_conc)

# Verify the masking
print("Original max value in regridded:", np.nanmax(regridded_seaice_conc))
print("Max value after masking in regridded_seaice_masked:", np.nanmax(regridded_seaice_filtered))

# Plot the sea ice concentration at time step 0 (you can change the time index)
plt.figure(figsize=(12, 6))
plt.pcolormesh(llon, llat, regridded_seaice_filtered[0], cmap='coolwarm', shading='auto')
plt.colorbar(label='Sea Ice Concentration')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(f'Sea Ice Concentration at Time {time_values[0]}')  # Use time label for title
plt.show()


# In[ ]:





# In[40]:


icemask = regridded_seaice_filtered[0].copy()/regridded_seaice_filtered[0]
icemask[icemask == 0] = np.nan


# In[41]:


#apply the icemask
# Apply the mask to the entire seaice_detrended dataset
seaice_masked = regridded_seaice_filtered * icemask[np.newaxis, :, :]


# In[ ]:





# In[ ]:





# In[42]:


# #detrend seaice conc - not using the masked seaice
# #KEEP THIS CELL BUT COMMENT OR UNCOMMENT THIS/CELL BELOW TO USE SEAICE_MASKED OR NOT ...
# # linear regression in time on each grid point
# n_si,slope_si,intercept_si,p_val_si,r_square_si,rmse_si = linregress_3D(regridded_seaice_filtered)
# # remove trend for dot
# xt_si = np.zeros(regridded_seaice_filtered.shape)
# for t in range(0,regridded_seaice_filtered.shape[0]):
#     xt_si[t] = t
# seaice_detrended = regridded_seaice_filtered - (slope_si*xt_si + intercept_si)
# slope_si = slope_si.squeeze()
#
# # Verify shapes
# print("Longitude shape:", lon_mesh.shape)
# print("Latitude shape:", lat_mesh.shape)
# print("Data shape (slope_si):", slope_si.shape)


# In[ ]:





# In[43]:


#detrend seaice conc - using seaice_masked
# linear regression in time on each grid point
n_si,slope_si,intercept_si,p_val_si,r_square_si,rmse_si = linregress_3D(seaice_masked)
# remove trend for dot
xt_si = np.zeros(seaice_masked.shape)
for t in range(0,seaice_masked.shape[0]):
    xt_si[t] = t
seaice_detrended = seaice_masked - (slope_si*xt_si + intercept_si)
slope_si = slope_si.squeeze()

# Verify shapes
print("Longitude shape:", lon_mesh.shape)
print("Latitude shape:", lat_mesh.shape)
print("Data shape (slope_si):", slope_si.shape)


# In[ ]:





# In[44]:


print('slope range: ', slope_si.min(), slope_si.max(), slope_si.mean(), slope_si.std())

print('seaice conc range:', regridded_seaice_filtered.min(), regridded_seaice_filtered.max())


# In[45]:


seaice_final = seaice_detrended.copy()


# In[46]:


plt.figure(figsize=(12,5))
plt.pcolormesh(llon, llat, slope_si, vmin=-6e-4, vmax=6e-4, cmap='coolwarm')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Slope of Sea Ice Concentration Slope')
plt.colorbar()


# In[47]:


#DESEASONAL


# In[48]:


#deseasonal for sea ice and dot
months = time_dot.astype('datetime64[M]').astype(int) % 12 + 1


# In[49]:


# find and remove seasonality from dot_detrended and seaice_detrended
dot_seas = np.zeros((12,dot.shape[1],dot.shape[2])); seaice_seas = np.zeros((12,dot.shape[1],dot.shape[2])); seaice_seas_test = np.zeros((12,dot.shape[1],dot.shape[2]))
for m in range(1,13):
    dot_seas[m-1] = np.nanmean(dot_detrended[np.where(months==m)], axis=0)
    seaice_seas[m-1] = np.nanmean(seaice_detrended[np.where(months==m)], axis=0)
    seaice_seas_test[m-1] = np.nanmean(seaice_masked[np.where(months == m)], axis=0)


# remove seasonality from dot and sea ice concetration
dot_final = np.zeros(dot_detrended.shape); seaice_final = np.zeros(seaice_detrended.shape);seaice_final_test = np.zeros(seaice_detrended.shape)
for i,m in enumerate(months):
    dot_final[i] = dot_detrended[i] - dot_seas[m-1]
    seaice_final[i] = seaice_detrended[i] - seaice_seas[m-1]
    seaice_final_test[i] = seaice_masked[i] - seaice_seas_test[m-1]


# In[50]:


seaice_final.min()


# In[51]:


print("dot_detrended shape:", dot_detrended.shape)
print("seaice_detrended shape:", seaice_detrended.shape)
print("months shape:", months.shape)


# In[52]:


fig = plt.figure(figsize=(10,5))
plt.subplot(3,2,1)
plt.plot(np.nanmean(dot, axis=(1,2)))
plt.ylabel('dot')
plt.subplot(3,2,2)
plt.plot(np.nanmean(regridded_seaice_filtered, axis=(1,2)))
plt.ylabel('seaice conc')
plt.subplot(3,2,3)
plt.plot(np.nanmean(dot_detrended, axis=(1,2)))
plt.ylabel('dot - trend')
plt.subplot(3,2,4)
plt.plot(np.nanmean(seaice_detrended, axis=(1,2)))
plt.ylabel('seaice conc - trend')
plt.subplot(3,2,5)
plt.plot(np.nanmean(dot_final, axis=(1,2)))
plt.ylabel('dot - trend - seas')
plt.subplot(3,2,6)
plt.plot(np.nanmean(seaice_final, axis=(1,2)))
plt.ylabel('seaice conc - trend - seas')

plt.tight_layout()


# In[53]:


fig = plt.figure(figsize=(10,5))
plt.subplot(3,2,1)
plt.plot(np.nanmean(dot, axis=(1,2)))
plt.ylabel('dot')
plt.subplot(3,2,2)
plt.plot(np.nanmean(regridded_seaice_filtered, axis=(1,2)))
plt.ylabel('seaice conc')
plt.subplot(3,2,3)
plt.plot(np.nanmean(dot_detrended, axis=(1,2)))
plt.ylabel('dot - trend')
plt.subplot(3,2,4)
plt.plot(np.nanmean(regridded_seaice_filtered, axis=(1,2)))
plt.ylabel('seaice conc (again)')
plt.subplot(3,2,5)
plt.plot(np.nanmean(dot_final, axis=(1,2)))
plt.ylabel('dot - trend - seas')
plt.subplot(3,2,6)
plt.plot(np.nanmean(seaice_final_test, axis=(1,2)))
plt.ylabel('seaice conc - trend - seas')

plt.tight_layout()


# In[54]:


# SSH trends without removing ocean variability due to sea ice concentration
plt.figure(figsize=(10,3))
n,slope,intercept,p_val,r_square,rmse = linregress_3D(dot_final)
plt.pcolormesh(lon_mesh, lat_mesh, slope, cmap='seismic')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('SSH Trend without removal of variability due to SIC')
cbar = plt.colorbar()
cbar.set_label('Trend')


# In[91]:


time_step = 0  # You can change this to any index from 0 to 195 (since the shape is (196, 64, 360))

# Extract and plot the masked data for the selected time step
# seaice_at_t = regridded_seaice_filtered[time_step]
seaice_at_t = seaice_final[time_step]

# Plot the data using pcolormesh
plt.figure(figsize=(10, 6))
plt.pcolormesh( lon_mesh, lat_mesh,seaice_at_t, cmap='coolwarm', shading='auto')
plt.colorbar(label='Sea Ice Concentration')
plt.title(f"Sea Ice Concentration at Time Step {time_step}")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()


# In[56]:


print(np.nanmin(regridded_seaice_filtered))
regridded_seaice_filtered = np.nan_to_num(regridded_seaice_filtered,0)
print(regridded_seaice_filtered.min())
print(regridded_seaice_filtered.shape)
# regridded_seaice_filtered_swap = np.swapaxes(regridded_seaice_filtered, 1, 2)
# print(regridded_seaice_filtered_swap.shape)


# In[57]:


regridded_sic_xa = xr.DataArray(regridded_seaice_filtered,coords =[time_dot, lat_values, lon_values], dims = ["time", "lat", "lon"], name = "regridded_SIC")

regridded_sic_xa.to_netcdf("regridded_SIC_xa.nc")


# In[58]:


dot_final.shape

# dot_final_swap = np.swapaxes(dot_final, 1, 2)
#seaice_final_swap = np.swapaxes(seaice_final, 1, 2)
# dot_final_swap.shape


# In[59]:


# MCA


# In[60]:


dot_xa = xr.DataArray(dot_final, coords=[time_dot, lat_values, lon_values], dims=["time", "lat", "lon"],)
#dot_swap_xa = xr.DataArray(dot_final_swap, coords=[time_dot, lon_values, lat_values], dims=["time", "lon", "lat"],)
seaice_xa = xr.DataArray(seaice_final, coords=[time_dot, lat_values, lon_values], dims=["time", "lat", "lon"],)
#seaice_swap_xa = xr.DataArray(seaice_final_swap, coords=[time_dot, lon_values, lat_values], dims=["time", "lon", "lat"],)




# print(dot_xa.values)
# print('/n---/n')
# print(seaice_xa.values)

min_value = seaice_xa.where(~seaice_xa.isnull()).min().item()
max_value = seaice_xa.where(~seaice_xa.isnull()).max().item()
print(f"Range of values (excluding NaN): Min = {min_value}, Max = {max_value}")


# In[61]:


total_nans = dot_xa.isnull().sum().item()
print(f"Total NaN values in dot_xa: {total_nans}")

total_nans_seaice = seaice_xa.isnull().sum().item()
print(f"Total NaN values in seaice_xa: {total_nans_seaice}")

total_points_dot = dot_xa.size
print(f"Total data points in dot_xa: {total_points_dot}")

total_points_seaice = seaice_xa.size
print(f"Total data points in seaice_xa: {total_points_seaice}")

print(seaice_xa.shape)
print(dot_xa.shape)



# In[62]:


#plot where there are nan values (this is where i created the mask for SIC > 1 earlier)
time_step = 0
# Create a mask where NaN values are True
seaice_nan_mask = np.isnan(regridded_seaice_filtered[time_step])

# Plot the mask using pcolormesh
plt.figure(figsize=(10, 6))
plt.pcolormesh(
    lon_mesh, lat_mesh, seaice_nan_mask,
    cmap='Greys', shading='auto'
)
plt.colorbar(label='NaN Mask (1 = NaN, 0 = Valid)')
plt.title(f"Sea Ice NaN Locations at Time Step {time_step}")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()


# In[63]:


print(dir(xe))


# In[93]:


#set all NaN values to = 0 in dot and seaice. this means there has been no change
dot_xa_masked = dot_xa.fillna(0)
#dot_xa_masked_swap = dot_swap_xa.fillna(0)
seaice_xa_masked = seaice_xa.fillna(0)
# regridded_sic_xa_masked = regridded_sic_xa.fillna(0)

print("Masked dot_xa shape:", dot_xa_masked.shape)
print("Masked seaice_xa shape:", seaice_xa_masked.shape)
# print('masked regridded shape', regridded_sic_xa_masked.shape)
#print("swapped axes dot shape", dot_xa_masked_swap.shape)


# In[94]:


print("Remaining NaNs in dot_xa_masked:", dot_xa_masked.isnull().sum().item())
print("Remaining NaNs in seaice_xa_masked:", seaice_xa_masked.isnull().sum().item())
print(seaice_xa_masked.min())
# print(regridded_sic_xa_masked.min())


# In[66]:


print(seaice_xa_masked.shape)


# In[67]:


# #create nc file from seaice_xa_masked ds
#
# # Create a new xarray DataArray with the seaice_xa_masked data
# seaice_xa_masked_da = xr.DataArray(
#     seaice_xa_masked,  # Use your masked sea ice data here
#     dims=("time", "lat", "lon"),  # Define the dimensions
#     coords={
#         "time": time_steps_list,  # Replace with your actual time steps
#         "lat": lat_target,  # Latitude coordinate (e.g., (64,))
#         "lon": lon_target,  # Longitude coordinate (e.g., (360,))
#     },
#     name="seaice_xa_masked",  # Give the variable a name
# )
#
# # Create a Dataset with the seaice_xa_masked DataArray
# seaice_xa_masked_ds = xr.Dataset(
#     {"seaice_xa_masked": seaice_xa_masked_da}
# )
#
# # Print dataset to verify its structure
# print(seaice_xa_masked_ds)
#
# # Save the dataset to a NetCDF file
# output_file = "seaice_xa_masked_nsdic.nc"
# seaice_xa_masked_ds.to_netcdf(output_file)
#
# print(f"NetCDF file created: {output_file}")




# In[ ]:





# In[68]:





# In[106]:


import inspect
model = xe.cross.MCA(n_modes=22, standardize=True)
#print(dir(model))
print(dot_xa_masked.coords)
print(seaice_xa_masked.coords)

model.fit(dot_xa_masked, seaice_xa_masked, dim='time')
comps1, comps2 = model.components()  # Singular vectors (spatial patterns)
scores1, scores2 = model.scores()  # Expansion coefficients (temporal patterns)


# In[107]:


s1 = scores1.values; s2 = scores2.values
c1 = comps1.values; c2 = comps2.values


# In[108]:


print(path)
print(gridtype)


# In[109]:


# Full path to the 'results' directory
results_dir = os.path.join(path, 'results')

# Check if the directory already exists
if os.path.exists(results_dir):
    print("Directory already exists")
else:
    # Create the directory if it doesn't exist
    os.makedirs(results_dir, exist_ok=True)
    print("Directory created!")


# In[110]:


# write PCs on files
# fname = path + 'results/PCs_dot_' + gridtype + '.txt'
fname = path + 'results/PCs_dot_' + gridtype + '_seas.txt'
with open(fname, 'w') as file:
    file.write('time\tPC1\tPC2\tPC3\n')
    for i in range(len(time_dot)):
        file.write(f'{time_dot[i]}\t{s1[0,i]}\t{s1[1,i]}\t{s1[2,i]}\n')
file.close()

# fname = path + 'results/PCs_ws_curl_' + gridtype + '.txt'
fname = path + 'results/PCs_si_concentration_' + gridtype + '_seas.txt'
with open(fname, 'w') as file:
    file.write('time\tPC1\tPC2\tPC3\n')
    for i in range(len(time_dot)):
        file.write(f'{time_dot[i]}\t{s2[0,i]}\t{s2[1,i]}\t{s2[2,i]}\n')
file.close()


# In[111]:


#how can you change from using cross.MCA to calling the subplots PCA1, PCA2, PCA3 etc ...


# In[112]:


r_pears = []; p_pears = []
r_spear = []; p_spear = []
for m in range(1,5):
    x = scores1.sel(mode=m)
    y = scores2.sel(mode=m)
    r, p = scipy.stats.pearsonr(x, y)
    r_pears.append(np.round(r, 2)); p_pears.append(np.round(p, 2))
    r, p = scipy.stats.spearmanr(x, y)
    r_spear.append(np.round(r, 2)); p_spear.append(np.round(p, 2))


# In[113]:


fig = plt.figure(figsize=(12,10))

for i in range(0,4):
    j = 3*i+1 # a = 2*(i+1) + i - 1
    
    plt.subplot(4,3,j)
    scores1.sel(mode=i+1).plot(label = 'SSH')
    scores2.sel(mode=i+1).plot(label = 'SIC')
    plt.title('r = ' +str(r_spear[i]))
    plt.xlabel('Year'); plt.ylabel('PC'+str(int(i+1)))
    plt.xticks(rotation=30) 
    plt.legend()
    
    plt.subplot(4,3,j+1)
    comps1.sel(mode=i+1).plot()
    plt.title('SSH')
    
    plt.subplot(4,3,j+2)
    comps2.sel(mode=i+1).plot()
    plt.title('Sea Ice Concentration')

fig.tight_layout()

fig.savefig(path + 'results/MCA_output_' + gridtype + '_seas.png')


# In[114]:


# compute squared covariance fraction (equivalent to explained variance in PCA)
sq_cov_frac = model.squared_covariance_fraction()
scf = sq_cov_frac.values
cum_sum = np.cumsum(scf)
print('Cumulative sum of squares:', cum_sum)
print('Squared covariance fraction:',scf)


# In[115]:


print(type(cum_sum))
print(len(cum_sum))
mode_percent = []
for i in range(len(cum_sum)):
    if i == 0:
        each_mode_percent =float(cum_sum[i])
        mode_percent.append(each_mode_percent)
    else:
        each_mode_percent = float(cum_sum[i]-cum_sum[i-1])
        mode_percent.append(each_mode_percent)
print('PCA fraction', mode_percent)


# In[ ]:





# In[ ]:





# In[120]:


import inspect

model = xe.cross.MCA(n_modes=22, standardize=True)
#print(dir(model))
# print(dot_xa_masked.coords)
# print(seaice_xa_masked.coords)

model.fit(dot_xa_masked, regridded_sic_xa, dim='time')
comps1, comps2 = model.components()  # Singular vectors (spatial patterns)
scores1, scores2 = model.scores()  # Expansion coefficients (temporal patterns)
s1 = scores1.values;
s2 = scores2.values
c1 = comps1.values;
c2 = comps2.values

#how can you change from using cross.MCA to calling the subplots PCA1, PCA2, PCA3 etc ...
r_pears = [];
p_pears = []
r_spear = [];
p_spear = []
for m in range(1, 5):
    x = scores1.sel(mode=m)
    y = scores2.sel(mode=m)
    r, p = scipy.stats.pearsonr(x, y)
    r_pears.append(np.round(r, 2));
    p_pears.append(np.round(p, 2))
    r, p = scipy.stats.spearmanr(x, y)
    r_spear.append(np.round(r, 2));
    p_spear.append(np.round(p, 2))
fig = plt.figure(figsize=(12, 10))

for i in range(0, 4):
    j = 3 * i + 1  # a = 2*(i+1) + i - 1

    plt.subplot(4, 3, j)
    scores1.sel(mode=i + 1).plot(label='SSH')
    scores2.sel(mode=i + 1).plot(label='SIC')
    plt.title('r = ' + str(r_spear[i]))
    plt.xlabel('Year');
    plt.ylabel('PC' + str(int(i + 1)))
    plt.xticks(rotation=30)
    plt.legend()

    plt.subplot(4, 3, j + 1)
    comps1.sel(mode=i + 1).plot()
    plt.title('SSH')

    plt.subplot(4, 3, j + 2)
    comps2.sel(mode=i + 1).plot()
    plt.title('Sea Ice Concentration')

fig.tight_layout()

# compute squared covariance fraction (equivalent to explained variance in PCA)
sq_cov_frac = model.squared_covariance_fraction()
scf = sq_cov_frac.values
cum_sum = np.cumsum(scf)
print('Cumulative sum of squares:', cum_sum)
print('Squared covariance fraction:', scf)
print(type(cum_sum))
print(len(cum_sum))
mode_percent = []
for i in range(len(cum_sum)):
    if i == 0:
        each_mode_percent = float(cum_sum[i])
        mode_percent.append(each_mode_percent)
    else:
        each_mode_percent = float(cum_sum[i] - cum_sum[i - 1])
        mode_percent.append(each_mode_percent)
print('PCA fraction', mode_percent)


# In[ ]:




