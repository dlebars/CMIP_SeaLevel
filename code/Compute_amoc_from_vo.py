###############################################################################
# Compute_amoc_from_vo.py:
# Compute time series of AMOC at different latitudes using the meridional 
# velocity fields from CMIP5 amd CMIP6
###############################################################################

import itertools
import os

import numpy as np
import xarray as xr

import mod_loc as loc

verbose = True # Print additional information
LAT_SEL = [26, 35] # Latitude to select
MIP = 'cmip6'

if MIP == 'cmip5':
    VAR = 'vo'
    EXP = ['historical', 'rcp26', 'rcp45','rcp85']
    dir_outputs = '/nobackup/users/bars/CMIP5_regridded/'
    
elif MIP == 'cmip6':
    VAR = 'vo' # msftmz, msftyz
    EXP = ['historical'] #, 'ssp126', 'ssp245', 'ssp585']

dir_inputs = '../inputs/'
dir_outputs = f'/nobackup/users/bars/{MIP.upper()}_regridded/'

def open_files(file_names):
    try:
        ds = xr.open_mfdataset(file_names, combine='by_coords', use_cftime=True)
    except:
        print(f'Open by_coords did not work for:{file_names}')
        print('Using nested option instead')
        try:
            ds = xr.open_mfdataset(file_names, combine='nested', 
                                   concat_dim='time', 
                                   use_cftime=True)
        except:
            print('Nested option did not work either. Decoding time manually.')
            
            try:
                ds = xr.open_mfdataset(file_names, combine='by_coords', decode_times = False)
                ds.time.attrs['units'] = 'days since 0000-01-01'
                ds['time'].data[:] = ds.time.data[:] + 365
                ds = xr.decode_cf(ds)
                
            except:
                print('ERROR: Reading dataset was not successful')

    return ds

def harmonize_lat_lon_names(ds):
    '''Harmonise the name of longitude and latitude variables in a CMIP dataset'''
    
    if ('lat' not in ds.coords):
        if 'latitude' in ds.coords:
            ds = ds.rename({'latitude':'lat'})
        elif 'rlat' in ds.coords:
            ds = ds.rename({'rlat':'lat'})
        elif 'nav_lat' in ds.coords:
            ds = ds.rename({'nav_lat':'lat'})

    if ('lon' not in ds.coords):
        if 'longitude' in ds.coords:
            ds = ds.rename({'longitude':'lon'})
        elif 'rlon' in ds.coords:
            ds = ds.rename({'rlon':'lon'})
        elif 'nav_lon' in ds.coords:
            ds = ds.rename({'nav_lon':'lon'})
            
    return ds

def change_longitude_values(ds): 
        
    print('Changing longitude from -180,0 to 180,360')
    
    new_lon = np.where(ds['lon'].values<0, ds['lon'].values+360, ds['lon'].values)
    ds['lon'].values = new_lon
    
    return ds

for exp in EXP:
    print(f'Working on VAR {VAR}, lat_sel {lat_sel}, exp {exp}')

    year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, exp)
    print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')

    ModelList = loc.read_model_list(dir_inputs, MIP, exp, VAR, False)

    print(ModelList)

    dimMod = len(ModelList.Model)
    time_all = np.arange(year_min, year_max ) + 0.5
    dimt = len(time_all)
    print(f'Number of years: {dimt}')

    name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

    da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                      dims=['model', 'time'])

    for i in range(dimMod):
        
        files = loc.select_files(MIP, exp, VAR, ModelList.loc[i], verbose)
        ds = open_files(files)
        
        print(ds)
        
        ds = harmonize_lat_lon_names(ds)
        
        if np.any(ds['lon']<0):
            ds = change_longitude_values(ds)
        
        # Remove all useless data
        #condition 2D latitude/longitude
        
        if (len(ds['lon'].shape)) == 2:
            print('Two dimensional lon/lat arrays')
        
            mask_lon = (ds['lon'] > 280) & (ds['lon'] < 355)
            mask_lat = (ds['lat'] >= 20) & (ds['lat'] <= 40)
        
            cropped_ds = ds.where(mask_lon & mask_lat, drop=True)
            
        elif (len(ds['lon'].shape)) == 1:
            print('One dimension lon/lat arrays')
            
            cropped_ds = ds.sel(lon=slice(280,355))
            
            for lat_sel in LAT_SEL:
                cropped_ds = cropped_ds.sel(lat=lat_sel, method='nearest')
                zonal_sum = cropped_ds.sum(axis='lon')
                
        
        
        
        print(cropped_ds)