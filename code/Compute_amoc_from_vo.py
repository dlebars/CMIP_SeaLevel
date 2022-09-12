###############################################################################
# Compute_amoc_from_vo.py:
# Compute time series of AMOC at different latitudes using the meridional 
# velocity fields from CMIP5 amd CMIP6
#
# A few assumptions are made about the grid:
# 1- Grid cells are equally spaces in the zonal direction at the LAT_SEL latitudes
# 2- Zonal direction is approximately in line with y indexes even for 2D 
# curvilinear grids
#
# Run time can take a while because of the interpollation (8mn per model)
# On linux use nohup, and python -u to avoid output buffering :
# nohup python -u Compute_amoc_from_vo.py > out_Compute_amoc_from_vo_cmip6.txt &
###############################################################################

import itertools
import os
import time

import numpy as np
import xarray as xr

import mod_loc as loc
    
verbose = True # Print additional information
LAT_SEL = [26, 35] # Latitude to select
MIP = 'cmip6'

lon_min, lon_max = 280, 355
lat_min, lat_max = 20, 40 # Only used to crop data

if MIP == 'cmip5':
    VAR = 'vo'
    EXP = ['historical', 'rcp26', 'rcp45','rcp85']
    
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

def harmonize_lat_lon_lev_names(ds):
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
            
    if 'olevel' in ds.coords:
        ds = ds.rename({'olevel':'lev'})
        
    if ('x' not in ds.coords):
        if 'i' in ds.coords:
            ds = ds.rename({'i':'x'})
        if 'nlon' in ds.coords:
            ds = ds.rename({'nlon':'x'})
            
    if ('bnds' not in ds.dims):
        if 'd2' in ds.dims:
            ds = ds.rename({'d2':'bnds'})
        if 'axis_nbounds' in ds.dims:
            ds = ds.rename({'axis_nbounds':'bnds'})
    
    if type(ds)==xr.core.dataset.Dataset:
        if ('lev_bnds' not in ds.data_vars):
            if 'lev_bounds' in ds.data_vars:
                ds = ds.rename({'lev_bounds':'lev_bnds'})
            
    return ds

def change_longitude_values(ds): 
        
    print('Changing longitude from -180,0 to 180,360')
    
    new_lon = np.where(ds['lon'].values<0, ds['lon'].values+360, ds['lon'].values)
    ds['lon'].values = new_lon
    
    return ds

def compute_amoc(lat_sec, lev_bnds_in):
    '''Compute the Antactic Meridional Overturning Circulation from meridional 
    velocities along a section'''

    RE = 6.371e6 # Radius of the Earth
    
    dlon = (lat_sec['lon'][1:].values-lat_sec['lon'][:-1].values)
    dlon = np.append(dlon, dlon[-1])
    zonal_length = np.cos(np.radians(lat_sec.lat.values))*RE*np.radians(dlon)

    # Remove time dependence of depth
    if 'time' in lev_bnds_in.coords:
        lev_bnds = lev_bnds_in.isel(time=0)
    else:
        lev_bnds = lev_bnds_in

    vertical_length = lev_bnds.isel(bnds=1)-lev_bnds.isel(bnds=0)
    
    print('Multiplying velocities and section area')
    vo_vol = lat_sec*zonal_length*vertical_length
    
    print('Computing sum along longitude')
    zonal_sum = vo_vol.sum(dim='lon')
    
    print('Computing vertical sum')
    vertical_cumsum = zonal_sum.cumsum(dim='lev')

    amoc = vertical_cumsum.max(dim='lev')
    
    return amoc

for exp in EXP:
    print(f'Working on VAR {VAR}, exp {exp}')

    year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, exp)
    print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')

    ModelList = loc.read_model_list(dir_inputs, MIP, exp, VAR, False)
    
    # Subselect models
    ModelList = ModelList.iloc[4:] # 2 and 3 don't work...

    print(ModelList)

    dimMod = len(ModelList.Model)
    time_all = np.arange(year_min, year_max ) + 0.5
    dimt = len(time_all)
    print(f'Number of years: {dimt}')

    for i in range(dimMod):
        print(' ')
        print(f'###### Working on {ModelList.Model.iloc[i]} ##################')
        print(' ')
        start_time = time.time()
        
        da = xr.DataArray(np.zeros([1, dimt, len(LAT_SEL)]), 
                          coords=[[ModelList.Model.iloc[i]], time_all, LAT_SEL], 
                          dims=['model', 'time', 'latitude'])
        
        files = loc.select_files(MIP, exp, VAR, ModelList.iloc[i], verbose)
        ds = open_files(files)
        
        print(ds)
        
        loc_da = ds[VAR].squeeze()
        
        loc_da = harmonize_lat_lon_lev_names(loc_da)
        ds = harmonize_lat_lon_lev_names(ds)
        
        if np.any(loc_da['lon']<0):
            loc_da = change_longitude_values(loc_da)
        
        if (len(loc_da['lon'].shape)) == 2:
            print('### Two dimensional lon/lat arrays ###')
        
            mask_lon = (loc_da['lon'] > lon_min) & (loc_da['lon'] < lon_max)
            mask_lat = (loc_da['lat'] >= lat_min) & (loc_da['lat'] <= lat_max)
        
            cropped_da = loc_da.where(mask_lon & mask_lat, drop=True)
            
            avg_lat = cropped_da.lat.mean(axis=1)
            std_lat = cropped_da.lat.std(axis=1)
            
            print(' ')
            print('Is it fine to average the latitude values?')
            print('Standard deviation of latitude along the x direction:')
            print(std_lat.values)
            print(' ')
            
            for indl, lat_sel in enumerate(LAT_SEL):
                print(f'Latitude {lat_sel}')
                ind_lat = np.abs(avg_lat-lat_sel).argmin().values
                lat_c = cropped_da[:,:,ind_lat].copy()
                
                print('Computing yearly average')
                lat_c = loc.yearly_mean(lat_c)
                lat_c = lat_c.swap_dims({'x': 'lon'})
                
                print('Computing amoc')
                amoc = compute_amoc(lat_c, ds['lev_bnds'])
                
                try:
                    da[0,:,indl] = amoc.sel(time=slice(year_min,year_max))
                except:
                    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    print('There seem to be missing time data, model not used')

            
        elif (len(loc_da['lon'].shape)) == 1:
            print('### One dimension lon/lat arrays ###')
            
            cropped_da = loc_da.sel(lon=slice(lon_min, lon_max))
            
            for indl, lat_sel in enumerate(LAT_SEL):
                lat_c = cropped_da.sel(lat=lat_sel, method='nearest')
                
                lat_c = loc.yearly_mean(lat_c)
                
                amoc = compute_amoc(lat_c, ds['lev_bnds'])
                
                try:
                    da[i,:,indl] = amoc.sel(time=slice(year_min,year_max))
                except:
                    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    print('There seem to be missing time data, model not used')
    
        da = da/1e6 # Convert m3/s to Sv

        da.attrs['long_name'] = 'AMOC volume transport'
        da.attrs['description'] = ('AMOC volume transport computed fom integrating the'+ 
                                   'meridional velocity')
        da.attrs['units'] = 'Sv'
        ds = xr.Dataset()
        ds['amoc'] = da

#         ### Remove missing models
#         ds = ds.where(ds['amoc']!=0)
#         ds = ds.dropna('model',how='all')

        ds = ds.sel(time=slice(year_min,year_max))

        print("### Export data to a NetCDF file ######################################")
        script_name = os.path.basename(__file__)
        name_output = f'{dir_outputs}{MIP}_amoc_{VAR}_{exp}_{ModelList.Model.iloc[i]}_{year_min}_{year_max-1}.nc'
        loc.export2netcdf(ds, name_output, script_name)
        
        print(' ')
        print('Run time for this model:')
        print(f'{(time.time() - start_time)/60} minutes')
        print(' ')