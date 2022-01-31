###############################################################################
# Compute_amoc_ts.py:
# Compute time series of AMOC at different latitudes using the msftmz and
# msftyz files
###############################################################################

import os

import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = True # Print additional information
VAR = 'msftyz' # msftmz, msftyz
lat_s = 26 # Latitude to select
MIP = 'cmip6'
# EXP available:
# cmip6: 'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'historical', 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'historical'

dir_outputs = '/nobackup/users/bars/CMIP6_regridded/'
dir_inputs = '../inputs/'

year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, EXP)
print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')

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
            print('ERROR: Nested option did not work either')

    return ds

ModelList = loc.read_model_list(dir_inputs, MIP, EXP, VAR, False)

print(ModelList)

dimMod = len(ModelList.Model)
time_all = np.arange(year_min, year_max ) + 0.5
dimt = len(time_all)
print(f'Number of years: {dimt}')

name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                  dims=['model', 'time'])

for i in range(dimMod):
    print(f'####### Working on model {i}, {ModelList.Model[i]}  ############')

    hist_files = loc.select_files(MIP, 'historical', VAR, ModelList.loc[i], verbose)
    
    if EXP != 'historical':
        sce_files = loc.select_files(MIP, EXP, VAR, ModelList.loc[i], verbose)

    hist_ds = open_files(hist_files)
    
    if EXP=='historical':
        all_ds = hist_ds
    else:
        sce_ds = open_files(sce_files)
        all_ds = xr.concat([hist_ds,sce_ds],'time')
        
    # Change the name of some variables
    if ('lat' not in all_ds.coords):
        if 'latitude' in all_ds.coords:
            all_ds = all_ds.rename({'latitude':'lat'})
        elif 'rlat' in all_ds.coords:
            all_ds = all_ds.rename({'rlat':'lat'})
        elif 'nav_lat' in all_ds.coords:
            all_ds = all_ds.rename({'nav_lat':'lat'})
        elif 'j-mean' in all_ds.coords:
            all_ds = all_ds.rename({'j-mean':'lat'})
        elif 'y' in all_ds.coords:
            all_ds = all_ds.rename({'y':'lat'})
    
    if '3basin' in all_ds.coords:
        all_ds = all_ds.rename({'3basin':'basin'})

    if 'region' in all_ds.coords:
        all_ds = all_ds.rename({'region':'sector'})
        
    if 'olevel' in all_ds.coords:
        all_ds = all_ds.rename({'olevel':'lev'})
    
    
    VAR1 = all_ds[VAR].squeeze()
    VAR1a = loc.yearly_mean(VAR1)
    print(VAR1a)
    
    try: 
        sector_names = VAR1a.sector.astype(dtype=str).load()
        sector_names = np.char.rstrip(sector_names, ' ')
        atl_ind = np.where(sector_names == 'atlantic_arctic_ocean')[0].item()
        
    except:
        
        try:
            sector_names = VAR1a.sector.astype(dtype=str).load()
            atl_ind = np.where(sector_names == 'a')[0].item()
            
        except:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('The atlantic_arctic_ocean sector was not found')
            
            if ModelList.Model[i] == 'IPSL-CM6A-LR':
                atl_ind = 1
            else:
                atl_ind = 0
                
            print(f'using index {atl_ind} instead')
            print('Check the basin variable metadata:')
            print(VAR1a.basin)
    
    if ModelList.Model[i] == 'IPSL-CM6A-LR':
        VAR1a = VAR1a.swap_dims({'y': 'lat'})
        VAR1a = VAR1a.isel(lat=slice(-1)) # Last latitude value is double
    
    VAR1a = VAR1a.isel(basin=atl_ind).sel(lat=lat_s, method='nearest')
    VAR1a = VAR1a.max(dim='lev')
    
    try:
        da[i,:] = VAR1a.sel(time=slice(year_min,year_max))
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('There seem to be missing time data, model not used')
    
if VAR == 'msftmz':
    da.attrs['long_name'] = 'Ocean Meridional Overturning Mass Streamfunction'
elif VAR == 'msftyz':
    da.attrs['long_name'] = 'Ocean Y Overturning Mass Streamfunction'

da.attrs['description'] = ('Overturning mass streamfunction arising from all advective '+
                           'mass transport processes, resolved and parameterized.')
da.attrs['units'] = 'kg s-1'
ds = xr.Dataset()
ds[VAR] = da

### Remove missing models
ds = ds.where(ds[VAR]!=0)
ds = ds.dropna('model',how='all')

ds = ds.sel(time=slice(year_min,year_max))

print("### Export data to a NetCDF file ######################################")
script_name = os.path.basename(__file__)
name_output = f'{dir_outputs}{MIP}_{VAR}_{lat_s}N_{EXP}_{year_min}_{year_max-1}.nc'
loc.export2netcdf(ds, name_output, script_name)
