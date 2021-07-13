###############################################################################
# ComputeGlobalThermalExpansion.py
###############################################################################

import os
import sys

import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = True # Print additional information
VAR = 'zostoga'
MIP = 'cmip6'
# EXP available:
# cmip6: 'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'historical', 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'ssp245'
dir_outputs = '../outputs/'
dir_inputs = '../inputs/'

year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, EXP)
print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')
print(f'using this reference period: {ref_p_min}-{ref_p_max-1}, including {ref_p_max-1}')

gap = 0.01 # Maximum gap authorized (in meters) when removing discontinuities

# Make sure that the years from the reference period are read and saved
year_min_min = min(year_min, ref_p_min)

def open_files(file_names):
    try:
        ds = xr.open_mfdataset(file_names, combine='by_coords', use_cftime=True)
    except:
        print(f'Open by_coords did not work for:{file_names}')
        print('Using nested option instead')
        try:
            ds = xr.open_mfdataset(hist_files, combine='nested', 
                                   concat_dim='time', 
                                   use_cftime=True)
        except:
            print('ERROR: Nested option did not work either')

    return ds

ModelList = loc.read_model_list(dir_inputs, MIP, EXP, VAR)
dimMod = len(ModelList.Model)
time_all = np.arange(year_min_min, year_max ) + 0.5
dimt = len(time_all)
print(f'Number of years: {dimt}')

ds = xr.Dataset()
name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                  dims=['model', 'time'])
trend_da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
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
        
    VAR1 = all_ds[VAR].squeeze()
    VAR1a = loc.yearly_mean(VAR1)

    # Remove discontinuites in some time series
    try:
        VAR1a = loc.remove_discontinuities(VAR1a, gap)
    except ValueError as err:
        print(err.args)
        continue
    
    # Compute the trend from the piControl simulations and save trend
    if verbose:
        pic.info_branching(hist_ds.attrs)
    
    try:
        # Convert the year from piControl to historical run
        if MIP == 'cmip5':
            conv_pic_hist = float(VAR1a.time[0]) - float(hist_ds.attrs['branch_time'])
        elif MIP == 'cmip6':
            attrs = {'units': hist_ds.attrs['parent_time_units']}
            time_flt = [float(hist_ds.attrs['branch_time_in_parent'])]
            time_ds = xr.Dataset({'time': ('time', time_flt, attrs)})
            time_ds = xr.decode_cf(time_ds, use_cftime=True)
            conv_pic_hist = float(VAR1a.time[0]) - time_ds.time.dt.year.values[0]
    except:
        # Pick a random large value that makes sure branching is not used in
        # trend_zos_pic_cmip5
        conv_pic_hist = -9999
        
    Trend_pic_coeff, branching_method = pic.trend_pic(
        MIP, VAR, ModelList.iloc[i], order=1, year_min=1850, year_max=2100,
        conv_pic_hist=conv_pic_hist, gap=gap, rmv_disc=True, verbose=verbose)
    
    # Build polynomial from coefficients
    Trend_pic = xr.polyval(coord=VAR1a.time, coeffs=Trend_pic_coeff)

    trend_da[i,:] = Trend_pic.sel(time=slice(year_min_min,year_max))
    da[i,:] = VAR1a.sel(time=slice(year_min_min,year_max))
    
ds[VAR+'_corrected'] = da
ds['trend_picontrol'] = trend_da

### Remove missing models
ds = ds.where(ds.zostoga_corrected!=0)
ds = ds.dropna('model',how='all')

if year_max == 2300: #TODO Old code bellow
    print('### List of models that run to 2300')
    for i in range(0,dimMod-1):
        tot_mis = np.sum(np.isnan(AVAR1[0,i,:]))
        if tot_mis <=50:
            print(ModelList.Model[i])

### Adjust for reference period
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

if verbose and EXP!='historical':
    print('### Before detrending:')
    loc.print_results_da(ds[VAR+'_corrected'])

### Correct for the trend in the piControl simulation
ds[VAR+'_corrected'] = ds[VAR+'_corrected'] -  ds['trend_picontrol']

# Re-correct for the reference period
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

ds = ds.sel(time=slice(year_min,year_max))

if verbose and EXP!='historical':
    print('### After detrending:')    
    loc.print_results_da(ds[VAR+'_corrected'])

print("### Export data to a NetCDF file ######################################")
script_name = os.path.basename(__file__)
name_output = f'{dir_outputs}{MIP}_{VAR}_{EXP}_{year_min}_{year_max-1}.nc'
loc.export2netcdf(ds, name_output, script_name)
