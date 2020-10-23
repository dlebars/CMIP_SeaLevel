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
# EXP available:
# cmip6: 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'rcp45'

# Select the mip that corresponds to the scenario
MIP_dic = {'ssp119':'cmip6',
           'ssp126':'cmip6',
           'ssp245':'cmip6',
           'ssp370':'cmip6',
           'ssp585':'cmip6', 
           'rcp26':'cmip5', 
           'rcp45':'cmip5',
           'rcp60':'cmip5',
           'rcp85':'cmip5'}
MIP = MIP_dic[EXP]

ref_p_min = 1986
ref_p_max = 2006

year_min = 1986  # Included
year_max = 2100  # Included

gap = 0.02 # Maximum gap authorized (in meters) when removing discontinuities

dir_outputs = '../outputs/'
dir_inputs = '../inputs/'

# Select the file containing the model list to analyse
if MIP == 'cmip5':
    col_names = ['Center','Model']
    ModelList = pd.read_csv(f'{dir_inputs}CMIP5modelSelection_{EXP}_{VAR}.txt', 
                            delim_whitespace=True, names=col_names,
                            comment='#')
elif MIP == 'cmip6':
    dir_SelectPath = '../SelectPaths_CMIP6/'
    ModelList = pd.read_csv(f'{dir_SelectPath}AvailableExperiments_{VAR}'+
                            f'_historical_piControl_{EXP}.csv')

dimMod = len(ModelList.Model)
time_all = np.arange(year_min, year_max ) + 0.5
dimt = len(time_all)
print(dimt)

ds = xr.Dataset()
name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                  dims=['model', 'time'])
trend_da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                  dims=['model', 'time'])

for i in range(dimMod):
    print(f'####### Working on model {i}, {ModelList.Model[i]}  ############')

    # Read paths and file names
    if MIP == 'cmip5':
        sce_files = loc.select_cmip5_files(EXP, VAR, ModelList.loc[i])
        hist_files = loc.select_cmip5_files('historical', VAR, ModelList.loc[i])
        
    elif MIP == 'cmip6':
        if ModelList.Model[i] == 'MPI-ESM1-2-HR':
            # For this model the scenarios are done at DKRZ while piControl 
            # and historical are done at MPI-M
            ModelList.Center[i] = 'DKRZ'

        sce_files = loc.select_cmip6_files(EXP, VAR, ModelList.iloc[i])

        # Read historical simulation as well
        if (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
            ModelList.Center[i] = 'MPI-M'

        hist_files = loc.select_cmip6_files('historical', VAR, ModelList.iloc[i])

    if len(sce_files+hist_files) > 0:
        if verbose:
            print('#### Using the following files: ####')
            [print(str(x)) for x in  (sce_files+hist_files)]
    else:
        sys.exit('ERROR: No file available at that location')

    # Open files
    try:
        hist_ds = xr.open_mfdataset(hist_files, combine='by_coords', use_cftime=True)
        sce_ds = xr.open_mfdataset(sce_files, combine='by_coords', use_cftime=True)
    except:
        print('Open by_coords did not work for:'+ ModelList.Center[i]+
              '/'+ModelList.Model[i]+'/'+EXP)
        print('Using nested option instead')
        hist_ds = xr.open_mfdataset(hist_files, combine='nested', concat_dim='time', 
                               use_cftime=True)
        sce_ds = xr.open_mfdataset(sce_files, combine='nested', concat_dim='time', 
                               use_cftime=True)
    
    all_ds = xr.concat([hist_ds,sce_ds],'time')
        
    VAR1 = all_ds[VAR].squeeze()
    VAR1a = loc.yearly_mean(VAR1)

    # Remove discontinuites in some time series
    VAR1a = loc.remove_discontinuities(VAR1a, gap)
    
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

    trend_da[i,:] = Trend_pic.sel(time=slice(year_min,year_max))
    da[i,:] = VAR1a.sel(time=slice(year_min,year_max))
    
ds[VAR+'_corrected'] = da
ds['trend_picontrol'] = trend_da

if year_max == 2300: #TODO Old code bellow
    print('### List of models that run to 2300')
    for i in range(0,dimMod-1):
        tot_mis = np.sum(np.isnan(AVAR1[0,i,:]))
        if tot_mis <=50:
            print(ModelList.Model[i])

### Adjust for reference period
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

if verbose:
    print('### Before detrending:')
    loc.print_results_da(ds[VAR+'_corrected'])

### Correct for the trend in the piControl simulation
ds[VAR+'_corrected'] = ds[VAR+'_corrected'] -  ds['trend_picontrol']

# Re-correct for the reference period
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

if verbose:
    print('### After detrending:')    
    loc.print_results_da(ds[VAR+'_corrected'])

print("### Export data to a NetCDF file ######################################")
script_name = os.path.basename(__file__)
name_output = f'{dir_outputs}{MIP}_SeaLevel_{EXP}_{VAR}_{year_min}_{year_max}.nc'
loc.export2netcdf(ds, name_output, script_name)
