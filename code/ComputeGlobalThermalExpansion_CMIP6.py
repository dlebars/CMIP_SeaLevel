###############################################################################
# ComputeGlobalMeanThermalExpansion_CMIP6.py: From CMIP6 data
###############################################################################

import os
import sys

import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = True
VAR = 'zostoga'
EXP = 'ssp585' # 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'

ref_p_min = 1986
ref_p_max = 2005 #TODO excluded so increase to 2006?

year_min = 1986  # Included
year_max = 2100  # Included

Dir_SelectPath = '../SelectPaths_CMIP6/'
dir_outputs = '../outputs/'
ModelList = pd.read_csv(Dir_SelectPath+'AvailableExperiments_'+str(VAR)+
                        '_historical_piControl_'+EXP+'.csv')

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

    if ModelList.Model[i] == 'MPI-ESM1-2-HR':
        # For this model the scenarios are done at DKRZ while piControl 
        # and historical are done at MPI-M
        ModelList.Center[i] = 'DKRZ'

    files1 = loc.select_cmip6_files(EXP, VAR, ModelList.iloc[i])

    # Add historical simulation as well
    if (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
        ModelList.Center[i] = 'MPI-M'

    files12 = loc.select_cmip6_files('historical', VAR, ModelList.iloc[i])
    files1 = files12 + files1

    if len(files1) > 0:
        if verbose:
            print('#### Using the following files: ####')
            [print(str(x)) for x in  files1]
    else:
        sys.exit('ERROR: No file available at that location')

    try:
        ds1 = xr.open_mfdataset(files1, combine='by_coords', use_cftime=True)
    except:
        print('Open by_coords did not work for:'+ ModelList.Center[i]+
              '/'+ModelList.Model[i]+'/'+EXP)
        print('Using nested option instead')
        ds1 = xr.open_mfdataset(files1, combine='nested', concat_dim='time', 
                               use_cftime=True)

    VAR1 = ds1[VAR].squeeze()
    VAR1a = loc.yearly_mean(VAR1)

    if ModelList.Model[i] == 'MRI-ESM2-0':
        VAR1a = loc.remove_discontinuities(VAR1a, 0.02)
        
    # Compute the trend from the piControl simulations and save trend
    if verbose:
        pic.info_branching(ds1.attrs)
    
    try:
        # Convert the year from piControl to historical run
        conv_pic_hist = float(ds1.time[0]) - float(ds1.attrs['branch_time_in_parent'])
    except:
        # Pick a random large value that makes sure branching is not used in
        # trend_zos_pic_cmip5
        conv_pic_hist = -9999
        
    Trend_pic_coeff = pic.trend_pic(VAR, ModelList.iloc[i], order=1, 
                                          year_min=1850, year_max=2100,
                                          conv_pic_hist=conv_pic_hist)
    
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
name_output = f'{dir_outputs}CMIP6_SeaLevel_{EXP}_{VAR}_{year_min}_{year_max}.nc'
loc.export2netcdf(ds, name_output, script_name)
