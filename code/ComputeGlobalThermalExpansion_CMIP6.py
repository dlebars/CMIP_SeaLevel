###############################################################################
# ComputeGlobalMeanThermalExpansion_CMIP6.py: From CMIP6 data
###############################################################################

from pathlib import Path
import os
import sys

import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal

import mod_loc as loc

verbose = True
VAR = 'zostoga'
EXP = ['ssp585','piControl'] # 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
MIP = ['ScenarioMIP','CMIP'] # MIP to which the above runs belong

ref_p_min = 1986
ref_p_max = 2005 #TODO excluded so increase to 2006?

year_min = 1986  # Included
year_max = 2100  # Included

DataDir  = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'
Dir_SelectPath = '../SelectPaths_CMIP6/'
dir_outputs = '../outputs/'
ModelList = pd.read_csv(Dir_SelectPath+'AvailableExperiments_'+str(VAR)+
                        '_historical_piControl_'+str(EXP[0])+'.csv')

dimMod = len(ModelList.Model)
time_all = np.arange(year_min, year_max ) + 0.5
dimt = len(time_all)
print(dimt)

ds = xr.Dataset()
name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

for j in range(2):
    da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                      dims=['model', 'time'])
    trend_da = 
    for i in range(dimMod):
        print(f'####### Working on model {i}, {Model[i]}  ######################')
        
        if (MIP[j] == 'ScenarioMIP') and  (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
            # For this model the scenarios are done at DKRZ while piControl 
            # and historical are done at MPI-M
            ModelList.Center[i] = 'DKRZ'
            
        files1 = loc.full_data_paths(DataDir, MIP[j], ModelList.iloc[i], EXP[j], VAR)

        if j==0:
            # Add historical simulation as well
            if (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
                ModelList.Center[i] = 'MPI-M'
                
            files12 = loc.full_data_paths(DataDir, 'CMIP', ModelList.iloc[i], 
                                      'historical', VAR)
            files1 = files12 + files1
            
        if len(files1) > 0:
            if verbose:
                print('#### Using the following files: ####')
                [print(str(x)) for x in  files1]
        else:
            sys.exit('ERROR: No file available at that location')
            
        try:
            f1 = xr.open_mfdataset(files1, combine='by_coords', use_cftime=True)
        except:
            print('Open by_coords did not work for:'+ ModelList.Center[i]+
                  '/'+ModelList.Model[i]+'/'+EXP[j])
            print('Using nested option instead')
            f1 = xr.open_mfdataset(files1, combine='nested', concat_dim='time', 
                                   use_cftime=True)
        
        VAR1   = f1[VAR].squeeze()
        VAR1a = loc.yearly_mean(VAR1)

        if EXP[j] == 'piControl':
            # Assumes piControl simulation starts in 1850
            # This is not the case for all models so it should be improved!
            VAR1a = VAR1a.assign_coords(time=(VAR1a.time- VAR1a.time[0]+ 1850.5))
    
        if ModelList.Model[i] == 'MRI-ESM2-0':
            VAR1a = loc.remove_discontinuities(VAR1a, 0.02)
        
        da[i,:] = VAR1a.sel(time=slice(year_min,year_max))
    ds[name_da[j]] = da

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
dtrend = signal.detrend(ds['trend_piControl'], axis=1, type='linear')
ds['trend_piControl'] = ds['trend_piControl'] - dtrend
ds[VAR+'_corrected'] = ds[VAR+'_corrected'] - ds['trend_piControl']

# Re-correct for the reference period
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

if verbose:
    print('### After detrending:')    
    loc.print_results_da(ds[VAR+'_corrected'])

print("### Export data to a NetCDF file ######################################")
script_name = os.path.basename(__file__)
name_output = f'{dir_outputs}CMIP6_SeaLevel_{EXP[0]}_{VAR}_{year_min}_{year_max}.nc'
loc.export2netcdf(ds, name_output, script_name)
