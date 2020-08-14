###############################################################################
# ComputeGlobalMeanThermalExpansion_CMIP6.py: From CMIP6 data
###############################################################################

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal
import os
import sys
import mod_loc as loc

verbose = True
VAR = 'zostoga'
EXP = ['ssp585','piControl'] # 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
MIP = ['ScenarioMIP','CMIP'] # MIP to which the above runs belong

ref_p_min = 1986
ref_p_max = 2005 #TODO excluded so increase to 2006

year_min = 1986  # Included
year_max = 2100  # Excluded (2101)

DataDir  = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'
Dir_SelectPath = '../SelectPaths_CMIP6/'
dir_outputs = '../outputs/'
ModelList = pd.read_csv(Dir_SelectPath+'AvailableExperiments_'+str(VAR)+
                        '_historical_piControl_'+str(EXP[0])+'.csv')

dimMod = len(ModelList.Model)
time_all = np.arange(year_min, year_max)+0.5
dimt = len(time_all)
print(dimt)

#AVAR1  = np.zeros([2,dimMod,dimt])
ds = xr.Dataset()
name_da = {0: VAR+'_corrected', 1: 'trend_piControl'}

for j in range(2):
    da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, time_all], 
                      dims=['model', 'time'])
    for i in range(dimMod):
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
            f1 = xr.open_mfdataset(files1,combine='by_coords')
        except:
            print('Open by_coords did not work for:'+ ModelList.Center[i]+
                  '/'+ModelList.Model[i]+'/'+EXP[j])
            print('Using nested option instead')
            f1 = xr.open_mfdataset(files1,combine='nested', concat_dim='time')
        VAR1   = f1[VAR].squeeze()

        time1  = f1.time
        dimt1  = len(time1)
        try:
            timeUT = time1.dt.year
        except:
            # Mix of time format, dt needs to be applied element wise
            timeUT = np.array([time1[i].dt.year.values.item() for i in range(len(time1))])
            timeUT = xr.DataArray(timeUT, coords=[timeUT], dims=['time'])

        VAR1a = loc.yearly_mean(VAR1)

        timeUTa = VAR1a.year
        if EXP[j] == 'piControl':
            # Assumes piControl simulation starts in 1850
            # This is not the case for all models so it should be improved!
            timeUTa = timeUTa - timeUTa[0] + 1850.5

        dimtl = len(timeUTa)
    
        print('Time vector timeUTa:')
        print(timeUTa[0])
        print(timeUTa[-1])
        indt = np.where((timeUTa >= year_min) & (timeUTa < year_max))[0]
        indt2 = np.where((time_all >= (timeUTa[0].values-0.1)) & 
                         (time_all <= (timeUTa[-1].values+0.1)))[0]
        print('These two lengths should be the same:')
        print('len(indt) '+str(len(indt)))
        print('len(indt2) '+str(len(indt2)))

        #AVAR1[j,i,indt2] = VAR1a[indt]
        da[i,:] = VAR1a[indt]
    ds[name_da[j]] = da

if year_max == 2300: #TODO
    print('### List of models that run to 2300')
    for i in range(0,dimMod-1):
        tot_mis = np.sum(np.isnan(AVAR1[0,i,:]))
        if tot_mis <=50:
            print(ModelList.Model[i])

### Adjust for reference period
# AVAR1c = np.zeros_like(AVAR1)
# indref = np.where((time_all > ref_p_min) & (time_all < ref_p_max))[0]
# for j in range(2):
#     for i in range(dimMod):
#         AVAR1c[j,i,:] = AVAR1[j,i,:] - AVAR1[j,i,indref].mean()
ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean()

### Correct for the trend in the piControl simulation
# AVAR1ct = np.zeros_like(AVAR1)
# AVAR1ct[1,:,:] = signal.detrend(AVAR1[1,:,:], axis=1, type='linear')
# Trend_piControl = AVAR1[1,:,:] - AVAR1ct[1,:,:]
# AVAR1ct[0,:,:] = AVAR1[0,:,:] - Trend_piControl

AVAR1ct = np.zeros([2,dimMod, dimt])
AVAR1ct[1,:,:] = signal.detrend(ds['trend_piControl'], axis=1, type='linear')
Trend_piControl = ds['trend_piControl'] - AVAR1ct[1,:,:]
AVAR1ct[0,:,:] = ds[VAR+'_corrected'] - Trend_piControl

# detrend_da = xr.DataArray(np.zeros([dimMod, dimt]), coords=[ModelList.Model, locs], 
#                           dims=['model', 'time'])
# detrend_da = 

# Re-correct for the reference period
for j in range(2):
    for i in range(dimMod):
        AVAR1ct[j,i,:] = AVAR1ct[j,i,:] - AVAR1ct[j,i,indref].mean()
#ds = ds - ds.sel(time=slice(ref_p_min,ref_p_max)).mean()

if verbose:
    loc.print_results(time_all, AVAR1c, AVAR1ct)

print("### Export data to a NetCDF file ######################################")

script_name = os.path.basename(__file__)

# Build a Dataset from the numpy array
da = xr.DataArray(AVAR1ct, coords=[ EXP, ModelList.Model, time_all], 
                  dims=['experiment', 'model', 'time'], name=VAR+'_detrended')

name_output = f'{dir_outputs}CMIP6_SeaLevel_{EXP[0]}_{VAR}_{year_min}_{year_max}.nc'

loc.export2netcdf(da, name_output, script_name)
