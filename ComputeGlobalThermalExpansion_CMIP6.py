###############################################################################
# ComputeGlobalMeanThermalExpansion_CMIP6.py: From CMIP6 data
###############################################################################

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal
import os
from datetime import datetime
import sys

### Functions definition ######################################################

def full_data_paths(DataDir, MIP, ModelList, EXP, VAR):
    '''Provides a list of full paths to data'''
    DataPath = (DataDir+MIP+'/'+ModelList.Center+'/'+ModelList.Model+
                '/'+EXP+'/'+ModelList.Ensemble+'/Omon/'+VAR+'/'+
                ModelList.Grid+'/'+ModelList[EXP+'_Version'])
    print('Looking for files there:')
    print(DataPath)
    p = Path(DataPath)
    all_files = sorted(p.glob('*'+VAR+'*.nc'))
    return all_files

###############################################################################

VAR = 'zostoga'
EXP = ['ssp585','piControl'] # 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
MIP = ['ScenarioMIP','CMIP'] # MIP to which the above runs are from

ref_p_min = 1986
ref_p_max = 2005

year_min = 1986  # Included
year_max = 2100  # Excluded (2101)

DataDir  = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'
Dir_SelectPath = './SelectPaths_CMIP6/'

ModelList = pd.read_csv(Dir_SelectPath+'AvailableExperiments_'+str(VAR)+
                        '_historical_piControl_'+str(EXP[0])+'.csv')

dimMod = len(ModelList.Model)
dimt = year_max-year_min
print(dimt)
time_all = np.arange(year_min, year_max)+0.5
print(len(time_all))

AVAR1  = np.zeros([2,dimMod,dimt])

for j in range(2):
    for i in range(dimMod):
        if (MIP[j] == 'ScenarioMIP') and  (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
            # For this model the scenrios are done at DKRZ while piControl 
            # and historical are done at MPI-M
            ModelList.Center[i] = 'DKRZ'
        files1 = full_data_paths(DataDir, MIP[j], ModelList.iloc[i], EXP[j], VAR)

        if j==0:
            # Add historical simulation as well
            if (ModelList.Model[i] == 'MPI-ESM1-2-HR'):
                ModelList.Center[i] = 'MPI-M'
            files12 = full_data_paths(DataDir, 'CMIP', ModelList.iloc[i], 
                                      'historical', VAR)
            files1 = files12 + files1
        if len(files1) > 0:
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

        print(VAR1)
        time1  = f1.time
        dimt1  = len(time1)
        try:
            timeUT = time1.dt.year
        except:
            # Mix of time format, dt needs to be applied element wise
            timeUT = np.array([time1[i].dt.year.values.item() for i in range(len(time1))])
            timeUT = xr.DataArray(timeUT, coords=[timeUT], dims=['time'])

        # Convert from month to year
        try:
            VAR1.coords['year'] = VAR1.time.dt.year
        except:
            years = np.array([VAR1.time[i].dt.year.values.item() for i in range(len(VAR1))])
            VAR1.coords['year'] = xr.DataArray(years, dims=['time'])
        VAR1a   = VAR1.groupby('year').mean(dim='time')
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

        AVAR1[j,i,indt2] = VAR1a[indt]
        
        print('AVAR1[j,i,indt2]')
        print(AVAR1[j,i,indt2])

if year_max == 2300:
    print('### List of models that run to 2300')
    for i in range(0,dimMod-1):
        tot_mis = np.sum(np.isnan(AVAR1[0,i,:]))
        if tot_mis <=50:
            print(ModelList.Model[i])

### Adjust for reference period
AVAR1c = np.zeros_like(AVAR1)
indref = np.where((time_all > ref_p_min) & (time_all < ref_p_max))[0]
for j in range(2):
    for i in range(dimMod):
        AVAR1c[j,i,:] = AVAR1[j,i,:] - AVAR1[j,i,indref].mean()

### Correct for the trend in the piControl simulation
AVAR1ct = np.zeros_like(AVAR1)
AVAR1ct[1,:,:] = signal.detrend(AVAR1[1,:,:], axis=1, type='linear')
Trend_piControl = AVAR1[1,:,:] - AVAR1ct[1,:,:]
AVAR1ct[0,:,:] = AVAR1[0,:,:] - Trend_piControl
              
# Re-correct for the reference period
for j in range(2):
    for i in range(dimMod):
        AVAR1ct[j,i,:] = AVAR1ct[j,i,:] - AVAR1ct[j,i,indref].mean()

### Compute mean and standard deviations of detrended and non-detrended time
# series
AVAR1c_m   = np.mean(AVAR1c,axis=1)
AVAR1c_sd  = np.std(AVAR1c,axis=1)
AVAR1c_95p = AVAR1c_m + 1.64*AVAR1c_sd
AVAR1c_05p = AVAR1c_m - 1.64*AVAR1c_sd

AVAR1ct_m   = np.mean(AVAR1ct,axis=1)
AVAR1ct_sd  = np.std(AVAR1ct,axis=1)
AVAR1ct_95p = AVAR1ct_m + 1.64*AVAR1ct_sd
AVAR1ct_05p = AVAR1ct_m - 1.64*AVAR1ct_sd
              
### Ouput values
indc = np.where(time_all == 2099.5)[0][0]
print("Mean and 5-95 percentile range before piC detrend: ")
print("Year 2099")
print(str(AVAR1c_m[0,indc])+' [ '+str(AVAR1c_05p[0,indc])+' - '
      +str(AVAR1c_95p[0,indc])+' ]')

indr = np.where((time_all >= 2081) & (time_all <= 2100))[0]
print('Year 2081-2099')
print(str(AVAR1c_m[0,indr].mean())+' [ '+str(AVAR1c_05p[0,indr].mean())+' - '
      +str(AVAR1c_95p[0,indr].mean())+' ]')

print("Mean and 5-95 percentile range after piC detrend: ")
print("Year 2099")
print(str(AVAR1ct_m[0,indc])+' [ '+str(AVAR1ct_05p[0,indc])+' - '
      +str(AVAR1ct_95p[0,indc])+' ]')

print("Year 2081-2099")
print(str(AVAR1ct_m[0,indr].mean())+' [ '+str(AVAR1ct_05p[0,indr].mean())+' - '
      +str(AVAR1ct_95p[0,indr].mean())+' ]')

print("### Export data to a NetCDF file ######################################")

# Build a DataSet
da = xr.DataArray(AVAR1ct, coords=[ EXP, ModelList.Model, time_all], 
                  dims=['experiment', 'model', 'time'])
MAT_OUT_ds = xr.Dataset({VAR+'_detrended': da})

MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                   'ComputeGlobalMeanThermalExpansion_CMIP6.py')

MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

NameOutput = 'CMIP6_SeaLevel_'+EXP[0]+'_'+VAR+'_'+str(year_min)+'-'+str(year_max)+'.nc'
if os.path.isfile(NameOutput):
    os.remove(NameOutput)
MAT_OUT_ds.to_netcdf(NameOutput) #mode='a' to append or overwrite
