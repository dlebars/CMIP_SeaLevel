###############################################################################
# ComputeGlobalMeanThermalExpansion_CMIP5.py: From CMIP5 data
#   Python translation of PlotThermalExp.ncl
###############################################################################

import pandas as pd
import numpy as np
import xarray as xr
from scipy import signal
import os
from datetime import datetime
import mod_loc as loc

### Functions definition ######################################################
def remove_discontinuities(da, gap):
    '''Remove discontinuities in a time series, numpy or data array.
    da: The input data, gap: the maximum gap allowed in the data above which 
    the discontinuity is removed'''
    da_out = da.copy()
    diff = np.array(da[1:]) - np.array(da[:-1])
    print(diff)
    indpb = np.where(np.abs(diff) > gap)[0]
    print("### Indices at which there are discontinuities: ####")
    print(indpb)
    for k in indpb:
        da_out[k+1:] = da[k+1:] - da[k+1] + da[k]
    return da_out

###############################################################################

VAR = 'zostoga' # zostoga, zossga, zosga
EXP = ['rcp85','piControl']

ref_p_min = 1986
ref_p_max = 2005

year_min = 1986  # Included
year_max = 2100  # Excluded

Dir_CMIP5_TE = '../../CMIP5_ThermalExp/'
Dir_outputs = '../outputs/'

col_names = ['Centers','Models']
if year_max == 2300:
    ModelList = pd.read_csv(Dir_CMIP5_TE+'CMIP5modelSelection_'+EXP[0]+'_'+VAR+
                            '_2300.txt', delim_whitespace=True, 
                            names=['Centers','Models'], comment='#')
else:
    ModelList = pd.read_csv(Dir_CMIP5_TE+'CMIP5modelSelection_'+EXP[0]+'_'+VAR+'.txt', 
                            delim_whitespace=True, names=['Centers','Models'], 
                            comment='#')

dimMod = len(ModelList.Models)
dimt = year_max-year_min
print(dimt)
time_all = np.arange(year_min, year_max)+0.5
print(len(time_all))

AVAR1  = np.zeros([2,dimMod,dimt])

for j in range(2):
    for i in range(dimMod):
        files1 = loc.select_cmip5_files(VAR, EXP[j], ModelList.Centers[i], 
                                    ModelList.Models[i])
        print('#### Using the following files: ####')
        [print(str(x)) for x in  files1]
        try:
            f1 = xr.open_mfdataset(files1,combine='by_coords')
        except:
            print('Open by_coords did not work for:'+ ModelList.Centers[i]+
                  '/'+ModelList.Models[i]+'/'+EXP[j])
            print('Using nested option instead')
            f1 = xr.open_mfdataset(files1,combine='nested', concat_dim='time')
        VAR1 = f1[VAR].squeeze()

        print(VAR1)
        time1  = f1.time
        dimt1  = len(time1)
        try:
            timeUT = time1.dt.year
        except:
            # Mix of time format, dt needs to be applied element wise
            timeUT = np.array([time1[i].dt.year.values.item() for i in range(len(time1))])
            timeUT = xr.DataArray(timeUT, coords=[timeUT], dims=['time'])

        if j == 0:
            # Add historical simulation as well
            files12 = loc.select_cmip5_files(VAR, 'historical', ModelList.Centers[i], 
                                         ModelList.Models[i])
            print('### Also using these historical files: ###')
            [print(str(x)) for x in  files1]
            f12 = xr.open_mfdataset(files12,combine='by_coords')
            timeUT11 = timeUT
            time12   = f12.time
            timeUT12 = time12.dt.year
                           
            # Need to remove overlap period
            ind12    = np.where(timeUT12 < timeUT11[0])
            print('ind12')
            print(ind12)
            VAR12     = f12[VAR].squeeze()

            timeUT = xr.concat([timeUT12[ind12],timeUT11], dim='time')
            VAR1 = xr.concat([VAR12[ind12],VAR1], dim='time')


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
            timeUTa = timeUTa - timeUTa[0] + 1850.5

        dimtl = len(timeUTa)

        if ( ModelList.Models[i] == 'bcc-csm1-1' or 
             ModelList.Models[i] == 'bcc-csm1-1-m' or
             ModelList.Models[i] == 'GISS-E2-R-CC' and 
            (VAR in ['zossga', 'zostoga']) and EXP[j] != 'piControl'):
            VAR1a = remove_discontinuities(VAR1a.values, 0.02)
    
        print('Time vector timeUTa:')
        print(timeUTa[0])
        print(timeUTa[-1])
        indt = np.where((timeUTa >= year_min) & (timeUTa < year_max))[0]
        indt2 = np.where((time_all >= (timeUTa[0].values-0.1)) & 
                         (time_all <= (timeUTa[-1].values+0.1)))[0]
        print('These two lengths should be the same:')
        print('len(indt) '+str(len(indt)))
        print('len(indt2) '+str(len(indt2)))

        if ([EXP[j], ModelList.Models[i], VAR] in 
              [['piControl', 'MIROC-ESM-CHEM', 'zossga'], 
               ['piControl', 'MIROC-ESM-CHEM', 'zostoga'], 
               ['rcp45', 'CMCC-CM', 'zossga'], 
               ['rcp85', 'HadGEM2-ES', 'zossga'], 
               ['rcp85', 'HadGEM2-ES', 'zostoga'], 
               ['rcp85', 'CNRM-CM5', 'zossga'], 
               ['rcp85', 'CNRM-CM5', 'zosga']]):
                # Missing year or double years so interpollate
                AVAR1[j,i,indt2] = np.interp(time_all[indt2], timeUTa[indt], VAR1a[indt])
        else:
            AVAR1[j,i,indt2] = VAR1a[indt]

        print('AVAR1[j,i,indt2]')
        print(AVAR1[j,i,indt2])

if year_max == 2300:
    print('### List of models that run to 2300')
    for i in range(0,dimMod-1):
        tot_mis = np.sum(np.isnan(AVAR1[0,i,:]))
        if tot_mis <=50:
            print(ModelList.Models[i])

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

### Compute mean and standard deviations of detrended and non-detrended time ##
#series
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
da = xr.DataArray(AVAR1ct, coords=[ EXP, ModelList.Models, time_all], 
                  dims=['experiment', 'model', 'time'])
MAT_OUT_ds = xr.Dataset({VAR+'_detrended': da})

MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                   'ComputeGlobalMeanThermalExpansion_CMIP5.py')

MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

NameOutput = Dir_outputs+'CMIP5_SeaLevel_'+EXP[0]+'_'+VAR+'_'+str(year_min)+'-'+str(year_max)+'.nc'
if os.path.isfile(NameOutput):
    os.remove(NameOutput)
MAT_OUT_ds.to_netcdf(NameOutput) #mode='a' to append or overwrite
