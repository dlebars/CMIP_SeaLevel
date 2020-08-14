# mod_loc.py
from pathlib import Path
from datetime import datetime
import os

import numpy as np
import xarray as xr
from scipy import signal

def select_cmip5_files(VAR, EXP, Center, Model):
    '''Return a list of paths to the CMIP5 data files'''
    DataDir  = '/nobackup/users/bars/synda/cmip5/output1/'
    p = Path(DataDir+Center+'/'+Model+
                 '/'+EXP+'/'+'mon')
    files = list(p.glob('*/*/*/*/'+VAR+'/*'+VAR+'*.nc'))
    # Select the last version of data: 
    vs = []
    for k in range(len(files)):
        part = files[k].parts
        vs.append(part[len(part)-3])
    vs.sort()
    files = sorted(p.glob('*/*/*/'+vs[-1]+'/'+VAR+'/*'+VAR+'*.nc'))
    return files

def yearly_mean(ds):
    '''Convert the data set or data array to year'''
    try:
        ds.coords['year'] = ds.time.dt.year
    except:
        years = np.array([ds.time[i].dt.year.values.item() for i in range(len(ds.time))])
        ds.coords['year'] = xr.DataArray(years, dims=['time'])
    y_ds   = ds.groupby('year').mean(dim='time')
    return y_ds
    
def trend_zos_pic_cmip5(ModelList, order, year_min, year_max, conv_pic_hist, 
                        verbose=False):
    '''Compute zos trend over the pre-industrial control model simulations'''
    # - Could it work for zos as well? It would make the ComputeGlobalThermalExpansion 
    #scripts simpler
    # - Possibility to compare linear and 2nd order detrending?
    # Works also on different experiments? No.
    # CMIP5 and CMIP6?
    
    VAR = "zos" # To remove if the script doesn't work for zostoga
    EXP = "piControl"
    tot_year = year_max - year_min + 1
    
    Models = ModelList.Models
    files = select_cmip5_files(VAR, 'piControl', ModelList.Centers, Models)

    if verbose:
        print("#### Using following files: ####")
        print(files)
    
    try:
        ds = xr.open_mfdataset(files,combine='by_coords')
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'Could not open data from {Models[i]}')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

    y_ds = yearly_mean(ds)
    new_year = np.array(y_ds.year) + conv_pic_hist
    overlap_years = len(np.where((new_year >= year_min) & (new_year <= year_max))[0])
    print(f'Number of overlapping years : {overlap_years}')
    
    # Require that at least 90% of the years are available
    if overlap_years >= tot_year*0.9:
        print('Using branching time')
        y_ds = y_ds.assign_coords(year=(y_ds.year + conv_pic_hist))
    else:
        print('Not using branching time for piControl')
        # Assumes piControl simulation starts in 1850
        y_ds = y_ds.assign_coords(year=(y_ds.year- y_ds.year[0]+ 1850.5))
        
    VAR1 = y_ds[VAR].sel(year=slice(year_min, year_max))
    VAR1_coeff = VAR1.polyfit(dim='year',deg=order)
    
    return VAR1_coeff.polyfit_coefficients

def rotate_longitude(ds):
    '''Rotate the longitude of an xarray dataset from [0,360] to [-180,180]'''
    
    if 'lon' in ds.dims:
        lon = 'lon'
    elif 'longitude' in ds.dims:
        lon = 'longitude'
    else:
        print('Longitude not found in dimensions')
    
    ds = ds.roll({lon:180}, roll_coords=True)
    ds[lon] = np.where(ds[lon]>180, ds[lon]-360, ds[lon])
    return ds

def full_data_paths(DataDir, MIP, ModelList, EXP, VAR):
    '''Provides a list of full paths to CMIP6 data'''
    DataPath = (DataDir+MIP+'/'+ModelList.Center+'/'+ModelList.Model+
                '/'+EXP+'/'+ModelList.Ensemble+'/Omon/'+VAR+'/'+
                ModelList.Grid+'/'+ModelList[EXP+'_Version'])
    print('Looking for files there:')
    print(DataPath)
    p = Path(DataPath)
    all_files = sorted(p.glob('*'+VAR+'*.nc'))
    return all_files

def export2netcdf(da, name_output, script_name):
    '''Convert a DataArray to a Dataset, add some general metadata, make sure 
    remove any file with the name_output and export as netcdf'''
    
    ds = xr.Dataset({da.name: da})
    ds.attrs['source_file'] = f'This NetCDF file was built from {script_name}'
    ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

    if os.path.isfile(name_output):
        os.remove(name_output)
    ds.to_netcdf(name_output)
    
def print_results(time_all, AVAR1c, AVAR1ct):
    '''Compute mean and standard deviations of detrended and non-detrended 
    time series'''
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