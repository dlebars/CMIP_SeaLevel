# mod_loc.py
from pathlib import Path
from datetime import datetime
import os
import sys

import numpy as np
import xarray as xr
import pandas as pd
from scipy import signal

def select_cmip5_files(EXP, VAR, ModelList):
    '''Return a list of paths to the CMIP5 data files'''
    
    if EXP == 'rcp60':
        data_dir = '/nobackup/users/bars/synda_data_bck/cmip5/output1/'
        path_string = f'*/*/*/*/*{VAR}*.nc'
        path_nb = 2
    else:
        data_dir = '/nobackup/users/bars/synda/cmip5/output1/'
        path_string = f'*/*/*/*/{VAR}/*{VAR}*.nc'
        path_nb = 3
    p = Path(data_dir+ModelList.Center+'/'+ModelList.Model+
             '/'+EXP+'/'+'mon')
    files = list(p.glob(path_string))
    # Select the last version of data: 
    vs = []
    for k in range(len(files)):
        part = files[k].parts
        vs.append(part[len(part)-path_nb])
    vs.sort()
    if EXP == 'rcp60':
        files = sorted(p.glob(f'*/*/*/{vs[-1]}/*{VAR}*.nc'))
    else:
        files = sorted(p.glob(f'*/*/*/{vs[-1]}/{VAR}/*{VAR}*.nc'))
        
    return files

def select_cmip6_files(EXP, VAR, ModelList):
    '''Provides a list of full paths to CMIP6 data'''
    
    MIP = {'piControl' : 'CMIP', 
           'historical' : 'CMIP',
           'ssp119' : 'ScenarioMIP', 
           'ssp126' : 'ScenarioMIP', 
           'ssp245' : 'ScenarioMIP', 
           'ssp370' : 'ScenarioMIP', 
           'ssp585' : 'ScenarioMIP'}
    
    realm = {'zos' : 'O', # O for Ocean and A for atmosphere
             'zostoga' : 'O',
             'ps' : 'A',
             'uas' : 'A',
             'vas' : 'A'}
    
    data_dir  = '/nobackup/users/bars/synda_data/CMIP6/'
    data_path = (data_dir+MIP[EXP]+'/'+ModelList.Center+'/'+ModelList.Model+
                '/'+EXP+'/'+ModelList[EXP+'_Variant']+'/'+realm[VAR]+'mon/'+VAR+'/'+
                ModelList.Grid+'/'+ModelList[EXP+'_Version'])
    print('Looking for files there:')
    print(data_path)
    p = Path(data_path)
    all_files = sorted(p.glob('*'+VAR+'*.nc'))
    
    return all_files

def yearly_mean(ds):
    '''Convert the data set or data array to year'''
    
    try:
        ds.coords['year'] = ds.time.dt.year
    except:
        years = np.array([ds.time[i].dt.year.values.item() 
                          for i in range(len(ds.time))])
        ds.coords['year'] = xr.DataArray(years, dims=['time'])
        
    y_ds = ds.groupby('year').mean(dim='time')
    y_ds = y_ds.rename({'year':'time'})
    # Center the time axis to the middle of the year
    y_ds = y_ds.assign_coords(time=(y_ds.time+0.5))
    
    return y_ds

def rotate_longitude(ds, name_lon):

    ds = ds.assign_coords({name_lon:(((ds[name_lon] + 180 ) % 360) - 180)})
    ds = ds.sortby(ds[name_lon])

    return ds

def export2netcdf(ds, name_output, script_name):
    '''Export a dataset as netcdf. 
    If the input is a DataArray, convert to Dataset. 
    Add some general metadata and remove any file with the name_output before
    exporting'''
    
    if isinstance(ds, xr.DataArray):
        ds = xr.Dataset({ds.name: ds})
    ds.attrs['source_file'] = f'This NetCDF file was built from {script_name}'
    ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

    if os.path.isfile(name_output):
        os.remove(name_output)
    ds.to_netcdf(name_output)
    
def print_results(time_all, AVAR1c):
    '''Compute mean and standard deviations of detrended and non-detrended 
    time series'''
    
    AVAR1c_m   = np.mean(AVAR1c,axis=1)
    AVAR1c_sd  = np.std(AVAR1c,axis=1)
    AVAR1c_95p = AVAR1c_m + 1.64*AVAR1c_sd
    AVAR1c_05p = AVAR1c_m - 1.64*AVAR1c_sd

    ### Ouput values
    indc = np.where(time_all == 2099.5)[0][0]
    print("Mean and 5-95 percentile range: ")
    print("Year 2099")
    print(str(AVAR1c_m[0,indc])+' [ '+str(AVAR1c_05p[0,indc])+' - '
          +str(AVAR1c_95p[0,indc])+' ]')

    indr = np.where((time_all >= 2081) & (time_all <= 2100))[0]
    print('Year 2081-2099')
    print(str(AVAR1c_m[0,indr].mean())+' [ '+str(AVAR1c_05p[0,indr].mean())+' - '
          +str(AVAR1c_95p[0,indr].mean())+' ]')

def print_results_da(da):
    '''Compute mean and standard deviations of detrended and non-detrended 
    time series'''
    
    AVAR1c_m   = da.mean(dim='model')
    AVAR1c_sd  = da.std(dim='model')
    AVAR1c_95p = AVAR1c_m + 1.64*AVAR1c_sd
    AVAR1c_05p = AVAR1c_m - 1.64*AVAR1c_sd

    ### Ouput values
    print("Mean and 5-95 percentile range: ")
    print("Year 2099")
    print(f'{float(AVAR1c_m.sel(time=2099.5))} [ '+
          f'{float(AVAR1c_05p.sel(time=2099.5))} - '+
          f'{float(AVAR1c_95p.sel(time=2099.5))} ]')

    print('Year 2081-2099')  
    print(f'{float(AVAR1c_m.sel(time=slice(2081,2101)).mean())} [ '+
          f'{float(AVAR1c_05p.sel(time=slice(2081,2101)).mean())} - '+
          f'{float(AVAR1c_95p.sel(time=slice(2081,2101)).mean())} ]')

def remove_discontinuities(da, gap):
    '''Remove discontinuities in a time series, numpy or data array.
    da: The input data
    gap: the maximum gap allowed in the data above which the 
    discontinuity is removed'''
    
    #da_out = da.copy() -> Not clear why this was here
    if isinstance(da, xr.DataArray):
        # Make sure to load the data, Dask arrays do not support item assigment
        da.load()
        diff = da.diff('time')
    elif isinstance(da, np.ndarray):
        diff = np.array(da[1:]) - np.array(da[:-1])
    else:
        print('ERROR: Input object type not supported')
        
    indpb = np.where(np.abs(diff) > gap)[0]
    
    if len(indpb) > 10:
        raise ValueError('Error', 'Too many many discontinuities, not using this model')
    
    if len(indpb) > 0:
        print("### Removing discontinuities at these indices: ####")
        print(indpb)
        for k in indpb:
            da[k+1:] = da[k+1:] - da[k+1] + da[k]
    
    return da

def start_end_ref_dates(MIP, EXP):
    
    # This is the reference period of AR5. It can easily be changed later 
    # during the data analysis. 
    ref_p_min = 1986  # Included. Beginning of reference period
    ref_p_max = 2006  # Excluded. End of reference period

    year_start_sce = {'cmip5': 2006, 'cmip6': 2015}
    # year_end_sce is excluded.
    # 2101 works for CMIP5, not for some models of CMIP6
    year_end_sce = {'cmip5': 2101, 'cmip6': 2101}

    if EXP == 'historical':
        year_min = 1850 # Could start from 1850? Used to be 1900
        year_max = year_start_sce[MIP]
    else:
        year_min = year_start_sce[MIP]
        year_max = year_end_sce[MIP]
    
    return year_min, year_max, ref_p_min, ref_p_max

def read_model_list(dir_inputs, MIP, EXP, VAR):
    '''Reads the list of models to use for the analysis'''
    
    if MIP == 'cmip5':
        col_names = ['Center','Model']
        ModelList = pd.read_csv(f'{dir_inputs}CMIP5modelSelection_{EXP}_{VAR}.txt', 
                                delim_whitespace=True, names=col_names,
                                comment='#')
    elif MIP == 'cmip6':
        dir_SelectPath = '../SelectPaths_CMIP6/'
        if EXP in ['piControl', 'historical']:
             ModelList = pd.read_csv(f'{dir_SelectPath}AvailableExperiments_{VAR}'+
                                    f'_historical_piControl.csv')       
        else:
            ModelList = pd.read_csv(f'{dir_SelectPath}AvailableExperiments_{VAR}'+
                                    f'_{EXP}_historical_piControl.csv')
            
    return ModelList

def select_files(MIP, EXP, VAR, ModelList_loc, verbose=False):
    '''Select file names and paths'''
    
    if MIP == 'cmip5':
        files = select_cmip5_files(EXP, VAR, ModelList_loc)
        
    elif MIP == 'cmip6': 
        # For this model the scenarios are done at DKRZ while piControl 
        # and historical are done at MPI-M
        if (ModelList_loc.Model == 'MPI-ESM1-2-HR'):
            if EXP in ['historical', 'piControl']:
                ModelList_loc.Center = 'MPI-M'
            else:
                ModelList_loc.Center = 'DKRZ'

        files = select_cmip6_files(EXP, VAR, ModelList_loc)
            
    if verbose:
        if len(files) > 0:
            print('#### Using the following files: ####')
            [print(str(x)) for x in files]
        else:
            sys.exit('ERROR: No file available at that location') 
        
    return files