# mod_loc.py
from pathlib import Path
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
    
def trend_zos_pic_cmip5(ModelList, order, verbose=False):
    '''Compute zos trend over the pre-industrial control model simulations'''
    # - Could it work for zos as well? It would make the ComputeGlobalThermalExpansion 
    #scripts simpler
    # - Possibility to compare linear and 2nd order detrending?
    # Works also on different experiments? No.
    # CMIP5 and CMIP6?
    
    VAR = "zos" # To remove if the script doesn't work for zostoga
    EXP = "piControl"

    year_min = 1986
    year_max = 2100

    Models   = ModelList.Models
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
    # Assumes piControl simulation starts in 1850 (some models start from 0)
    # Could be improved by checking the branch year attribute
    y_ds = y_ds.assign_coords(year=(y_ds.year- y_ds.year[0]+ 1850.5))

    VAR1 = y_ds[VAR].sel(year=slice(year_min, year_max))
    VAR1_coeff = VAR1.polyfit(dim='year',deg=order)
    
    return VAR1_coeff.polyfit_coefficients