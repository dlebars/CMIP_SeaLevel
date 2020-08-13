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