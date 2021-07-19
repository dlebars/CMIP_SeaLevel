# mod_postpro.py
import numpy as np
import xarray as xr
import pandas as pd

def define_area(reg):
    '''Provides box coordinates given a region name'''
    
    if reg == 'dutch_coast':
        lon_min, lon_max = 3, 7
        lat_min, lat_max = 51, 54
    elif reg == 'north_sea':
        lon_min, lon_max = -2, 9
        lat_min, lat_max = 48, 60
    elif reg == 'knmi14_reg':
        lon_min, lon_max = -3.5, 7.5
        lat_min, lat_max = 51, 60
    
    return lon_min, lon_max, lat_min, lat_max

def ds2df(cmip_ds, lon_min, lon_max, lat_min, lat_max, start_year, end_year):
    '''Transform a dataset to a dataframe by averaging sea level over a region'''
    
    sel_ds = cmip_ds.sel(lon=slice(lon_min,lon_max), lat=slice(lat_min,lat_max))
    sel_da = sel_ds['CorrectedReggrided_zos'].mean(dim=['lon','lat'])

    df = pd.DataFrame(dict(time=np.arange(start_year,end_year)+0.5))
    df = df.set_index('time')

    for mod in sel_da.model.values:
        df[mod] = sel_da.sel(model=mod).drop('model').to_dataframe()
        
    return df

def read_zos_ds(data_dir, mip, sce):
    '''Read both historical and scenario datasets, select the intersecting 
    models and concate the two datasets'''
    
    hist_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zos_historical/{mip}_zos_historical_*.nc')
    sce_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zos_{sce}/{mip}_zos_{sce}_*.nc')
    model_intersection = list(set(hist_ds.model.values) & 
                              set(sce_ds.model.values))
    model_intersection.sort()
    tot_ds = xr.concat([hist_ds,sce_ds],'time').sel(model=model_intersection)
    
    return tot_ds