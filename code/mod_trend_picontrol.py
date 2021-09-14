# mod_trend_picontol.py
import xarray as xr
import numpy as np

import mod_loc as loc

def try_attr(attr, ds_attrs):
    '''Check if an attribute is present and print its value'''
    try:
        print(f'{attr} : {ds_attrs[attr]}')
    except:
        print(f'No {attr} attribute')

def info_branching(ds_attrs):
    print('Check info about piControl branching:')
    
    # Check which MIP it is
    try:
        mip = ds_attrs['mip_era']
    except:
        mip = ds_attrs['project_id']

    try_attr('parent_experiment_id', ds_attrs)
        
    if mip == 'CMIP5':
        try_attr('parent_experiment_rip', ds_attrs)
        try_attr('branch_time', ds_attrs)
    elif mip == 'CMIP6':
        try_attr('parent_variant_label', ds_attrs)
        try_attr('branch_time_in_child', ds_attrs)
        try_attr('branch_time_in_parent', ds_attrs)
        try_attr('parent_time_units', ds_attrs)
        
def trend_pic(MIP, VAR, ModelList, order, year_min, year_max, conv_pic_hist, 
              gap, rmv_disc, verbose=False):
    '''Compute zos or zostoga trend over the pre-industrial control model simulations.
    Use rmv_disc=True to remove discontinuities in time. Only works for time 
    series (zostoga) not for 3D data (zos)'''
    
    EXP = 'piControl'
    tot_year = year_max - year_min + 1
        
    files = loc.select_files(MIP, EXP, VAR, ModelList, verbose)
    
    ds = xr.open_mfdataset(files, combine='by_coords', use_cftime=True)

    y_ds = loc.yearly_mean(ds)
    y_ds = y_ds.load()
    
    if ModelList.Model=='BCC-CSM2-MR' and VAR=='zos':
        y_ds = y_ds.rename({'lat':'rlat', 'lon':'rlon'})
    
    new_year = np.array(y_ds.time) + conv_pic_hist
    overlap_years = len(np.where((new_year >= year_min) & (new_year <= year_max))[0])
    print(f'Number of overlapping years : {overlap_years}')
    
    # Require that at least 90% of the years are available
    if overlap_years >= tot_year*0.9:
        print('Using branching time')
        branching_method = 'unsing_branching_time'
        y_ds = y_ds.assign_coords(time=(y_ds.time + conv_pic_hist))
    else:
        print('Not using branching time for piControl')
        branching_method = 'not_unsing_branching_time'
        # Assumes piControl simulation starts in 1850
        y_ds = y_ds.assign_coords(time=(y_ds.time - y_ds.time[0] + 1850.5))
        
    VAR1 = y_ds[VAR].squeeze()
    VAR1 = VAR1.sel(time=slice(year_min, year_max))
    if rmv_disc:
        VAR1 = loc.remove_discontinuities(VAR1, gap)
    VAR1_coeff = VAR1.polyfit(dim='time',deg=order)
    
    return VAR1_coeff.polyfit_coefficients, branching_method