###############################################################################
# PreparePlaneVariables.py: 
# - Read 2D space variable from local CMIP5/CMIP6 files synchronized from ESGF 
# nodes using synda
# - Optionally correct the data by removing the trend from piCcontrol simulations
# - Regrid data to a common 1x1 lat/lon grid
# - Export the result as a NetCDF file
#
# Run time can take a while because of the interpollation, around 10 hours for
# a scenario with 20 model
# On linux use nohup, and python -u to avoid output buffering :
# nohup python -u ComputeOceanDynamicSeaLevel.py > out_rcp60.txt &
###############################################################################

from datetime import datetime
import os

import numpy as np
import xarray as xr
import xesmf as xe

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = True
VAR = 'zos' # 'zos', 'ps', 'uas', 'vas'
MIP = 'cmip6' # cmip5 or cmip6
# EXP available:
# cmip6: 'piControl', 'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'piControl', 'historical', 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'historical'

detrend = True # Detrend using piControl simulation (does not work for piControl)
trend_order = 1 # Order of the polynomial fit used to detrend the data based on
                # the piControl simulation

dir_outputs = '/nobackup/users/bars/CMIP6_regridded/' #'../outputs/'
dir_inputs = '../inputs/'

####### End of user defined parameters ########################################

sce_list = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585', 
            'rcp26', 'rcp45', 'rcp60','rcp85']

# Determine wether the anomaly compared to a reference period needs to be taken
anom_dic = {'zos' : True,
            'ps' : False,
            'uas' : False,
            'vas' : False}

ModelList = loc.read_model_list(dir_inputs, MIP, EXP, VAR)

# Remove a model for which the analysis fails
ModelList = ModelList.loc[ModelList.Model!='CNRM-CM6-1-HR']
# Remove a model that is on an unstructured grid for which regridding fails
ModelList = ModelList.loc[ModelList.Model!='AWI-CM-1-1-MR']

Model = ModelList.Model

# Read the regular 1*1 grid to use for regridded outputs
mask_ds = xr.open_dataset(dir_inputs+'reference_masks.nc')

# The mask closes the Mediteranean sea which is not necessary so here I modify
# it
mask_ds['mask'].loc[dict(lat=slice(34,36), lon=354.5)] = 1

mask_ds = loc.rotate_longitude(mask_ds, 'lon')

weights = np.cos(np.deg2rad(mask_ds.lat))
weights.name = 'weights'

# Make a dataset to regrid with xESMF
ds_out = xr.Dataset({'lat': (['lat'], mask_ds.lat),
                     'lon': (['lon'], mask_ds.lon)})

print('Model used:')
print(Model)

for i in range(len(Model)):
    print(f'####### Working on model {i}, {Model.iloc[i]}  ######################')

    if EXP in ['piControl', 'historical']:
        files = loc.select_files(MIP, EXP, VAR, ModelList.iloc[i], verbose)
        y_ds = loc.open_files(files)

    elif EXP in sce_list:
        hist_files = loc.select_files(MIP, 'historical', VAR, ModelList.iloc[i], verbose)
        sce_files = loc.select_files(MIP, EXP, VAR, ModelList.iloc[i], verbose)
        
        hist_y_ds = loc.open_files(hist_files)
        sce_y_ds = loc.open_files(sce_files)
        y_ds = xr.concat([hist_y_ds,sce_y_ds],'time')
        
    else:
        print(f'ERROR: EXP {EXP} not supported')

    if Model.iloc[i] == 'BCC-CSM2-MR':
        y_ds = y_ds.rename({'lat':'rlat', 'lon':'rlon'})
    
    if ('lat' not in y_ds.coords) and ('lon' not in y_ds.coords):
        if 'latitude' and 'longitude' in y_ds.coords:
            y_ds = y_ds.rename({'latitude':'lat', 'longitude':'lon'})
        
        elif 'nav_lat' and 'nav_lon' in y_ds.coords:
            y_ds = y_ds.rename({'nav_lat':'lat', 'nav_lon':'lon'})
    
    # Build array of years
    # For piControl it is read from input data since models use different time 
    #references.
    # For historical and scenarios it is fixed to have uniform output size that
    #helps data analysis.
    
    if EXP=='piControl':
        years = y_ds.time.values
        year_min = int(years[0]) # Only used for output name
        year_max = int(years[-1])+1
        ref_p_min = year_min
        ref_p_max = year_min+20
        
    else:
        year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, EXP)
        years = np.arange(year_min,year_max) + 0.5

    print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')

    if anom_dic[VAR]:
        print(f'using this reference period: {ref_p_min}-{ref_p_max-1}, including {ref_p_max-1}')
    
    # Build regridder with xESMF
    try:
        reg_method = 'bilinear'
        regridder = xe.Regridder(y_ds, ds_out, reg_method, periodic=True)

    except:
        try:
            reg_method = 'nearest_s2d'
            regridder = xe.Regridder(y_ds, ds_out, reg_method, periodic=True)
        except:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(f'Regridding did not work for {Model.iloc[i]}')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            continue

    if verbose:
        print(f'Using {reg_method} for regridding')
        #print(regridder) # Use to debug
    
    if anom_dic[VAR]:
        ref_da = y_ds[VAR].sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

    if VAR=='zos':
        # Use the reference field to mask the small seas that are not connected to the
        # ocean and areas where sea ice is included on the ocean load
        ref_da_corr = ref_da - ref_da.mean()
        ref_da_mask  = np.where((ref_da_corr>=2) | (ref_da_corr<=-2),np.nan,1)

    if detrend:
        if EXP == 'historical':
            attrs = y_ds.attrs
        else:
            attrs = hist_y_ds.attrs
            
        Trend_pic, branching_method = pic.trend_pic_ts(
            y_ds, attrs, MIP, VAR, ModelList.iloc[i], trend_order, 
            rmv_disc=False, verbose=verbose)
        
        # Remove the average over the reference period
        Trend_pic = Trend_pic - Trend_pic.sel(time=slice(ref_p_min,ref_p_max)
                                             ).mean(dim='time')
    
    da_full = y_ds[VAR]
    
    if Model.iloc[i] == 'FGOALS-g3':
        # The historical file is a bit too long for this model
        da_full = da_full.drop_duplicates(dim='time', keep='last')
    
    
    if (Model.iloc[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                          'MRI-CGCM3']): 
        da_full = np.where(da_full==0,np.nan,da)
    
    if anom_dic[VAR]:
        da_full = da_full - ref_da
        
    if VAR=='zos':
        da_full = da_full*ref_da_mask
        
    if detrend:
        da_full = da_full - Trend_pic.sel(time=da_full.time)
        
    # Regrid to the reference 1*1 degree grid
    reg_da = regridder(da_full)

    if VAR=='zos':
        # Mask problematic regions here
        if Model.iloc[i] in ['MIROC5', 'GFDL-ESM2M', 'GFDL-CM3','GISS-E2-R', 
                         'GISS-E2-R-CC']:
            reg_da  = reg_da*mask_ds.mask_med
        else:
            reg_da  = reg_da*mask_ds.mask

        # Fill the NaN values
        reg_da = reg_da.interpolate_na('lon', method='nearest', 
                                       fill_value='extrapolate')
        reg_da = reg_da*mask_ds.mask

        # Remove spatial mean
        area_mean = reg_da.weighted(weights).mean(('lon', 'lat'))
        print(f'Removing spatial mean of:{np.round(area_mean.values,2)} m')
        reg_da = reg_da - area_mean
            
    if len(reg_da.time) == len(years):
        MAT_CorrectedZOS_reg = reg_da
    elif len(reg_da.time) > len(years):
        print('ERROR: There is an issue with the number of years, probably some'+ 
              'duplicates to remove')
    elif len(reg_da.time) < len(years):
        print('WARNING: There are some missing years.'+
              f' Data has {reg_da.time} years, should be {len(years)}'+
              'interpolating the missing years')
        MAT_CorrectedZOS_reg = reg_da.intep(time=years)

    print("### Export data to a NetCDF file ######################################")
    
    # Build a data array from the numpy array
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg, 
                                        coords=[years, mask_ds.lat, mask_ds.lon], 
                                        dims=['time', 'lat', 'lon'])
    if VAR=='zos':
        # Convert from m to cm
        MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg*100
        MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
        MAT_CorrectedZOS_reg.attrs['long_name'] = 'Ocean dynamic sea level'
    elif VAR=='ps':
        MAT_CorrectedZOS_reg.attrs['units'] = 'Pa'
        MAT_CorrectedZOS_reg.attrs['long_name'] = 'Surface Air Pressure'
    elif VAR=='uas':
        MAT_CorrectedZOS_reg.attrs['units'] = 'm s-1'
        MAT_CorrectedZOS_reg.attrs['long_name'] = 'Eastward Near-Surface Wind'
    elif VAR=='vas':
        MAT_CorrectedZOS_reg.attrs['units'] = 'm s-1'
        MAT_CorrectedZOS_reg.attrs['long_name'] = 'Northward Near-Surface Wind' 
    else:
        print(f'ERROR: Variable {VAR} not supported')
    
    if anom_dic[VAR]:
        MAT_CorrectedZOS_reg.attrs['ref_period'] = (
            f'The folowing reference period {ref_p_min}-{ref_p_max-1},'+
            f' including {ref_p_max-1}, was used to compute anomalies')
    
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Model.iloc[i]]},0)
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    MAT_CorrectedZOS_reg.attrs['variant'] = ModelList[f'{EXP}_Variant'].iloc[i]
    
    if detrend:
        MAT_CorrectedZOS_reg.attrs['branching_method'] = (
            'This dataset was detrended using piControl trend')
        MAT_CorrectedZOS_reg.attrs['branching_method'] = branching_method
        MAT_CorrectedZOS_reg.attrs['detrending_order'] = f'{trend_order}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})
    MAT_OUT_ds.attrs['emission_scenario'] = EXP
    
    script_name = os.path.basename(__file__)
    name_output = f'{dir_outputs}{MIP}_{VAR}_{EXP}_{Model.iloc[i]}_{year_min}_{year_max-1}.nc'
    loc.export2netcdf(MAT_OUT_ds, name_output, script_name)

    
    
