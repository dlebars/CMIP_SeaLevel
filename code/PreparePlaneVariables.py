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
# cmip6: 'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'historical', 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'ssp585'

detrend = True # Detrend using piControl simulation
trend_order = 1 # Order of the polynomial fit used to detrend the data based on
                # the piControl simulation

year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, EXP)

print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')
print(f'using this reference period: {ref_p_min}-{ref_p_max-1}, including {ref_p_max-1}')

dir_outputs = '../outputs/'
dir_inputs = '../inputs/'

ModelList = loc.read_model_list(dir_inputs, MIP, EXP, VAR)

# Remove a model for which the analysis fails
ModelList = ModelList.loc[ModelList.Model!='CNRM-CM6-1-HR']
# Remove a model that is on an unstructured grid for which regridding fails
ModelList = ModelList.loc[ModelList.Model!='AWI-CM-1-1-MR']

Model = ModelList.Model

# Build array of years
years = np.arange(year_min,year_max) + 0.5

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
    
    hist_files = loc.select_files(MIP, 'historical', VAR, ModelList.iloc[i], verbose)
    
    if EXP != 'historical':
        sce_files = loc.select_files(MIP, EXP, VAR, ModelList.iloc[i], verbose)       
    
    # Open files
    try:
        hist_ds = xr.open_mfdataset(hist_files, combine='by_coords', 
                                    use_cftime=True)
        hist_y_ds = loc.yearly_mean(hist_ds)
        hist_y_ds = hist_y_ds.load()
        
        if EXP != 'historical':
            sce_ds = xr.open_mfdataset(sce_files, combine='by_coords', 
                                       use_cftime=True)
            sce_y_ds = loc.yearly_mean(sce_ds)
            sce_y_ds = sce_y_ds.load()
            
    except:
        print(f'!!!!!!!!! Could not open data from {Model.iloc[i]}!!!!!!!!!!!!!!!')
        print('Try the function open_mfdataset with the option combine="nested" ')
        continue
 
    if EXP != 'historical':
        y_ds = xr.concat([hist_y_ds,sce_y_ds],'time')
    else:
        y_ds = hist_y_ds

    if Model.iloc[i] == 'BCC-CSM2-MR':
        y_ds = y_ds.rename({'lat':'rlat', 'lon':'rlon'})
    
    if 'latitude' and 'longitude' in y_ds.coords:
        y_ds = y_ds.rename({'latitude':'lat', 'longitude':'lon'})
    elif 'nav_lat' and 'nav_lon' in y_ds.coords:
        y_ds = y_ds.rename({'nav_lat':'lat', 'nav_lon':'lon'})
    
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
    
    ref_da = y_ds[VAR].sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

    if VAR=='zos':
        # Use the reference field to mask the small seas that are not connected to the
        # ocean and areas where sea ice is included on the ocean load
        ref_da_corr = ref_da - ref_da.mean()
        ref_da_mask  = np.where((ref_da_corr>=2) | (ref_da_corr<=-2),np.nan,1)

    if detrend:
        Trend_pic, branching_method = pic.trend_pic_ts(
            y_ds, hist_ds.attrs, MIP, VAR, ModelList.iloc[i], trend_order, 
            rmv_disc=False, verbose=verbose)
        
        # Remove the average over the reference period
        Trend_pic = Trend_pic - Trend_pic.sel(time=slice(ref_p_min,ref_p_max)
                                             ).mean(dim='time')
    
    MAT_CorrectedZOS_reg = np.zeros([len(years), len(mask_ds.lat), len(mask_ds.lon)])
    
    ##### Loop on the years ######################################
    for idx, year in enumerate(years):
        print(f'Working on year: {year}')

        da = y_ds[VAR].sel(time=year)

        if (Model.iloc[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                          'MRI-CGCM3']): 
            da = np.where(da==0,np.nan,da)

        anom_da = da - ref_da
        
        if VAR=='zos':
            anom_da = anom_da*ref_da_mask

        if detrend:
            anom_da = anom_da - Trend_pic.sel(time=year)

        # Regrid to the reference 1*1 degree grid
        reg_da = regridder(anom_da)

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
            
        MAT_CorrectedZOS_reg[idx,:,:] = reg_da

    print("### Export data to a NetCDF file ######################################")
    
    # Convert from m to cm
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg*100, 
                                        coords=[years, mask_ds.lat, mask_ds.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Model.iloc[i]]},0)
    MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    MAT_CorrectedZOS_reg.attrs['variant'] = ModelList[f'{EXP}_Variant'].iloc[i]
    
    if detrend:
        MAT_CorrectedZOS_reg.attrs['branching_method'] = branching_method
        MAT_CorrectedZOS_reg.attrs['detrending_order'] = f'{trend_order}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})
    MAT_OUT_ds.attrs['emission_scenario'] = EXP
    
    script_name = os.path.basename(__file__)
    name_output = f'{dir_outputs}{MIP}_{VAR}_{EXP}_{Model.iloc[i]}_{year_min}_{year_max-1}.nc'
    loc.export2netcdf(MAT_OUT_ds, name_output, script_name)

    
    