###############################################################################
# ComputeOceanDynamicSeaLevel.py: 
# - Read zos variable from local CMIP5/CMIP6 files synchronized from ESGF nodes using 
# synda
# - Correct the data by removing the trend from piCcontrol simulations
# - Regrid all fields to a common 1x1 lat/lon grid
# - Export the result as a NetCDF file
# Equivalent to the former PrepThermalExpMapsTS.ncl script
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
import pandas as pd
import xesmf as xe

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = True
VAR = 'zos'
MIP = 'cmip6' # cmip5 or cmip6
# EXP available:
# cmip6: 'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
# cmip5: 'historical', 'rcp26', 'rcp45', 'rcp60','rcp85'
EXP = 'ssp585'
trend_order = 1 # Order of the polynomial fit used to detrend the data based on
                # the piControl simulation

year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, EXP)
#year_min = 2097 # Specify a shorter time range for tests
#year_max = 2099
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
        hist_ds = hist_ds.load()
        
        if EXP != 'historical':
            sce_ds = xr.open_mfdataset(sce_files, combine='by_coords', 
                                       use_cftime=True)
            sce_ds = sce_ds.load()
    except:
        print(f'!!!!!!!!! Could not open data from {Model.iloc[i]}!!!!!!!!!!!!!!!')
        print('Try the function open_mfdataset with the option combine="nested" ')
        continue
    
    if EXP != 'historical':
        ds = xr.concat([hist_ds,sce_ds],'time')
    else:
        ds = hist_ds
        
    y_ds = loc.yearly_mean(ds)
    
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
        # Used to take this filename as input:
        #filename=f'{reg_method}_{Model.iloc[i]}_{VAR}_{EXP}')
    except:
        try:
            reg_method = 'nearest_s2d'
            regridder = xe.Regridder(y_ds, ds_out, reg_method, periodic=True)
        except:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(f'Regridding did not work for {Model.iloc[i]}')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            continue
    
    if verbose:
        print(regridder)
    
    RefVAR1 = y_ds[VAR].sel(time=slice(ref_p_min,ref_p_max)).mean(dim='time')

    # Use the reference field to mask the small seas that are not connected to the
    # ocean and areas where sea ice is included on the ocean load
    RefVAR1_corr = RefVAR1 - RefVAR1.mean() #TODO weigthed mean
    MaskRefVAR1  = np.where((RefVAR1_corr>=2) | (RefVAR1_corr<=-2),np.nan,1)

    if verbose:
        pic.info_branching(hist_ds.attrs)
            
    try:
        # Convert the year from piControl to historical run
        if MIP == 'cmip5':
            conv_pic_hist = float(y_ds.time[0]) - float(hist_ds.attrs['branch_time'])
        elif MIP == 'cmip6':
            attrs = {'units': hist_ds.attrs['parent_time_units']}
            time_flt = [float(hist_ds.attrs['branch_time_in_parent'])]
            time_ds = xr.Dataset({'time': ('time', time_flt, attrs)})
            time_ds = xr.decode_cf(time_ds, use_cftime=True)
            conv_pic_hist = float(y_ds.time[0]) - time_ds.time.dt.year.values[0]
    except:
        # Pick a random large value that makes sure branching is not used in
        # trend_pic
        conv_pic_hist = -9999

    Trend_pic_coeff, branching_method = pic.trend_pic(
        MIP, VAR, ModelList.iloc[i], order=trend_order, year_min=1850, 
        year_max=2100,conv_pic_hist=conv_pic_hist, gap=None, rmv_disc=False, 
        verbose=verbose)
    
    try:
        # This breaks when Trend_pic_coeff does not contain values.
        # It hapens for BCC-CSM2-MR for which polyfit does not return 
        # coefficients but does not crash
        test = Trend_pic_coeff.values
        
        # Build polynomial from coefficients
        Trend_pic = xr.polyval(coord=y_ds.time, coeffs=Trend_pic_coeff)

        # Remove the average over the reference period
        Trend_pic = Trend_pic - Trend_pic.sel(time=slice(ref_p_min,ref_p_max)
                                             ).mean(dim='time')
    except:
        print('!!! WARNING: Detrending from piControl for this model does not'+ 
              ' work !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        branching_method = 'no_detrending'
        Trend_pic = xr.DataArray(np.zeros(len(years)), coords=[years], dims=["time"])
    
    MAT_CorrectedZOS_reg = np.zeros([len(years), len(mask_ds.lat), len(mask_ds.lon)])
    
    ##### Loop on the years ######################################
    for idx, year in enumerate(years):
        print(f'Working on year: {year}')

        VAR1 = y_ds[VAR].sel(time=year)

        if (Model.iloc[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                          'MRI-CGCM3']): 
            VAR1 = np.where(VAR1==0,np.nan,VAR1)

        AnomVAR1 = (VAR1 - RefVAR1)*MaskRefVAR1

        DTrendVAR1 = AnomVAR1 - Trend_pic.sel(time=year)

        # Regrid to the reference 1*1 degree grid
        DTrendVAR1_reg = regridder(DTrendVAR1)

        # Mask problematic regions here
        if Model.iloc[i] in ['MIROC5', 'GFDL-ESM2M', 'GFDL-CM3','GISS-E2-R', 
                         'GISS-E2-R-CC']:
            DTrendVAR1_reg  = DTrendVAR1_reg*mask_ds.mask_med
        else:
            DTrendVAR1_reg  = DTrendVAR1_reg*mask_ds.mask

        # Fill the NaN values
        DTrendVAR1_reg = DTrendVAR1_reg.interpolate_na('lon', method='nearest', 
                                                       fill_value='extrapolate')
        DTrendVAR1_reg = DTrendVAR1_reg*mask_ds.mask
        area_mean = DTrendVAR1_reg.weighted(weights).mean(('lon', 'lat'))
        print(f'Removing spatial mean of:{np.round(area_mean.values,2)} m')
        MAT_CorrectedZOS_reg[idx,:,:] = DTrendVAR1_reg - area_mean

#    regridder.clean_weight_file()

    ### Export to a NetCDF file
    # Convert from m to cm
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg*100, 
                                        coords=[years, mask_ds.lat, mask_ds.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Model.iloc[i]]},0)
    MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    MAT_CorrectedZOS_reg.attrs['branching_method'] = branching_method
    MAT_CorrectedZOS_reg.attrs['detrending_order'] = f'{trend_order}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})

    MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                       'ComputeOceanDynmicSeaLevel_CMIP5.py')
    MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')
    MAT_OUT_ds.attrs['emission_scenario'] = EXP

    NameOutput = f'{dir_outputs}{MIP}_{VAR}_{EXP}_{Model.iloc[i]}_{year_min}_{year_max-1}.nc'
    
    if os.path.isfile(NameOutput):
        os.remove(NameOutput)
    MAT_OUT_ds.to_netcdf(NameOutput)
