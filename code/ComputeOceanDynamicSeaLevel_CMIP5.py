###############################################################################
# ComputeOceanDynmicSeaLevel_CMIP5.py: 
# - Read zos variable from local CMIP5 files synchronized from ESGF nodes using 
# synda
# - Correct the data by removing the trend from piCcontrol simulations
# - Regrid all fields to a common 1x1 lat/lon grid
# - Export the result as a NetCDF file
# Equivalent to the former PrepThermalExpMapsTS.ncl script
###############################################################################

from datetime import datetime
import os

import numpy as np
import xarray as xr
import pandas as pd
import xesmf as xe

import mod_loc as loc
import mod_trend_picontrol as pic

verbose = False
VAR = 'zos'
EXP = 'historical' # historical, rcp45, rcp85

ref_p_min = 1986  # Included. Beginning of reference period
ref_p_max = 2006  # Excluded. End of reference period

if EXP == 'historical':
    year_min = 1900 # Could start from 1850
    year_max = 2006
else:
    year_min = 2006
    year_max = 2101 
    
dir_outputs = '../outputs/'
dir_inputs = '../inputs/'

if EXP == 'historical':
    EXPm = 'rcp85'
else:
    EXPm = EXP

ModelList = pd.read_csv(dir_inputs+'CMIP5modelSelection_'+EXPm+'_'+VAR+'.txt', 
                        delim_whitespace=True, names=['Center','Model'], 
                        comment='#')
Model = ModelList.Model

###### Start and end of each period 
years_s = np.arange(year_min,year_max) + 0.5

#Read the regular 1*1 grid to use for regridded outputs
mask_ds = xr.open_dataset(dir_inputs+'reference_masks.nc')

weights = np.cos(np.deg2rad(mask_ds.lat))
weights.name = 'weights'

# Make a dataset to regrid with xESMF
ds_out = xr.Dataset({'lat': (['lat'], mask_ds.lat),
                     'lon': (['lon'], mask_ds.lon)})

print('Model used:')
print(Model)

for i in range(len(Model)):
    print(f'####### Working on model {i}, {Model[i]}  ######################')
    files_hist = loc.select_cmip5_files('historical', VAR, ModelList.Center[i], 
                                        Model[i])
    if EXP != 'historical':
        files_sce = loc.select_cmip5_files(EXP, VAR, ModelList.Center[i], 
                                           Model[i])
    if verbose:
        print('#### Using the following historical files: ####')
        print(files_hist)
        if EXP != 'historical':
            print('#### Using the following scenario files: ####')
            print(files_sce)

    try:
        hist_ds = xr.open_mfdataset(files_hist,combine='by_coords')
        if EXP != 'historical':
            sce_ds = xr.open_mfdataset(files_sce,combine='by_coords')
    except:
        print(f'!!!!!!!!! Could not open data from {Model[i]}!!!!!!!!!!!!!!!')
        continue
    
    if EXP != 'historical':
        ds = xr.concat([hist_ds,sce_ds],'time')
    else:
        ds = hist_ds
        
    y_ds = loc.yearly_mean(ds)
    
    if len(y_ds.lat.shape) == 1:
        name_lat = 'lat'
        name_lon = 'lon'
    elif len(y_ds.lat.shape) == 2:
        name_lat = 'rlat'
        name_lon = 'rlon'        
    
    if 'i' and 'j' in y_ds.coords:
        y_ds = y_ds.rename({'j':name_lat, 'i':name_lon})
    
    try:
        reg_method = 'bilinear'
        regridder = xe.Regridder(y_ds, ds_out, reg_method, periodic=True)
    except:
        try:
            reg_method = 'nearest_s2d'
            regridder = xe.Regridder(y_ds, ds_out, reg_method, periodic=True)
        except:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(f'Regridding did not work for {Model[i]}')
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
        conv_pic_hist = float(y_ds.time[0]) - float(hist_ds.attrs['branch_time'])
    except:
        # Pick a random large value that makes sure branching is not used in
        # trend_zos_pic_cmip5
        conv_pic_hist = -9999
    
    Trend_pic_coeff = loc.trend_zos_pic_cmip5(ModelList.iloc[i], order=1, 
                                              year_min=1850, year_max=2100,
                                              conv_pic_hist=conv_pic_hist)
    
    # Build polynomial from coefficients and convert from m to cm
    Trend_pic = xr.polyval(coord=y_ds.time, coeffs=Trend_pic_coeff)*100
    
    # Remove the average over the reference period
    Trend_pic = Trend_pic - Trend_pic.sel(time=slice(ref_p_min,ref_p_max)
                                         ).mean(dim='time')
    
    Trend_pic = Trend_pic.rename({Trend_pic.dims[1]:name_lat, 
                              Trend_pic.dims[2]:name_lon})
    
    MAT_CorrectedZOS_reg = np.zeros([len(years_s), len(mask_ds.lat), len(mask_ds.lon)])
    
    ##### Loop on the years ######################################
    for idx, year in enumerate(years_s):
        print(f'Working on year: {year}')

        VAR1 = y_ds[VAR].sel(time=year)

        if (Model[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                          'MRI-CGCM3']): 
            VAR1 = np.where(VAR1==0,np.nan,VAR1)

        AnomVAR1 = (VAR1 - RefVAR1)*100 # Convert m to cm
        AnomVAR1 = AnomVAR1*MaskRefVAR1

        # Effective number of years to detrend: year of interest 
        # minus mean of reference period
        nbyears = year+0.5 - (ref_p_max+ref_p_min)/2

        DTrendVAR1 = AnomVAR1 - Trend_pic.sel(time=year)

        # Regrid to the reference 1*1 degree grid            
        DTrendVAR1_reg = regridder(DTrendVAR1)
    
        # Mask other problematic regions here
        if Model[i] in ['MIROC5', 'GFDL-ESM2M', 'GFDL-CM3','GISS-E2-R', 
                         'GISS-E2-R-CC']:
            DTrendVAR1_reg  = DTrendVAR1_reg*mask_ds.mask_med
        else:
            DTrendVAR1_reg  = DTrendVAR1_reg*mask_ds.mask

        # Fill the NaN values
        DTrendVAR1_reg = DTrendVAR1_reg.interpolate_na('lon')

        DTrendVAR1_reg  = DTrendVAR1_reg*mask_ds.mask

        area_mean       = DTrendVAR1_reg.weighted(weights).mean(('lon', 'lat'))
        print(f'Removing area mean of:{np.round(area_mean.values,2)} cm')
        MAT_CorrectedZOS_reg[idx,:,:] = DTrendVAR1_reg - area_mean

    regridder.clean_weight_file()

    ### Export to a NetCDF file
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg, 
                                        coords=[years_s, mask_ds.lat, mask_ds.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Model[i]]},0)
    MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})

    MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                       'ComputeOceanDynmicSeaLevel_CMIP5.py')
    MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')
    MAT_OUT_ds.attrs['emission_scenario'] = EXP

    NameOutput = f'{dir_outputs}CMIP5_{VAR}_{EXP}_{Model[i]}_{year_min}_{year_max}.nc'
    if os.path.isfile(NameOutput):
        os.remove(NameOutput)
    MAT_OUT_ds.to_netcdf(NameOutput)
