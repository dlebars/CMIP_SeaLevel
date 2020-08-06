###############################################################################
# ComputeOceanDynmicSeaLevel_CMIP5.py: 
# - Read zos variable from local CMIP5 files synchronized from ESGF nodes using 
# synda
# - Correct the data by removing the trend from piCcontrol simulations
# - Regrid all fields to a common 1x1 lat/lon grid
# - Export the result as a NetCDF file
# Equivalent to the former PrepThermalExpMapsTS.ncl script
###############################################################################

import numpy as np
import xarray as xr
import mod_loc as loc
import pandas as pd
import xesmf as xe
from datetime import datetime
import os

verbose = False
VAR = 'zos'
EXP = 'historical' # historical, rcp45, rcp85

year_min_ref = 1986  # Included. Beginning of reference period
year_max_ref = 2006  # Excluded. End of reference period

if EXP == 'historical':
    year_min = 1900 # Could start from 1850
    year_max = 2006
else:
    year_min = 2006
    year_max = 2101 
    
dir_outputs = '../outputs/'
dir_inputs = '../inputs/'

col_names = ['Centers','Models']
if EXP == 'historical':
    EXPm = 'rcp85'
else:
    EXPm = EXP
ModelList = pd.read_csv(dir_inputs+'CMIP5modelSelection_'+EXPm+'_'+VAR+'.txt', 
                        delim_whitespace=True, names=['Centers','Models'], 
                        comment='#')
Models = ModelList.Models

###### Start and end of each period 
years_s = np.arange(year_min,year_max)

# Star's project, historical period
#years_s = ispan(1900,2005,1)

#Read the regular 1*1 grid to use for regridded outputs
mask_ds = xr.open_dataset(dir_inputs+'reference_masks.nc')

dimLonOut   = len(mask_ds.lon)
dimLatOut   = len(mask_ds.lat)
weights = np.cos(np.deg2rad(mask_ds.lat))
weights.name = 'weights'

# Make a dataset to regrid with xESMF
ds_out = xr.Dataset({'lat': (['lat'], mask_ds.lat),
                     'lon': (['lon'], mask_ds.lon)})

print("Models used:")
print(Models)

for i in range(len(Models)):
    print(f'####### Working on model {i}, {Models[i]}  ######################')
    files_hist = loc.select_cmip5_files(VAR, 'historical', ModelList.Centers[i], 
                                        Models[i])
    if EXP != 'historical':
        files_sce = loc.select_cmip5_files(VAR, EXP, ModelList.Centers[i], 
                                           Models[i])
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
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'Could not open data from {Models[i]}')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        continue
    
    if EXP != 'historical':
        ds = xr.concat([hist_ds,sce_ds],'time')
    else
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
            print(f'Regridding did not work for {Models[i]}')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            continue
    
    if verbose:
        print(regridder)
    
    RefVAR1 = y_ds[VAR].sel(year=slice(year_min_ref,year_max_ref)).mean(dim='year')

    # Use the reference field to mask the small seas that are not connected to the
    # ocean and areas where sea ice is included on the ocean load
    RefVAR1_corr = RefVAR1 - RefVAR1.mean() #TODO weigthed mean
    MaskRefVAR1  = np.where((RefVAR1_corr>=2) | (RefVAR1_corr<=-2),np.nan,1)

    if verbose:
        print('Check info about piControl branching:')
        try:
            print(f"parent_experiment_id : {hist_ds.attrs['parent_experiment_id']}")
        except:
            print('No parent_experiment_id attribute')
        try:
            print(f"parent_experiment_rip : {hist_ds.attrs['parent_experiment_rip']}")
        except:
            print('No parent_experiment_rip attribute')
        try:
            print(f"branch_time : {hist_ds.attrs['branch_time']}")
        except:
            print('No branch_time attribute')
    try:
        conv_pic_hist = float(y_ds.year[0]) - float(hist_ds.attrs['branch_time'])
    except:
        # Pick a random large value that makes sure branching is not used in
        # trend_zos_pic_cmip5
        conv_pic_hist = -9999
    Trend_pic_coeff = loc.trend_zos_pic_cmip5(ModelList.iloc[i], order=1, 
                                              year_min=1850, year_max=2100,
                                              conv_pic_hist=conv_pic_hist)
    # Build polynmial from coefficients and convert from m to cm per year
    Trend_pic = xr.polyval(coord=y_ds.year, coeffs=Trend_pic_coeff)*100
    # Remove the average over the reference period
    Trend_pic = Trend_pic - Trend_pic.sel(year=slice(year_min_ref,year_max_ref)
                                         ).mean(dim='year')
    Trend_pic = Trend_pic.rename({Trend_pic.dims[1]:name_lat, 
                              Trend_pic.dims[2]:name_lon})
    
    MAT_CorrectedZOS_reg = np.zeros([len(years_s),dimLatOut,dimLonOut])
    ##### Loop on the years ######################################
    for idx, year in enumerate(years_s):
        print(f'Working on year: {year}')

        VAR1 = y_ds[VAR].sel(year=year)

        if (Models[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                          'MRI-CGCM3']): 
            VAR1 = np.where(VAR1==0,np.nan,VAR1)

        AnomVAR1 = (VAR1 - RefVAR1)*100 # Convert m to cm
        AnomVAR1 = AnomVAR1*MaskRefVAR1

        # Effective number of years to detrend: year of interest 
        # minus mean of reference period
        nbyears = year+0.5 - (year_max_ref+year_min_ref)/2

        TrendVAR1 = Trend_pic.sel(year=year)

        DTrendVAR1 = AnomVAR1 - TrendVAR1

        # Regrid to the reference 1*1 degree grid            
        DTrendVAR1_reg = regridder(DTrendVAR1)
    
        # Mask other problematic regions here
        if Models[i] in ['MIROC5', 'GFDL-ESM2M', 'GFDL-CM3','GISS-E2-R', 
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
                                        coords=[years_s+0.5, mask_ds.lat, mask_ds.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Models[i]]},0)
    MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})

    MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                       'ComputeOceanDynmicSeaLevel_CMIP5.py')
    MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')
    MAT_OUT_ds.attrs['emission_scenario'] = EXP

    NameOutput = f'{dir_outputs}CMIP5_{VAR}_{EXP}_{Models[i]}_{year_min}_{year_max}.nc'
    if os.path.isfile(NameOutput):
        os.remove(NameOutput)
    MAT_OUT_ds.to_netcdf(NameOutput)
