###############################################################################
# ComputeOceanDynmicSeaLevel_CMIP5.py: 
# - Read zos variable from local CMIP5 files synchronized from ESGF nodes using 
# synda.
# - Correct the data using a common reference period and removing the trend 
# from PIcontrol simulations
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
EXP = 'rcp85'

year_min_ref = 1986  # Included. Beginning of reference period
year_max_ref = 2006  # Excluded. End of reference period
year_min = 2098 #2006
year_max = 2100 #2101

Dir_outputs = '../outputs/'
Dir_CMIP5_TE = '../../CMIP5_ThermalExp/'

col_names = ['Centers','Models']
ModelList = pd.read_csv(Dir_CMIP5_TE+'CMIP5modelSelection_'+EXP+'_'+VAR+'.txt', 
                        delim_whitespace=True, names=['Centers','Models'], 
                        comment='#')
Models = ModelList.Models

###### Start and end of each period 
years_s = np.arange(year_min,year_max)

# Star's project, historical period
#years_s = ispan(1900,2005,1)

#Read the regular 1*1 grid to use for regridded outputs
DIRgrid = '/nobackup/users/bars/SeaLevelFromHylke/CMIP5_OCEAN/Fingerprints/' # TODO: Add this standard grid somewhere to the project
rg = xr.open_dataset(DIRgrid+'Relative_icesheets.nc')
rg = rg.rename({'latitude':'lat', 'longitude':'lon'})

dimLonOut   = len(rg.lon)
dimLatOut   = len(rg.lat)
weights = np.cos(np.deg2rad(rg.lat))
weights.name = 'weights'

#Build a mask for the new grid
MaskOut = rg.DYN_ANT
MaskOut = MaskOut.where(MaskOut == 0, 1)
MaskOut = MaskOut.where(MaskOut != 0)
# Mask the Caspian sea
MaskOut.loc[dict(lat=slice(35,51), lon=slice(45,56))] = np.nan

# Mask that includes the Mediterranean and Black Sea
MaskOut_Med = MaskOut.copy()
MaskOut_Med.loc[dict(lat=slice(20.5,41.5), lon=slice(354,360))] = np.nan
MaskOut_Med.loc[dict(lat=slice(20.5,41.5), lon=slice(0,44))] = np.nan
MaskOut_Med.loc[dict(lat=slice(41,47.5), lon=slice(2,44))] = np.nan

# Mask that includes the Black Sea
MaskOut_BS = MaskOut.copy()
MaskOut_BS.loc[dict(lat=slice(39,46), lon=slice(26,42))] 

# Make a dataset for easy regridding with xESMF
ds_out = xr.Dataset({'lat': (['lat'], rg.lat),
                     'lon': (['lon'], rg.lon)})

print("Models used:")
print(Models)

ftrend = xr.open_dataset(f'{Dir_CMIP5_TE}TrendZOS_ForEXP{EXP}.nc')

for i in range(len(Models)):
    print(f'####### Working on model {i}, {Models[i]}  ######################')
    files_hist = loc.select_cmip5_files(VAR, 'historical', ModelList.Centers[i], 
                                Models[i])
    files_sce = loc.select_cmip5_files(VAR, EXP, ModelList.Centers[i], 
                                Models[i])
    if verbose:
        print('#### Using the following historical files: ####')
        print(files_hist)
        print('#### Using the following scenario files: ####')
        print(files_sce)

    try:
        hist_ds = xr.open_mfdataset(files_hist,combine='by_coords')
        sce_ds = xr.open_mfdataset(files_sce,combine='by_coords')
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'Could not open data from {Models[i]}')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        continue
    
    ds = xr.concat([hist_ds,sce_ds],'time')
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

    Trend_pic_coeff = loc.trend_zos_pic_cmip5(ModelList.iloc[i], order=1)
    # Build polynmial from coefficients and convert from m to cm per year
    Trend_pic = xr.polyval(coord=y_ds.year, coeffs=Trend_pic_coeff)*100
    Trend_pic = Trend_pic - Trend_pic.sel(year=slice(year_min_ref,year_max_ref)
                                         ).mean(dim='year')
    
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

#         TrendVAR1 = ftrend[Models[i]]*nbyears*100 # Call piCOntrol trend here
        TrendVAR1 = Trend_pic.sel(year=year)

        TrendVAR1 = TrendVAR1.rename({TrendVAR1.dims[0]:name_lat, 
                                      TrendVAR1.dims[1]:name_lon})

        DTrendVAR1 = AnomVAR1 - TrendVAR1

        # Regrid to the reference 1*1 degree grid            
        DTrendVAR1_reg = regridder(DTrendVAR1)
    
        # Mask other problematic regions here
        if Models[i] in ['MIROC5', 'GFDL-ESM2M', 'GFDL-CM3','GISS-E2-R', 
                         'GISS-E2-R-CC']:
            DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut_Med
        elif Models[i] in ['NorESM1-M', 'NorESM1-ME']:
            DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut_BS
        else:
            DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut           

        # Fill the NaN values
        DTrendVAR1_reg = DTrendVAR1_reg.interpolate_na('lon')

        DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut

        area_mean       = DTrendVAR1_reg.weighted(weights).mean(('lon', 'lat'))
        print(f'Removing area mean of:{np.round(area_mean.values,2)} cm')
        MAT_CorrectedZOS_reg[idx,:,:] = DTrendVAR1_reg - area_mean

    regridder.clean_weight_file()

    ### Export to a NetCDF file
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg, 
                                        coords=[years_s+0.5, rg.lat, rg.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Models[i]]},0)
    MAT_CorrectedZOS_reg.attrs['units'] = 'cm'
    MAT_CorrectedZOS_reg.attrs['regridding_method'] = f'xESMF package with {reg_method}'
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}': MAT_CorrectedZOS_reg})

    MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                       'ComputeOceanDynmicSeaLevel_CMIP5.py')
    MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')
    MAT_OUT_ds.attrs['emission_scenario'] = EXP

    NameOutput = f'{Dir_outputs}TEST2_CMIP5_{VAR}_{EXP}_{Models[i]}_{year_min}_{year_max}.nc'
    if os.path.isfile(NameOutput):
        os.remove(NameOutput)
    MAT_OUT_ds.to_netcdf(NameOutput)
