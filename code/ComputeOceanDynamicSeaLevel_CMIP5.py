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

def read_ModelNames(ds):
    '''Read model names, this need a little function because characters are not
    read properly by xarray'''
    mod_names = []
    for i in range(len(ds.ModelNames)):
        mod_names.append(str(ds.ModelNames[i].values)[2:-1])
    return mod_names

VAR = 'zos'
EXP = 'rcp85'

year_min_ref = 1986  # Included. Beginning of reference period
year_max_ref = 2006  # Excluded. End of reference period
year_min = 2006
year_max = 2101

Freq     = 'mon'  # Frequency of time output: mon or fx (only to read inputs)
DataDir  = '/nobackup/users/bars/synda/cmip5/output1/'

Dir_outputs = '../outputs/'
Dir_CMIP5_TE = '../../CMIP5_ThermalExp/'

col_names = ['Centers','Models']
ModelList = pd.read_csv(Dir_CMIP5_TE+'CMIP5modelSelection_'+EXP+'_'+VAR+'.txt', 
                        delim_whitespace=True, names=['Centers','Models'], 
                        comment='#')
Models = ModelList.Models

###### Start and end of each period 
# Erwin's project:
years_s = np.arange(year_min,year_max)
years_e = years_s+1
# Star's project, historical period
#years_s = ispan(1900,2005,1)
#years_e = ispan(1900+1,2005+1,1)

mid = (years_e + years_s)/2

#Read the regular 1*1 grid to use for regridded outputs
DIRgrid = '/nobackup/users/bars/SeaLevelFromHylke/CMIP5_OCEAN/Fingerprints/' # TODO: Add this standard grid somewhere to the project
rg = xr.open_dataset(DIRgrid+'Relative_icesheets.nc')
rg = rg.rename({'latitude':'lat', 'longitude':'lon'})

dimLonOut   = len(rg.lon)
dimLatOut   = len(rg.lat)
cLatOut     = np.cos(np.deg2rad(rg.lat))
#Build a mask for the new grid
MaskOut = rg.DYN_ANT
MaskOut = MaskOut.where(MaskOut == 0, 1)
MaskOut = MaskOut.where(MaskOut != 0)
# Mask the Caspian sea
MaskOut.loc[dict(lat=slice(35,51), lon=slice(45,56))] = np.nan
# TODO: Should the center of Australia be also masked?
# New mask that includes the mediterranean
MaskOut_Med = MaskOut
MaskOut_Med.loc[dict(lat=slice(20.5,41.5), lon=slice(354,360))] = np.nan
MaskOut_Med.loc[dict(lat=slice(20.5,41.5), lon=slice(0,44))] = np.nan
MaskOut_Med.loc[dict(lat=slice(41,47.5), lon=slice(2,44))] = np.nan

# Make a dataset for easy regridding with xESMF
ds_out = xr.Dataset({'lat': (['lat'], rg.lat),
                     'lon': (['lon'], rg.lon)})

print("Models used:")
print(Models)

#TODO later
fref   = xr.open_dataset(f'{Dir_CMIP5_TE}ReferenceZOS_ForEXP{EXP}_'
                         f'{year_min_ref}_{year_max_ref}.nc')
ftrend = xr.open_dataset(f'{Dir_CMIP5_TE}TrendZOS_ForEXP{EXP}.nc')

#Read the average zos fields to discount from the model sea level
# ;fzos_avg = addfile("CMIP5_SeaLevel_"+EXP+"_zos_avg_1950-2100.nc","r") ;For Sanne's project
fzos_avg = xr.open_dataset(f'{Dir_CMIP5_TE}CMIP5_SeaLevel_{EXP}_zos_avg_{year_min_ref}-2100.nc')
# ;fzos_avg = addfile("CMIP5_SeaLevel_"+EXP+"_zos_avg_1900-2006.nc","r") ;For Star's historical files

zos_avg_ModelNames = read_ModelNames(fzos_avg)
zos_avg  = fzos_avg.AverageSeaLevel
time_zos_avg = fzos_avg.time
print('Check the time vector:')
print(time_zos_avg)
indzosref = np.where((time_zos_avg >= year_min_ref) & (time_zos_avg < year_max_ref))[0]
print(indzosref)
zos_avg_ref =  zos_avg[:,indzosref].mean(axis=1)

for i in range(len(ModelList.Models)):
    print(f'####### Working on model {i}, {Models[i]}  #####')
    #### Read scenario data
    files1 = loc.select_cmip5_files(VAR, EXP, ModelList.Centers[i], 
                                ModelList.Models[i])
    print('#### Using following files: ####')
    print(files1)
    try:
        f1 = xr.open_mfdataset(files1,combine='by_coords')
    except:
        try:
            print('Open by_coords did not work for:'+ ModelList.Centers[i]+
                  '/'+ModelList.Models[i]+'/'+EXP)
            print('Using nested option instead')
            f1 = xr.open_mfdataset(files1,combine='nested', concat_dim='time')
        except:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(f'Could not open data from {Models[i]}')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            continue
    VAR1 = f1[VAR].squeeze()

    print(VAR1)
    time1  = f1.time
    dimt1  = len(time1)
    try:
        timeUT = time1.dt.year
    except:
        # Mix of time format, dt needs to be applied element wise
        timeUT = np.array([time1[i].dt.year.values.item() for i in range(len(time1))])
        timeUT = xr.DataArray(timeUT, coords=[timeUT], dims=['time'])

    #### Read historical data    
    files2 = loc.select_cmip5_files(VAR, 'historical', ModelList.Centers[i], 
                                    ModelList.Models[i])
    
    print('### Also using these historical files: ###')
    print(files2)
    f2 = xr.open_mfdataset(files2,combine='by_coords')
    time2 = f2.time
    timeUT2 = time2.dt.year
    
    if len(f1.lat.shape) == 1:
        name_lat = 'lat'
        name_lon = 'lon'
    elif len(f1.lat.shape) == 2:
        name_lat = 'rlat'
        name_lon = 'rlon'        
    
    if Models[i] in ['IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','NorESM1-M', 
                     'NorESM1-ME','ACCESS1-0','CCSM4','MPI-ESM-LR','MPI-ESM-MR']:
        f1 = f1.rename({'j':name_lat, 'i':name_lon})
        f2 = f2.rename({'j':name_lat, 'i':name_lon})
    
    try:
        regridder = xe.Regridder(f1, ds_out, 'bilinear',periodic=True)
        print(regridder)
    except:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'Regridding did not work for {Models[i]}')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        continue
    
    # TODO: This lon might not be necessary, only used in interpollation
#     if (Models(i).eq."bcc-csm1-1").or.(Models(i).eq."bcc-csm1-1-m").or. \
#      (Models(i).eq."GFDL-ESM2G").or.(Models(i).eq."GFDL-ESM2M").or. \
#      (Models(i).eq."GFDL-CM3") then
#         lon = np.where(lon<0,lon+360,lon)

    RefVAR1 = fref[Models[i]]
    RefVAR1 = RefVAR1.rename({RefVAR1.dims[0]:name_lat,
                              RefVAR1.dims[1]:name_lon})
    # Use the reference field to mask the small seas that are not connected to the
    # ocean and areas where sea ice is included on the ocean load
    RefVAR1_corr = RefVAR1 - RefVAR1.mean()
    MaskRefVAR1  = np.where((RefVAR1_corr>=2) | (RefVAR1_corr<=-2),np.nan,1)

    MAT_CorrectedZOS_reg = np.zeros([len(years_s),dimLatOut,dimLonOut])
    
    ##### Loop on the years ######################################
    for y in [0,1]: #range(len(years_s)):
        print(f'Working on period: {years_s[y]}-{years_e[y]}')
        indzossel   = np.where((time_zos_avg>=years_s[y]) & (time_zos_avg<years_e[y]))[0]
        if indzossel.size == 0:
            print(f'No data during the period: {years_s[y]}-{years_e[y]}')
            MAT_CorrectedZOS_reg[y,:,:] = np.nan
        else:
            if len(indzossel) == 1:
                zos_avg_sel = zos_avg[:,indzossel]
            else:
                zos_avg_sel = zos_avg[:,indzossel].mean(axis=1)

            if years_s[y]>=2006:
                ind_time_sel = np.where((timeUT>=years_s[y]) & (timeUT<=years_e[y]))[0]
                VAR1 = f1[VAR][ind_time_sel,:,:]
            else:
                ind_time_sel = np.where((timeUT2>=years_s[y]) & (timeUT2<=years_e[y]))[0]
                VAR1 = f2[VAR][ind_time_sel,:,:]

            if (Models[i] in ['MIROC5', 'GISS-E2-R', 'GISS-E2-R-CC', 'EC-EARTH', 
                              'MRI-CGCM3']): 
                VAR1 = np.where(VAR1==0,np.nan,VAR1)

            VAR1avg = VAR1.mean(axis=0) # Compute time average
            ind_zos_avg = zos_avg_ModelNames.index(Models[i])
            AnomVAR1 = (VAR1avg - RefVAR1)*100 # Convert m to cm
            AnomVAR1 = AnomVAR1*MaskRefVAR1

            ZOS_AVG_CORR = (zos_avg_sel[ind_zos_avg] - 
                            zos_avg_ref[ind_zos_avg])*100

            # Effective number of years to detrend: year of interest 
            # minus mean of reference period
            nbyears = mid[y] - (year_max_ref+year_min_ref)/2

            TrendVAR1 = ftrend[Models[i]]*nbyears*100
            TrendVAR1 = TrendVAR1.rename({TrendVAR1.dims[0]:name_lat, 
                                          TrendVAR1.dims[1]:name_lon})
            DTrendVAR1 = AnomVAR1 - TrendVAR1 - ZOS_AVG_CORR           
            print(DTrendVAR1)
            # Regrid to the reference 1*1 degree grid            
            DTrendVAR1_reg = regridder(DTrendVAR1.isel(time=0)) 
            #!!! lat/lon should be right most dimensions

            # Mask other problematic regions here
#           DTrendVAR1_reg@_FillValue = 1e+20
#           if (Models(i).eq."MIROC5").or.(Models(i).eq."GFDL-ESM2M").or. \
#              (Models(i).eq."GISS-E2-R").or.(Models(i).eq."GISS-E2-R-CC") then
#             DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut_Med
#             else
#               DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut
#           end if
            
            # Fill the NaN values TODO              
            # NCL command: poisson_grid_fill(DTrendVAR1_reg,True,1,100,1,0.5,0)
            DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut
            

#             area_mean       = wgt_areaave_Wrap(DTrendVAR1_reg, cLatOut, 1.0, 0) TODO
#             print("Removing area mean of:" + area_mean + " cm")
            MAT_CorrectedZOS_reg[y,:,:] = DTrendVAR1_reg #- area_mean

    regridder.clean_weight_file()

    ### Export in NetCDF file
    MAT_CorrectedZOS_reg = xr.DataArray(MAT_CorrectedZOS_reg, 
                                        coords=[mid, rg.lat, rg.lon], 
                                        dims=['time', 'lat', 'lon'])
    MAT_CorrectedZOS_reg = MAT_CorrectedZOS_reg.expand_dims({'model': [Models[i]]},0)
    
    MAT_OUT_ds = xr.Dataset({f'CorrectedReggrided_{VAR}_{EXP}': MAT_CorrectedZOS_reg})

    MAT_OUT_ds.attrs['source_file'] = ('This NetCDF file was built from '+ 
                                       'ComputeOceanDynmicSeaLevel_CMIP5.py')
    MAT_OUT_ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

    NameOutput = f'{Dir_outputs}CMIP5_{VAR}_{EXP}_{Models[i]}_{year_min}_{year_max}.nc'
    if os.path.isfile(NameOutput):
        os.remove(NameOutput)
    MAT_OUT_ds.to_netcdf(NameOutput)
