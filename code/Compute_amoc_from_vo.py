###############################################################################
# Compute_amoc_from_vo.py:
# Compute time series of AMOC at different latitudes using the meridional 
# velocity fields from CMIP5 amd CMIP6
#
# A few assumptions are made about the grid:
# 1- Grid cells are equally spaces in the zonal direction at the LAT_SEL latitudes
# 2- Zonal direction is approximately in line with y indexes even for 2D 
# curvilinear grids
###############################################################################

import itertools
import os

import numpy as np
import xarray as xr

import mod_loc as loc

verbose = True # Print additional information
LAT_SEL = [26, 35] # Latitude to select
MIP = 'cmip6'

lon_min, lon_max = 280, 355
lat_min, lat_max = 20, 40 # Only used to crop data

# Compute the width of the Atlantic along selected latitudes in meters
RE = 6.371e6 # Radius of the Earth
LAT_SEL_DIST = np.cos(np.radians(LAT_SEL))*RE*(lon_max-lon_min)

if MIP == 'cmip5':
    VAR = 'vo'
    EXP = ['historical', 'rcp26', 'rcp45','rcp85']
    
elif MIP == 'cmip6':
    VAR = 'vo' # msftmz, msftyz
    EXP = ['historical'] #, 'ssp126', 'ssp245', 'ssp585']

dir_inputs = '../inputs/'
dir_outputs = f'/nobackup/users/bars/{MIP.upper()}_regridded/'

def open_files(file_names):
    try:
        ds = xr.open_mfdataset(file_names, combine='by_coords', use_cftime=True)
    except:
        print(f'Open by_coords did not work for:{file_names}')
        print('Using nested option instead')
        try:
            ds = xr.open_mfdataset(file_names, combine='nested', 
                                   concat_dim='time', 
                                   use_cftime=True)
        except:
            print('Nested option did not work either. Decoding time manually.')
            
            try:
                ds = xr.open_mfdataset(file_names, combine='by_coords', decode_times = False)
                ds.time.attrs['units'] = 'days since 0000-01-01'
                ds['time'].data[:] = ds.time.data[:] + 365
                ds = xr.decode_cf(ds)
                
            except:
                print('ERROR: Reading dataset was not successful')

    return ds

def harmonize_lat_lon_lev_names(ds):
    '''Harmonise the name of longitude and latitude variables in a CMIP dataset'''
    
    if ('lat' not in ds.coords):
        if 'latitude' in ds.coords:
            ds = ds.rename({'latitude':'lat'})
        elif 'rlat' in ds.coords:
            ds = ds.rename({'rlat':'lat'})
        elif 'nav_lat' in ds.coords:
            ds = ds.rename({'nav_lat':'lat'})

    if ('lon' not in ds.coords):
        if 'longitude' in ds.coords:
            ds = ds.rename({'longitude':'lon'})
        elif 'rlon' in ds.coords:
            ds = ds.rename({'rlon':'lon'})
        elif 'nav_lon' in ds.coords:
            ds = ds.rename({'nav_lon':'lon'})
            
    if 'olevel' in ds.coords:
        ds = ds.rename({'olevel':'lev'})
            
    return ds

def change_longitude_values(ds): 
        
    print('Changing longitude from -180,0 to 180,360')
    
    new_lon = np.where(ds['lon'].values<0, ds['lon'].values+360, ds['lon'].values)
    ds['lon'].values = new_lon
    
    return ds

for exp in EXP:
    print(f'Working on VAR {VAR}, exp {exp}')

    year_min, year_max, ref_p_min, ref_p_max = loc.start_end_ref_dates(MIP, exp)
    print(f'Generating a file for this period: {year_min}-{year_max-1}, including {year_max-1}')

    ModelList = loc.read_model_list(dir_inputs, MIP, exp, VAR, False)

    print(ModelList)

    dimMod = len(ModelList.Model)
    time_all = np.arange(year_min, year_max ) + 0.5
    dimt = len(time_all)
    print(f'Number of years: {dimt}')

    da = xr.DataArray(np.zeros([dimMod, dimt, len(LAT_SEL)]), 
                      coords=[ModelList.Model, time_all, LAT_SEL], 
                      dims=['model', 'time', 'latitude'])

    for i in range(dimMod):
        print('###### Working on new model ##################################')
        
        files = loc.select_files(MIP, exp, VAR, ModelList.loc[i], verbose)
        ds = open_files(files)
        
        print(ds)
        
        loc_da = ds[VAR].squeeze()
        #y_da = loc.yearly_mean(loc_da)
        
        loc_da = harmonize_lat_lon_lev_names(loc_da)
        ds = harmonize_lat_lon_lev_names(ds)
        
        if np.any(loc_da['lon']<0):
            loc_da = change_longitude_values(loc_da)
        
        if (len(loc_da['lon'].shape)) == 2:
            print('### Two dimensional lon/lat arrays ###')
        
            mask_lon = (loc_da['lon'] > lon_min) & (loc_da['lon'] < lon_max)
            mask_lat = (loc_da['lat'] >= lat_min) & (loc_da['lat'] <= lat_max)
        
            cropped_da = loc_da.where(mask_lon & mask_lat, drop=True)
            
            #yearly avg
            
        elif (len(loc_da['lon'].shape)) == 1:
            print('### One dimension lon/lat arrays ###')
            
            cropped_da = loc_da.sel(lon=slice(lon_min, lon_max))
            
            for indl, lat_sel in enumerate(LAT_SEL):
                lat_c = cropped_da.sel(lat=lat_sel, method='nearest')
                
                lat_c = loc.yearly_mean(lat_c)
                
                dlon = (lat_c['lon'][1:].values-lat_c['lon'][:-1].values)
                dlon = np.append(dlon, dlon[-1])
                zonal_length = np.cos(np.radians(lat_sel))*RE*dlon
                
                # Remove time dependence of depth
                if 'time' in ds['lev_bnds'].coords:
                    lev_bnds = ds['lev_bnds'].isel(time=0)
                else:
                    lev_bnds = ds['lev_bnds']
                    
                vertical_length = lev_bnds.isel(bnds=1)-lev_bnds.isel(bnds=0)
                    
                vo_vol = lat_c*zonal_length*vertical_length
                
                zonal_sum = vo_vol.sum(dim='lon')
                vertical_cumsum = zonal_sum.cumsum(dim='lev')
                
                amoc = vertical_cumsum.max(dim='lev')
                
                print(amoc)
                
                try:
                    da[i,:,indl] = amoc.sel(time=slice(year_min,year_max))
                except:
                    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    print('There seem to be missing time data, model not used')
    
    da = da/10e6 # Convert m3/s to Sv
    
    da.attrs['long_name'] = 'AMOC volume transport'
    da.attrs['description'] = ('AMOC volume transport computed fom integrating the'+ 
                               'meridional velocity')
    da.attrs['units'] = 'Sv'
    ds = xr.Dataset()
    ds['amoc'] = da

    ### Remove missing models
    ds = ds.where(ds[VAR]!=0)
    ds = ds.dropna('model',how='all')

    ds = ds.sel(time=slice(year_min,year_max))

    print("### Export data to a NetCDF file ######################################")
    script_name = os.path.basename(__file__)
    name_output = f'{dir_outputs}{MIP}_amoc_{VAR}_{exp}_{year_min}_{year_max-1}.nc'
    loc.export2netcdf(ds, name_output, script_name)