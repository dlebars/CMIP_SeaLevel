# ExtractReferenceGrid.py: 
# - Extract reference grid from AR5 
# - Make general land masks from it
# - Export in a new NetCDF file for further use

import xarray as xr
import numpy as np
from datetime import datetime
import os

# Read the regular 1*1 grid to use for regridded outputs
dir_grid = '/Users/dewilebars/Projects/Project_ProbSLR/Data_Proj/Data_AR5/Fingerprints/'
dir_inputs = '../inputs/'
rg = xr.open_dataset(dir_grid+'Relative_icesheets.nc')
rg = rg.rename({'latitude':'lat', 'longitude':'lon'})

#Build a mask for the new grid
mask = rg.DYN_ANT
mask = mask.where(mask == 0, 1)
mask = mask.where(mask != 0)

# Mask the Caspian sea
mask.loc[dict(lat=slice(35,51), lon=slice(45,56))] = np.nan

# Mask the Black Sea
mask.loc[dict(lat=slice(39,48), lon=slice(26,42))] = np.nan

mask.attrs['long_name'] = ('Reference land mask that includes the Caspian Sea '+
                           'and the Black Sea')

# Mask that includes the Mediterranean sea as well
mask_med = mask.copy()
mask_med.loc[dict(lat=slice(20.5,41.5), lon=slice(354,360))] = np.nan
mask_med.loc[dict(lat=slice(20.5,41.5), lon=slice(0,44))] = np.nan
mask_med.loc[dict(lat=slice(41,47.5), lon=slice(2,44))] = np.nan

mask_med.attrs['long_name'] = ('Land mask that includes the '+
                           'Mediterranean Sea as well')

ds = xr.Dataset({'mask': mask, 'mask_med': mask_med})
ds.attrs['source_file'] = ('This NetCDF file was built from the script '+
                           'ExtractReferenceGrid.py')
ds.attrs['creation_date'] = datetime.now().strftime('%Y-%m-%d %H:%M')

NameOutput = dir_inputs+'reference_masks.nc'
if os.path.isfile(NameOutput):
    os.remove(NameOutput)
ds.to_netcdf(NameOutput)