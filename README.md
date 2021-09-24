# CMIP_SeaLevel

Prepare CMIP5 and CMIP6 data downloaded using synda to easily analyse it. It provides an easy way to clean CMIP data, compute yearly averages, enforce a spatial mean of zero (for zos variable), detrend the data based on piControl simulations (both glbal mean data like zostoga and plane data like zos) and reggrid different models on a common grid.

The code works only for data in the real of Amon and Omon.

### How to use this code

The first step is to download CMIP data. To lern more about CMIP, [this guide](https://pcmdi.llnl.gov/CMIP6/Guide/dataUsers.html) is a good place to start . Downloading is best done using the [Synda Tranfer Module](https://prodiguer.github.io/synda/sdt/sdt.html). In particular, use the [ESGF replication](https://prodiguer.github.io/synda/sdt/replication.html) that recreates the original ESGF path tree locally. 

This is a complex structure so first run SelectPaths.py in /SelectPaths_CMIP6 to create .csv files listing the datasets that were downloaded locally by Synda.

The code reads these .csv files to search through the local CMIP data.

The notebooks provide example of how the resulting datasets can be analysed using xarray.

### Dependencies:
- Python 3
- numpy
- xarray
- xESMF from [Pangeo-data](https://github.com/pangeo-data/xESMF) for the reggriding


