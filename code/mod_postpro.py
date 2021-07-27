# mod_postpro.py
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

def define_area(reg):
    '''Provides box coordinates given a region name'''
    
    if reg == 'dutch_coast':
        lon_min, lon_max = 3, 7
        lat_min, lat_max = 51, 54
    elif reg == 'north_sea':
        lon_min, lon_max = -2, 9
        lat_min, lat_max = 48, 60
    elif reg == 'knmi14_reg':
        lon_min, lon_max = -3.5, 7.5
        lat_min, lat_max = 51, 60
    
    return lon_min, lon_max, lat_min, lat_max

def ds2df(cmip_ds, lon_min, lon_max, lat_min, lat_max, start_year, end_year):
    '''Transform a dataset to a dataframe by averaging sea level over a region'''
    
    sel_ds = cmip_ds.sel(lon=slice(lon_min,lon_max), lat=slice(lat_min,lat_max))
    sel_da = sel_ds['CorrectedReggrided_zos'].mean(dim=['lon','lat'])

    df = pd.DataFrame(dict(time=np.arange(start_year,end_year)+0.5))
    df = df.set_index('time')

    for mod in sel_da.model.values:
        df[mod] = sel_da.sel(model=mod).drop('model').to_dataframe()
        
    return df

def read_zos_ds(data_dir, mip, sce):
    '''Read both historical and scenario datasets, select the intersecting 
    models and concatenate the two datasets'''
    
    hist_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zos_historical/{mip}_zos_historical_*.nc')
    sce_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zos_{sce}/{mip}_zos_{sce}_*.nc')

    model_intersection = list(set(hist_ds.model.values) & 
                              set(sce_ds.model.values))
    model_intersection.sort()
    tot_ds = xr.concat([hist_ds,sce_ds],'time').sel(model=model_intersection)
    
    return tot_ds

def read_zostoga_ds(data_dir, mip, sce):
    '''Read both historical and scenario datasets, select the intersecting 
    models and concatenate the two datasets'''
    
    hist_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zostoga/{mip}_zostoga_historical_*.nc')
    sce_ds = xr.open_mfdataset(
        f'{data_dir}/{mip}_zostoga/{mip}_zostoga_{sce}_*.nc')

    model_intersection = list(set(hist_ds.model.values) & 
                              set(sce_ds.model.values))
    model_intersection.sort()
    tot_ds = xr.concat([hist_ds,sce_ds],'time').sel(model=model_intersection)
    
    return tot_ds

def plot_cmip5_cmip6_ensembles(cmip5_sce, cmip5_df, cmip6_sce, cmip6_df, 
                               lower_bound, upper_bound, ra):
    '''Takes two dataframes from model ensembles as inputs and make a figure of 
    median and some quantiles to represente the divergence between models of the 
    ensemble.'''
    
    sce = f'{cmip5_sce}_{cmip6_sce}'
    
    cmip5_med = cmip5_df.quantile(0.5, axis=1).rolling(ra, center=True).mean()
    cmip5_lb = cmip5_df.quantile(lower_bound, axis=1).rolling(ra, center=True).mean()
    cmip5_ub = cmip5_df.quantile(upper_bound, axis=1).rolling(ra, center=True).mean()

    cmip6_med = cmip6_df.quantile(0.5, axis=1).rolling(ra, center=True).mean()
    cmip6_lb = cmip6_df.quantile(lower_bound, axis=1).rolling(ra, center=True).mean()
    cmip6_ub = cmip6_df.quantile(upper_bound, axis=1).rolling(ra, center=True).mean()

    fig, ax = plt.subplots(figsize=(6,6))
    ax.fill_between(cmip5_ub.index, 
                    cmip5_ub, 
                    cmip5_lb, 
                    color='blue', alpha=0.5,
                    label=f'cmip5, {int(lower_bound*100)}-{int(upper_bound*100)} percentiles')

    ax.fill_between(cmip6_ub.index, 
                    cmip6_ub, 
                    cmip6_lb, 
                    color='orange', alpha=0.5,
                    label=f'cmip6, {int(lower_bound*100)}-{int(upper_bound*100)} percentiles')

    ax.plot(cmip5_med, color='blue', label='cmip5 median')
    ax.plot(cmip6_med, color='orange', label='cmip6 median')

    plt.xlabel('time')
    plt.ylabel('sea level (cm)')
    plt.title(f'Compare zos projections from cmip5 and cmip6 for {sce} \n'+
              f'with reference period 1986-2005 and running average of {ra} years')
    ax.grid(True)
    plt.legend(loc='upper left')
    
def summary_fig_and_table(ax, df, colors=None, vlines=False):
    
    mi = 0.6 # Max color intensity
    
    # Get some pastel shades for the colors
    if not(colors):
        colors = plt.cm.Oranges(np.linspace(0, mi, len(df.index)))
        rowColours = colors
        
        # Expand the array
        ones = np.ones(len(df.columns))
        colors = colors[np.newaxis,:,:] * ones[:, np.newaxis, np.newaxis]
        
    elif colors=='alternate':
        colors1 = plt.cm.Oranges(np.linspace(0, mi, len(df.index)))
        colors2 = plt.cm.Blues(np.linspace(0, mi, len(df.index)))
        colors = np.zeros([len(df.columns), len(df.index), 4])
        colors[::2] = colors1
        colors[1::2] = colors2
        
        rowColours = plt.cm.Greys(np.linspace(0, mi, len(df.index)))

    # Start from white color
    colors[:,0,:] = 0
    
    index = np.arange(len(df.columns))
    bar_width = 0.6

    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(df.columns))

    # Plot bars and create text labels for the table
    cell_text = []
    for row in range(len(df.index)):
        ax.bar(index, 
               df.iloc[row]-y_offset, 
               bar_width, 
               bottom=y_offset, 
               color=colors[:,row,:])
        
        y_offset = df.iloc[row]
        cell_text.append(['%1.1f' % x for x in df.iloc[row]])
    
    ax.set_xlim(-0.5,index[-1]+0.5)

    # Add a table at the bottom of the axes
    ax.table(cellText=cell_text[::-1],
             rowLabels=df.index[::-1],
             rowColours=rowColours[::-1],
             colColours=colors[:,2,:],
             colLabels=df.columns,
             loc='bottom')

    ax.set_xticks([])
    
    if vlines:
        xcoords = index[:-1]+0.5
        xcoords = xcoords[::2]
        for xc in xcoords:
            plt.axvline(x=xc, color='black', linewidth=0.5, linestyle='--')
    
    return ax