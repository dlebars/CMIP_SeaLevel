###############################################################################
# SelectPaths.py: Interaogate the local CMIP6 database and selects the paths
# to relevent data.
# To Do: Solve the di issue
###############################################################################

import sys
import os
import pandas as pd
import itertools
from pathlib import Path

# Options: 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'
exp_id = ['historical', 'piControl', 'ssp119']
var = ['zostoga'] # zostoga, zos

CMIP6_path = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'

# The depth of paths are sometimes used later, these refer to the CMIP6 folder
st = CMIP6_path.split('/')
st = list(filter(None, st)) # Remove empty '' strings
di = len(st)

### Function definitions ######################################################

def all_paths(data_dir)
    '''List all paths from a given directory'''
    list_paths = []
    for root, dirs, files in os.walk(data_dir):
        list_paths.append(root)
    return list_paths
    
def select_paths(data_dir, experiment_id, variable, ens1=False):
    '''Select all path with data for a given experiment_id and variable. 
    Outout results in a list'''
    list_paths = []
    for root, dirs, files in os.walk(data_dir):
        if files:
            st = root.split('/')
            st = list(filter(None, st)) # Remove empty '' strings
            if ens1:
                if (st[di+3] == experiment_id) and (st[di+6] == variable) \
                and (st[di+4] == 'r1i1p1f1'):
                    list_paths.append(root)
            else:
                if (st[di+3] == experiment_id) and (st[di+6] == variable):
                    list_paths.append(root)
    return list_paths

def ind_mod(list_paths):
    '''Takes a list of paths as input and provides a set of all the individual
    models available'''
    list_mod = []
    for path in list_paths:
        st = path.split('/')
        st = list(filter(None, st)) # Remove empty '' strings
        list_mod.append(st[di+2])
    return set(list_mod)

def select_models_intersection(experiment_id, var):
    '''Select the models that have available data for the input variable(s) and
    experiments'''
    list_com = list(itertools.product(experiment_id, var))
    for idx, val in enumerate(list_com):
        paths = select_paths(val[0], val[1], True)
        ind_mods = ind_mod(paths)
        if idx > 0:
            ind_mods = ind_mods.intersection(ind_mods_prev)
        ind_mods_prev = ind_mods
    return ind_mods

def make_info_df(list_path):
    '''Read list of path to data and store info into pandas dataframe '''
    interm_df = pd.DataFrame(columns=['Center', 'Model', 'Ensemble', 'Grid', 'Version'])
    for i in range(len(list_path)):
        st = list_path[i].split('/')
        st = list(filter(None, st)) # Remove empty '' strings
        interm_df.iloc[i] = [ st[di+1], st[di+2], st[di+4], st[di+7], st[di+8]]
    return interm_df

def make_final_info_df(info_df, ind_mods):
    '''Reads dataframe output from function make_info_df and select the final
    model grid and version to use'''
    nb_mod = len(ind_mods)
    final_df = pd.DataFrame(columns=['Center', 'Model', 'Ensemble', 'Grid', 'Version'])
    for i in range(nb_mod):
        info_sel_df = info_df[info_df['Model'] == ind_mod[i]]
        if len(set(info_sel_df['Ensemble'])) != 1:
            sys.exit('ERROR: More/less than 1 ensemble name for '+ind_mod[i])
        if len(set(info_sel_df['Grid'])) > 1:
            print('More than one grid available for '+ind_mod[i])
            print(set(info_sel_df['Grid']))
            print('Using gr as default')
            Grid = 'gr'
        else:
            Grid = info_sel_df['Grid'][0]
        All_Versions = set(info_sel_df['Version'])
        if len(All_Versions) > 1:
            print('More than one version available for '+ind_mod[i])
            print(All_Versions)
            print('Using lastest version:')
            Version = sorted(All_Versions)[-1]
            print(Version)
        else:
            Version = All_Versions
        final_df.loc[i] = [info_sel_df['Center'].iloc[0], ind_mod[i], 
                           info_sel_df['Ensemble'].iloc[0], Grid, Version]
    return final_df

# def make_full_paths(data_dir, df):
#     '''Read a pandas data frame and build the total path to data folders'''
#     dimMod = len(df.Model)
#     full_paths = []
#     for i in range(dimMod):
#         p = Path(data_dir+'CMIP/'+df.Center[i]+'/'+df.Model[i]+'/'+EXP[j]+
#                  '/'+Freq)
#         files1 = list(p.glob('*/*/*/*/'+VAR+'/*'+VAR+'*.nc'))
        
#         full_paths.append(?)
        
#     return full_paths
    
###############################################################################

ind_mods = select_models_intersection(exp_id, var)
list_all_paths = all_paths(CMIP6_path)
info_df = make_info_df(list_all_paths)
final_info_df = make_final_info_df(info_df, ind_mods)

final_info_df.to_csv(index=False)