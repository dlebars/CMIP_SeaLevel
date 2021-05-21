###############################################################################
# SelectPaths.py: Interogate the local CMIP6 database and select the paths
# to relevent data.
###############################################################################

import sys
import os
import pandas as pd
import itertools
from pathlib import Path

### Function definitions ######################################################

def depth_path(data_dir):
    st = data_dir.split('/')
    st = list(filter(None, st)) # Remove empty '' strings
    return len(st)

def select_paths(data_dir, experiment_id, variable, ens1=False):
    '''Select all path with data for a given experiment_id and variable. 
    Outout results in a list'''
    depth = depth_path(data_dir)
    list_paths = []
    for root, dirs, files in os.walk(data_dir):
        if files:
            st = root.split('/')
            st = list(filter(None, st)) # Remove empty '' strings
            if ens1:
                if (st[depth+3] == experiment_id) and (st[depth+6] == variable) \
                and (st[depth+4] == 'r1i1p1f1'):
                    list_paths.append(root)
            else:
                if (st[depth+3] == experiment_id) and (st[depth+6] == variable):
                    list_paths.append(root)
    return list_paths

def select_ind_mod(list_paths, depth):
    '''Takes a list of paths as input and provides a set of all the individual
    models available'''
    list_mod = []
    for path in list_paths:
        st = path.split('/')
        st = list(filter(None, st)) # Remove empty '' strings
        list_mod.append(st[depth+2])
    return set(list_mod)

def select_models_intersection(data_dir, experiment_id, var, ens1):
    '''Select the models that have available data for the input variable(s) and
    experiments'''
    depth = depth_path(data_dir)
    # Make sure inputs are lists because itertools.product needs list inputs
    if not isinstance(var, list):
        var = [var]
    if not isinstance(experiment_id, list):
        experiment_id = [experiment_id]
    list_com = list(itertools.product(experiment_id, var))
    for idx, val in enumerate(list_com):
        paths = select_paths(data_dir, val[0], val[1], ens1)
        ind_mods = select_ind_mod(paths, depth)
        if idx > 0:
            ind_mods = ind_mods.intersection(ind_mods_prev)
        ind_mods_prev = ind_mods
    return list(ind_mods)

def make_info_df(list_path, depth):
    '''Read list of path to data and store info into pandas dataframe '''
    interm_df = pd.DataFrame(columns=['Center', 'Model', 'Ensemble', 'Grid', 'Version'])
    for i in range(len(list_path)):
        st = list_path[i].split('/')
        st = list(filter(None, st)) # Remove empty '' strings
        interm_df.loc[i] = [ st[depth+1], st[depth+2], st[depth+4], st[depth+7], st[depth+8]]
    return interm_df

def make_final_info_df(info_df, ind_mods):
    '''Reads dataframe output from function make_info_df and select the final
    model grid and version to use'''
    nb_mod = len(ind_mods)
    final_df = pd.DataFrame(columns=['Center', 'Model', 'Ensemble', 'Grid', 'Version'])
    for i in range(nb_mod):
        info_sel_df = info_df[info_df['Model'] == ind_mods[i]]
        print(info_sel_df)
        if len(set(info_sel_df['Ensemble'])) > 1:
            print('More than one ensemble available for '+ind_mods[i])
        if 'r1i1p1f1' in info_sel_df['Ensemble'].values:
            Ensemble = 'r1i1p1f1'
        elif 'r1i1p1f2' in info_sel_df['Ensemble'].values:
            Ensemble = 'r1i1p1f2'
        elif 'r1i1p1f3' in info_sel_df['Ensemble'].values:
            Ensemble = 'r1i1p1f3'
        elif 'r1i1p2f1' in info_sel_df['Ensemble'].values:
            Ensemble = 'r1i1p2f1'
        elif 'r4i1p1f1' in info_sel_df['Ensemble'].values:
            Ensemble = 'r4i1p1f1'
        else:
            print(set(info_sel_df['Ensemble']))
            sys.exit('ERROR: Standard ensemble not available see list above')
        print(f'Using ensemble {Ensemble}')
        info_sel_df = info_sel_df[info_sel_df['Ensemble'] == Ensemble]

        if len(set(info_sel_df['Grid'])) > 1:
            print('More than one grid available for '+ind_mods[i])
            print(set(info_sel_df['Grid']))
            print('Using gr as default or gr1 as default')
            if 'gr' in info_sel_df['Grid'].values:
                Grid = 'gr'
            elif 'gr1' in info_sel_df['Grid'].values:
                Grid = 'gr1'
            else:
                sys.exit('ERROR: Standard grid not available see list above')
        else:
            Grid = info_sel_df['Grid'].iloc[0]
        info_sel_df = info_sel_df[info_sel_df['Grid'] == Grid]
            
        All_Versions = set(info_sel_df['Version'])
        if len(All_Versions) > 1:
            print('More than one version available for '+ind_mods[i])
            print(All_Versions)
            print('Using lastest version:')
            Version = sorted(All_Versions)[-1]
            print(Version)
        else:
            Version = list(All_Versions)[0]
        final_df.loc[i] = [info_sel_df['Center'].iloc[0], ind_mods[i], 
                           Ensemble, Grid, Version]
    return final_df
    
###############################################################################

CMIP6_path = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'
depth = depth_path(CMIP6_path)
ens1 = True # True to select only r1i1p1f1, False otherwise

for var in ['zostoga', 'zos']:
    var = [var]
    for sce in ['historical','ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        print('####### Working on '+str(var)+', '+str(sce)+'#################'+
             '###############################################################')
        exp_id = [sce, 'historical', 'piControl']
        ind_mods = select_models_intersection(CMIP6_path, exp_id, var, ens1)
        print('Models available for this combination:')
        print(ind_mods)
        
        for idx, ei in enumerate(exp_id):
            list_all_paths = select_paths(CMIP6_path, ei, var[0], ens1)
            #print('\n'.join(list_all_paths))
            info_df = make_info_df(list_all_paths, depth)
            #print('Info before final selection:')
            #print(info_df)
            if idx == 0:
                final_info_df = make_final_info_df(info_df, ind_mods)
                final_info_df = final_info_df.rename(columns={'Version':ei+'_Version'})
            else:
                v_info_df = make_final_info_df(info_df, ind_mods)
                final_info_df[ei+'_Version'] = v_info_df.Version
        
        final_info_df.sort_values(by='Model', inplace=True)
        final_info_df.reset_index(drop=True, inplace=True)
        print('Final info to be saved as csv file:')
        print(final_info_df)
        
        if (len(var) == 1) & (len(exp_id) == 3):
            if sce == 'historical':
                file_name = (f'AvailableExperiments_{var[0]}_{exp_id[0]}_'+
                             f'{exp_id[1]}.csv')
            else:
                file_name = (f'AvailableExperiments_{var[0]}_{exp_id[0]}_'+
                             f'{exp_id[1]}_{exp_id[2]}.csv')
            final_info_df.to_csv(file_name, index=False)
        else:
            print('ERROR: Output file name not compatible with length of var'+
                  'or exp_id')