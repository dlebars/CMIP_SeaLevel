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

def select_paths(data_dir, experiment_id, variable, var):
    '''Select all path with data for a given experiment_id and variable. 
    Outout results in a list'''
    depth = depth_path(data_dir)
    list_paths = []
    for root, dirs, files in os.walk(data_dir):
        if files:
            st = root.split('/')
            st = list(filter(None, st)) # Remove empty '' strings
            if var=='any':
                if (st[depth+3] == experiment_id) and (st[depth+6] == variable):
                    list_paths.append(root)
            else:
                if (st[depth+3] == experiment_id) and (st[depth+6] == variable) \
                and (st[depth+4] == var):
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

def select_models_intersection(data_dir, experiment_id, variable, var):
    '''Select the models that have available data for the input variable(s) and
    experiments'''
    depth = depth_path(data_dir)
    # Make sure inputs are lists because itertools.product needs list inputs
    if not isinstance(variable, list):
        variable = [variable]
    if not isinstance(experiment_id, list):
        experiment_id = [experiment_id]
    list_com = list(itertools.product(experiment_id, variable))
    for idx, val in enumerate(list_com):
        paths = select_paths(data_dir, val[0], val[1], var)
        ind_mods = select_ind_mod(paths, depth)
        if idx > 0:
            ind_mods = ind_mods.intersection(ind_mods_prev)
        ind_mods_prev = ind_mods
    return list(ind_mods)

def make_info_df(list_path, depth):
    '''Read list of path to data and store info into pandas dataframe '''
    interm_df = pd.DataFrame(columns=['Center', 'Model', 'Variant', 'Grid', 'Version'])
    for i in range(len(list_path)):
        st = list_path[i].split('/')
        st = list(filter(None, st)) # Remove empty '' strings
        interm_df.loc[i] = [ st[depth+1], st[depth+2], st[depth+4], st[depth+7], st[depth+8]]
    return interm_df

def make_final_info_df(info_df, ind_mods, info_exp, exp):
    '''Reads dataframe output from function make_info_df and select the final
    model grid and version to use'''
    nb_mod = len(ind_mods)
    final_df = pd.DataFrame(columns=['Center', 'Model', 'Variant', 'Grid', 'Version'])
    for i in range(nb_mod):
        info_sel_df = info_df[info_df['Model'] == ind_mods[i]]
        print(info_sel_df)
        
        if len(set(info_sel_df['Variant'])) > 1:
            print('More than one variant available for '+ind_mods[i])
            
        if info_exp is not None:
            info_exp_sel_df = info_exp[info_exp['Model'] == ind_mods[i]]
            Variant = info_exp_sel_df[f'{exp}_Variant'].values[0]
        else:
            if 'r1i1p1f1' in info_sel_df['Variant'].values:
                Variant = 'r1i1p1f1'
            elif 'r1i1p1f2' in info_sel_df['Variant'].values:
                Variant = 'r1i1p1f2'
            elif 'r1i1p1f3' in info_sel_df['Variant'].values:
                Variant = 'r1i1p1f3'
            elif 'r1i1p2f1' in info_sel_df['Variant'].values:
                Variant = 'r1i1p2f1'
            elif 'r4i1p1f1' in info_sel_df['Variant'].values:
                Variant = 'r4i1p1f1'
            else:
                print(set(info_sel_df['Variant']))
                sys.exit('ERROR: Standard variant not available see list above')
                
        print(f'Using variant {Variant}')
        info_sel_df = info_sel_df[info_sel_df['Variant'] == Variant]

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
                           Variant, Grid, Version]
    return final_df
    
###############################################################################

CMIP6_path = '/nobackup_1/users/bars/synda_cmip6/CMIP6/'
depth = depth_path(CMIP6_path)
# var can either the name of a specific variant like r1i1p1f1
# or 'any' to automatically pick any variant available
var = 'any'

for variable in ['zostoga', 'zos']:
    variable = [variable]
    for sce in ['historical','ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        print('####### Working on '+str(variable)+', '+str(sce)+'#################'+
             '###############################################################')
        exp_id = [sce, 'historical', 'piControl']
        ind_mods = select_models_intersection(CMIP6_path, exp_id, variable, var)
        print('Models available for this combination:')
        print(ind_mods)
        
        for idx, ei in enumerate(exp_id):
            print(f'####### Selecting {ei} files with idx={idx}######')
            list_all_paths = select_paths(CMIP6_path, ei, variable[0], var)
            #print('\n'.join(list_all_paths))
            info_df = make_info_df(list_all_paths, depth)
            #print('Info before final selection:')
            #print(info_df)
            if idx == 0:
                final_info_df = make_final_info_df(info_df, ind_mods, None, None)
                final_info_df = final_info_df.rename(
                    columns={'Version':ei+'_Version', 'Variant':ei+'_Variant'})
            else:
                if var=='any' and idx==1 and sce!='historical':
                    # Historical should use the same variant as scenario
                    print('!!!!!! Using new case !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    v_info_df = make_final_info_df(info_df, ind_mods, 
                                                   final_info_df, sce)
                else:
                    v_info_df = make_final_info_df(info_df, ind_mods, None, None)
                    
                final_info_df[ei+'_Variant'] = v_info_df.Variant
                final_info_df[ei+'_Version'] = v_info_df.Version
        
        final_info_df.sort_values(by='Model', inplace=True)
        final_info_df.reset_index(drop=True, inplace=True)
        print('Final info to be saved as csv file:')
        print(final_info_df)
        
        if (len(variable) == 1) & (len(exp_id) == 3):
            if sce == 'historical':
                file_name = (f'AvailableExperiments_{variable[0]}_{exp_id[0]}_'+
                             f'{exp_id[1]}.csv')
            else:
                file_name = (f'AvailableExperiments_{variable[0]}_{exp_id[0]}_'+
                             f'{exp_id[1]}_{exp_id[2]}.csv')
            final_info_df.to_csv(file_name, index=False)
        else:
            print('ERROR: Output file name not compatible with length of variable'+
                  'or exp_id')