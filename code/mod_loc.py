# mod_loc.py
from pathlib import Path

def select_cmip5_files(VAR, EXP, Center, Model):
    '''Return a list of paths to the CMIP5 data files'''
    DataDir  = '/nobackup/users/bars/synda/cmip5/output1/'
    p = Path(DataDir+Center+'/'+Model+
                 '/'+EXP+'/'+'mon')
    files = list(p.glob('*/*/*/*/'+VAR+'/*'+VAR+'*.nc'))
    # Select the last version of data: 
    vs = []
    for k in range(len(files)):
        part = files[k].parts
        vs.append(part[len(part)-3])
    vs.sort()
    files = sorted(p.glob('*/*/*/'+vs[-1]+'/'+VAR+'/*'+VAR+'*.nc')) # Switched list to sorted
    return files