# mod_loc.py
from pathlib import Path
import numpy as np
import xarray as xr

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
    files = sorted(p.glob('*/*/*/'+vs[-1]+'/'+VAR+'/*'+VAR+'*.nc'))
    return files

def yearly_mean(ds):
    '''Convert the data set or data array to year'''
    try:
        ds.coords['year'] = ds.time.dt.year
    except:
        years = np.array([ds.time[i].dt.year.values.item() for i in range(len(ds.time))])
        ds.coords['year'] = xr.DataArray(years, dims=['time'])
    y_ds   = ds.groupby('year').mean(dim='time')
    return y_ds
    

# def trend_zos_pic_cmip5(ForEXP, ModelList):
#     '''Compute zos trend over the pre-industrial control model simulations'''
#     # - Could it work for zos as well? It would make the ComputeGlobalThermalExpansion 
#     #scripts simpler
#     # - Possibility to compare linear and 2nd order detrending?
#     # Work also on different experiments? No.
#     # CMIP5 and CMIP6?
    
#     VAR = "zos" # To remove if the script doesn't work for zostoga
#     EXP = "piControl"

#     year_min = 1986
#     year_max = 2100

# Freq     = "mon"  ; Frequency of time output: mon or fx (only to read inputs)
# DataDir  = "/nobackup/users/bars/synda/cmip5/output1/"

# ModelList = readAsciiTable("CMIP5modelSelection_"+ForEXP+"_"+VAR+".txt",1,"string",1)

# delim    = " "
# Centers  = rm_single_dims(str_get_field(ModelList, 1, delim))
# Models   = rm_single_dims(str_get_field(ModelList, 2, delim))

# dimMod   = dimsizes(Models)

# do i=0,dimMod-1
#   print("####### Working on model "+i+","+Models(i)+"  #####")

#     files_hist = loc.select_cmip5_files(VAR, 'historical', ModelList.Centers[i], 
#                                 ModelList.Models[i])
    
#   file_name = DataDir+Centers(i)+"/"+Models(i)+"/"+EXP+"/"+Freq
#   files1 = systemfunc("ls "+file_name+"/*/*/*/*/"+VAR(0)+"/*"+VAR(0)+"*.nc")
#   ;Use last version of data:
#   vs1    = str_get_field(files1, 14,"/")
#   delete(files1)
#   files1 = systemfunc("ls "+file_name+"/*/*/*/"+vs1(dimsizes(vs1)-1)+"/"+VAR+"/*"+VAR+"*.nc")
#   print("#### Using following files: ####")
#   print(files1)
#   f1      = addfiles(files1,"r")
#   f1s     = addfile(files1(0),"r") ; Need this trick to avoid issues of adding
#                                     ; lon or lat coordinate arrays from addfiles
#   time1  = f1[:]->time
#   dimt1  = dimsizes(time1)

#   if time1@calendar.eq."proleptic_gregorian" then
#     time1@calendar = "gregorian"
#   end if
#   timeUT = cd_calendar(time1, 4)

#       ;Assumes piControl simulation starts in 1850 (some models start from 0)
#   timeUT = timeUT - timeUT(0) +1850.5

#   ind_time_sel = ind((timeUT.ge.year_min).and.(timeUT.le.year_max))
#   timeUTsel    = timeUT(ind_time_sel)

#   lat    = f1s->lat
#   lon    = f1s->lon
#   printVarSummary(lon)
#   printVarSummary(lat)
#   dimlat = dimsizes(dimsizes(lat))
#   dimlon = dimsizes(dimsizes(lon))

#   VAR1    = f1[:]->$VAR(0)$(ind_time_sel,:,:)

#   printVarSummary(VAR1)
#   dimVAR1 = dimsizes(VAR1)

#   DTrend     = dtrend_msg_n(timeUTsel,VAR1(:,:,:),True,True,0)
#   slope2D    = onedtond(DTrend@slope, (/dimVAR1(1),dimVAR1(2)/) )
#   printVarSummary(slope2D)

#     # Adjust attibutes of the dataset that is returned
#     Name_NetCDF = "TrendZOS_ForEXP"+ForEXP+".nc"
#     system("/bin/rm -f "+Name_NetCDF)    ; remove any pre-existing file
#     ncdf = addfile(Name_NetCDF ,"c")  ; open output netCDF file

#     fAtt               = True            ; assign file attributes
#     fAtt@title         = "Storage of linear trend of zos between "+year_min+"-"+year_max+ \
#                          " for PlotThermalExpMaps.ncl script. Units are in cm/year"
#     fAtt@creation_date = systemfunc ("date")
#     fileattdef( ncdf, fAtt )            ; copy file attributes
#     ;### Export in NetCDF file
#     ncdf->$Models(i)$ = slope2D  ; Unit in m/year