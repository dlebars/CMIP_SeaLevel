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

VAR = 'zos'
EXP = 'rcp85'

year_min_ref = 1986  # Included. Beginning of reference period
year_max_ref = 2006  # Excluded. End of reference period

Freq     = "mon"  # Frequency of time output: mon or fx (only to read inputs)
DataDir  = "/nobackup/users/bars/synda/cmip5/output1/"

DirOut   = "CorrectedZOS_TS_yearly"
ModelList = readAsciiTable("CMIP5modelSelection_"+EXP+"_"+VAR+".txt",1,"string",1)

;##### Start and end of each period 
; Sanne's project:
;years_s = (/ 1950, 1965, 1980, 1995, 2006, 2020, 2035 /)
;years_e = (/ 1965, 1980, 1995, 2006, 2020, 2035, 2050 /)
; Erwin's project:
years_s = ispan(2006,2100,1)
years_e = ispan(2006+1,2101,1)
; Star's project
;years_s = ispan(1900,2005,1)
;years_e = ispan(1900+1,2005+1,1)

mid = (years_e + years_s)/2
dimy = dimsizes(years_s)

;Read the regular 1*1 grid to use for regridded outputs
DIRgrid     = "/nobackup/users/bars/SeaLevelFromHylke/CMIP5_OCEAN/Fingerprints/"
fgrid       = addfile(DIRgrid+"Relative_icesheets.nc","r")
LonOut      = fgrid->longitude
LatOut      = fgrid->latitude
dimLonOut   = dimsizes(LonOut)
dimLatOut   = dimsizes(LatOut)
cLatOut     = cos(LatOut*rad) # rad not defined
;Build a mask for the new grid
MaskOut     = fgrid->DYN_ANT
MaskOut     = (/where(MaskOut.eq.0,1e+20,1)/)
MaskOut@_FillValue = 1e+20
; Mask the Caspian sea
MaskOut({35:50},{45:56}) = MaskOut@_FillValue
; New mask that includes the mediterranean
MaskOut_Med              = MaskOut
MaskOut_Med({20.5:41.5},{354:360}) = MaskOut_Med@_FillValue
MaskOut_Med({20.5:41.5},{0:44})    = MaskOut_Med@_FillValue
MaskOut_Med({41:47.5},{2:44})      = MaskOut_Med@_FillValue

delim    = " "
Centers  = rm_single_dims(str_get_field(ModelList, 1, delim))
Models   = rm_single_dims(str_get_field(ModelList, 2, delim))

print("Models used:")
print(Models)

dimMod   = dimsizes(Models)

fref   = addfile("ReferenceZOS_ForEXP"+EXP+"_"+year_min_ref+"_"+year_max_ref+".nc", "r")
ftrend = addfile("TrendZOS_ForEXP"+EXP+".nc", "r")

;Read the average zos fields to discount from the model sea level
;fzos_avg = addfile("CMIP5_SeaLevel_"+EXP+"_zos_avg_1950-2100.nc","r") ;For Sanne's project
fzos_avg = addfile("CMIP5_SeaLevel_"+EXP+"_zos_avg_"+year_min_ref+"-2100.nc","r")
;fzos_avg = addfile("CMIP5_SeaLevel_"+EXP+"_zos_avg_1900-2006.nc","r") ;For Star's historical files

zos_avg_ModelNames = tostring(fzos_avg->ModelNames)
zos_avg  = fzos_avg->AverageSeaLevel
time_zos_avg = fzos_avg->time
print("Check the time vector:")
print(time_zos_avg)
indzosref   = ind((time_zos_avg.ge.year_min_ref).and.(time_zos_avg.lt.year_max_ref))
zos_avg_ref =  dim_avg_n(zos_avg(:,indzosref),1)

do i=30,30 ;0,dimMod-1
  print("####### Working on model "+i+","+Models(i)+"  #####")

  ; #### Read scenario data
  file_name = DataDir+Centers(i)+"/"+Models(i)+"/"+EXP+"/"+Freq
  files1 = systemfunc("ls "+file_name+"/*/*/*/*/"+VAR(0)+"/*"+VAR(0)+"*.nc")
  ;Use last version of data:
  vs1    = str_get_field(files1, 14,"/")
  delete(files1)
  files1 = systemfunc("ls "+file_name+"/*/*/*/"+vs1(dimsizes(vs1)-1)+"/"+VAR+"/*"+VAR+"*.nc")
  print("#### Using following files: ####")
  print(files1)
  f1      = addfiles(files1,"r")
  f1s     = addfile(files1(0),"r") ; Need this trick to avoid issues of adding
                                   ; lon or lat coordinate arrays from addfiles
  time1  = f1[:]->time

  if time1@calendar.eq."proleptic_gregorian" then
    time1@calendar = "gregorian"
  end if
  timeUT = cd_calendar(time1, 4)

  ; #### Read historical data
  file_name  = DataDir+Centers(i)+"/"+Models(i)+"/historical/"+Freq
  files2    = systemfunc("ls "+file_name+"/*/*/*/*/"+VAR+"/*"+VAR+"*.nc")
  vs2        = str_get_field(files2, 14,"/")
  delete(files2)
  files2 = systemfunc("ls "+file_name+"/*/*/*/"+vs2(dimsizes(vs2)-1)+"/"+VAR+"/*"+VAR+"*.nc")
  print("### Also using these historical files: ###")
  print(files2)
  f2      = addfiles(files2,"r")
  time2   = f2[:]->time
  if time2@calendar.eq."proleptic_gregorian" then
    time2@calendar = "gregorian"
  end if
  timeUT2 = cd_calendar(time2, 4)

  lat    = f1s->lat
  lon    = f1s->lon
  printVarSummary(lon)
  printVarSummary(lat)
  dimlat = dimsizes(dimsizes(lat))
  dimlon = dimsizes(dimsizes(lon))

  if (Models(i).eq."bcc-csm1-1").or.(Models(i).eq."bcc-csm1-1-m").or. \
     (Models(i).eq."GFDL-ESM2G").or.(Models(i).eq."GFDL-ESM2M").or. \
     (Models(i).eq."GFDL-CM3") then
    lon = where(lon.lt.0,lon+360,lon)
  end if

  RefVAR1    = fref->$Models(i)$
  ; Use the reference field to mask the small seas that are not connected to the
  ; ocean and areas where sea ice is included on the ocean load
  RefVAR1_corr = RefVAR1 - dim_avg_n(dim_avg_n(RefVAR1,0),0)
  MaskRefVAR1  = where((RefVAR1_corr.ge.2).or.(RefVAR1_corr.le.-2),1e+20,1)
  MaskRefVAR1@_FillValue = 1e+20

  MAT_CorrectedZOS_reg = new((/dimy,dimLatOut,dimLonOut/),float)
    
  ;##### Loop on the years ######################################
  do y=0,dimy-1
    print("Workgin on period: "+years_s(y)+"-"+years_e(y))
    indzossel   = ind((time_zos_avg.ge.years_s(y)).and.(time_zos_avg.lt.years_e(y)))
    if ismissing(indzossel) then
      print("No data during the period: "+years_s(y)+"-"+years_e(y))
      MAT_CorrectedZOS_reg(y,:,:) = MAT_CorrectedZOS_reg@_FillValue
    else
      dim_indzossel = dimsizes(indzossel)
      if dim_indzossel.eq.1 then
        zos_avg_sel =  zos_avg(:,indzossel)
      else
        zos_avg_sel =  dim_avg_n(zos_avg(:,indzossel),1)
      end if

      if years_s(y).ge.2006 then
        ind_time_sel = ind((timeUT.ge.years_s(y)).and.(timeUT.le.years_e(y)))
        VAR1    = f1[:]->$VAR(0)$(ind_time_sel,:,:)
      else
        ind_time_sel = ind((timeUT2.ge.years_s(y)).and.(timeUT2.le.years_e(y)))
        VAR1    = f2[:]->$VAR(0)$(ind_time_sel,:,:)
      end if

      if (Models(i).eq."MIROC5").or.(Models(i).eq."GISS-E2-R").or. \
         (Models(i).eq."GISS-E2-R-CC").or.(Models(i).eq."EC-EARTH").or. \
         (Models(i).eq."MRI-CGCM3") then
;        VAR1@_FillValue = 0
         VAR1 = where(VAR1.eq.0,1e+20,VAR1)
      end if
      if Models(i).eq."EC-EARTH" then
         VAR1@_FillValue = 1e+20
      end if

      dimVAR1 = dimsizes(VAR1)
      VAR1avg = dim_avg_n_Wrap(VAR1,0) ; Compute time average
      ind_zos_avg = ind(zos_avg_ModelNames.eq.Models(i))
      AnomVAR1   = VAR1avg ; Read metadata
      AnomVAR1   = (/(VAR1avg - RefVAR1)/)*100
      AnomVAR1   = (/AnomVAR1*MaskRefVAR1/)

      ZOS_AVG_CORR = tofloat((zos_avg_sel(ind_zos_avg) - zos_avg_ref(ind_zos_avg))*100)

      ; Effective number of years to detrend: year of interest 
      ; minus mean of reference period
      nbyears = mid(y) - (year_max_ref+year_min_ref)/2

      TrendVAR1  = AnomVAR1 ; Read metadata
      TrendVAR1  = (/(ftrend->$Models(i)$)*nbyears/)*100

      DTrendVAR1 = AnomVAR1 ; Read metadata 
      DTrendVAR1 = (/AnomVAR1 - TrendVAR1 - ZOS_AVG_CORR/)
        
;Regrid to the reference 1*1 degree grid
      if dimlat.eq.1 then
        DTrendVAR1_reg = linint2(lon,lat,DTrendVAR1,True,LonOut,LatOut,0)
        else if dimlat.eq.2 then
           DTrendVAR1_reg = rcm2rgrid(lat,lon,DTrendVAR1,LatOut,LonOut,1)
        end if
      end if
      ;Can mask other problematic regions here
      DTrendVAR1_reg@_FillValue = 1e+20
      if (Models(i).eq."MIROC5").or.(Models(i).eq."GFDL-ESM2M").or. \
         (Models(i).eq."GISS-E2-R").or.(Models(i).eq."GISS-E2-R-CC") then
        DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut_Med
        else
          DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut
      end if
      poisson_grid_fill(DTrendVAR1_reg,True,1,100,1,0.5,0)
      DTrendVAR1_reg  = DTrendVAR1_reg*MaskOut

      area_mean       = wgt_areaave_Wrap(DTrendVAR1_reg, cLatOut, 1.0, 0)
      print("Removing area mean of:" + area_mean + " cm")

      MAT_CorrectedZOS_reg(y,:,:) = DTrendVAR1_reg - area_mean
      delete(indzossel)
      delete(VAR1)
      delete(VAR1avg)
      delete(ind_time_sel)
      delete(AnomVAR1)
      delete(TrendVAR1)
      delete(DTrendVAR1)
      delete(DTrendVAR1_reg)
      delete(area_mean)
    end if
  end do

  MAT_CorrectedZOS_reg!0 = "time"
  MAT_CorrectedZOS_reg!1 = "latitude"
  MAT_CorrectedZOS_reg!2 = "longitude"
  MAT_CorrectedZOS_reg&time = mid
  MAT_CorrectedZOS_reg&latitude = LatOut
  MAT_CorrectedZOS_reg&longitude = LonOut

  ;### Export in NetCDF file
  Name_NetCDF = DirOut+"/CorrectedZOS_EXP"+EXP+"_"+Models(i)+".nc"
  system("/bin/rm -f "+Name_NetCDF)    ; remove any pre-existing file
  ncdf = addfile(Name_NetCDF ,"c")     ; open output netCDF file

  fAtt               = True            ; assign file attributes
  fAtt@title         = "Storage of zos state" + \
                     " with reference the period "+year_min_ref+"-"+year_max_ref+ \
                     " corrected for pre-industrial control trend and for global zos change."
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes

  ncdf->CorrectedZOS_reg = MAT_CorrectedZOS_reg
  ncdf->year_min_ref     = year_min_ref
  ncdf->year_max_ref     = year_max_ref

  delete(vs1)
  delete(files1)
  delete(f1)
  delete(time1)
  delete(timeUT)
  delete(vs2)
  delete(files2)
  delete(f2)
  delete(time2)
  delete(timeUT2)
  delete(lat)
  delete(lon)
  delete(RefVAR1)
  delete(MaskRefVAR1)
  delete(RefVAR1_corr)
  delete(ZOS_AVG_CORR)
  delete(ncdf)
  delete(MAT_CorrectedZOS_reg)
end do