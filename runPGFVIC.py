import datetime
import os
import grads
import numpy as np
import fileinput
import netCDF4 as nc
#import pyhdf.SD as sd
#import library_f90
import scipy.stats as ss
import subprocess
import dateutil.relativedelta as relativedelta
import time
import random
import cPickle as pickle
from calendar import monthrange 
#import matplotlib as mpl
#mpl.use('GTKAgg')
#import matplotlib.pyplot as plt
grads_exe = '/tigress/nwanders/Programs/opengrads-2.1.a2.oga.1.princeton/opengrads'
ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)

def datetime2gradstime(date):
 str = date.strftime('%HZ%d%b%Y')
 return str

def gradstime2datetime(str):
 date = datetime.datetime.strptime(str,'%HZ%d%b%Y')
 return date

def readNC(ncFile,varName, DOY=1):
    # Get netCDF file and variable name:
    f = nc.Dataset(ncFile)
    #print "New: ", ncFile
    
    #print ncFile
    #f = nc.Dataset(ncFile)  
    varName = str(varName)
    if DOY == "all":
      outputData = f.variables[varName][:,:,:]       # still original data    
    else:
      outputData = f.variables[varName][DOY,:,:]       # still original data
    
    f.close()
    
    return(outputData)


def Prepare_VIC_Global_Parameter_File(idate,fdate,dims,dataset):
 
 file = '/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/Global_Parameter.txt'
 fp = open(file,'w')
 
 #Write the VIC parameters to file
 
 # Define Global Parameters
 fp.write('NLAYER          3       # number of layers\n')
 fp.write('TIME_STEP       24       # model time step in hours (= 24 for water balance)\n')
 fp.write('STARTYEAR       %d      # year model simulation starts\n' % idate.year)
 fp.write('STARTMONTH      %d      # month model simulation starts\n' % idate.month)
 fp.write('STARTDAY        %d      # day model simulation starts\n' % idate.day)
 fp.write('STARTHOUR       0       # hour model simulation starts\n')
 fp.write('ENDYEAR         %d      # year model simulation ends\n' % fdate.year)
 fp.write('ENDMONTH        %d      # month model simulation ends\n' % fdate.month) 
 fp.write('ENDDAY          %d      # day model simulation ends\n' % fdate.day) 
 fp.write('SKIPYEAR        0       # no. of startup yrs to skip before writing output\n')
 fp.write('WIND_H          10.0    # height of wind speed measurement\n')
 fp.write('MEASURE_H       2.0     # height of humidity measurement\n')
 fp.write('NODES           10       # number of soil thermal nodes\n')
 fp.write('MAX_SNOW_TEMP   0.5     # maximum temperature at which snow can fall\n')
 fp.write('MIN_RAIN_TEMP   -0.5    # minimum temperature at which rain can fall\n')
 
 # Define Global Parameters
 fp.write('FULL_ENERGY     FALSE    # calculate full energy balance\n')
 #fp.write('FROZEN_SOIL     TRUE    # calculate frozen soils\n')
 fp.write('DIST_PRCP       TRUE        # use distributed precipitation\n')
 fp.write('COMPRESS        FALSE       # compress input and output files when done\n')
 fp.write('CORRPREC        FALSE       # correct precipitation for gauge undercatch\n')
 fp.write('GRID_DECIMAL    3           # number of decimals to use in gridded file names\n')
 fp.write('PRT_SNOW_BAND   FALSE   # print snow variables\n')
 fp.write('SNOW_STEP       1        # time step in hours to solve snow bands\n')
 fp.write('ROOT_ZONES      2               # number of root zones in veg parameter file\n')
 fp.write('BINARY_OUTPUT   TRUE   # default is ASCII, unless LDAS format\n')
 fp.write('MIN_WIND_SPEED  0.1     # minimum allowable wind speed\n')
 fp.write('PREC_EXPT       0.6             # fraction of grid cell receiving\n')
 fp.write('GRND_FLUX       FALSE # true for full energy, false for water balance\n')
 fp.write('QUICK_FLUX      FALSE   # true uses Liang (1999), false uses finite diff.\n')
 fp.write('NOFLUX          FALSE  # false uses const. T at damping depth\n')
 
 # Define (Meteorological) Forcing Files
 fp.write('FORCING1        /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/forcing/forcing_\n')
 fp.write('N_TYPES         4\n')
 fp.write('FORCE_TYPE      PREC    UNSIGNED 10\n')
 fp.write('FORCE_TYPE      TMAX    SIGNED  10\n')
 fp.write('FORCE_TYPE      TMIN    SIGNED  10\n')
 fp.write('FORCE_TYPE      WIND    UNSIGNED 10\n')
 fp.write('FORCE_FORMAT    BINARY \n')
 fp.write('FORCE_ENDIAN    LITTLE      # LITTLE for PC arch., BIG for Sun or HP-UX\n')
 fp.write('FORCE_DT        24            # time step of two input met files\n')
 fp.write('FORCEYEAR   %d   # year model meteorological forcing files start\n' % idate.year)
 fp.write('FORCEMONTH   %d   # month model meteorological forcing files start\n' % idate.month)
 fp.write('FORCEDAY        %d                   # day meteorological forcing files start\n' % idate.day)
 fp.write('FORCEHOUR       00                  # hour meteorological forcing files start\n')
 
 
 # INPUT and OUTPUT TYPE from PRINCETON  (mpan and lluo)
 fp.write('INPUT_GRID_DEF %d %d %.3f %.3f %.3f %.3f\n' % (dims['nlon'],dims['nlat'],dims['minlon'],dims['minlat'],dims['res'],dims['res']))
 fp.write('OUTPUT_GRID_DEF %d %d %.3f %.3f %.3f %.3f\n' % (dims['nlon'],dims['nlat'],dims['minlon'],dims['minlat'],dims['res'],dims['res']))
 #fp.write('INPUT_STEP_PER_FILE     1        # number of timesteps per input file\n')
 fp.write('OUTPUT_PER_STEP TRUE # number of timesteps per output file\n')
 #fp.write('INPUT_GZIP              FALSE       # true if forcing file gzipped\n')
 #fp.write('OUTPUT_GZIP             FALSE        # true for writing gzipped output\n')
 fp.write('GRID_INPUT          TRUE #true for reading the input in GrADS binary, default is false\n')
 fp.write('GRID_OUTPUT         TRUE  #true for writing the output in GrADS binary,default is false\n')
 fp.write('REGULAR_OUTPUT      FALSE  #true for writing the output in standard version, default is false\n')
 
 # Define Input and Output Data Files
 fp.write('SNOW_BAND       /tigress/nwanders/Scripts/VIC/1.0deg/global_lai.txt\n')
 fp.write('ARC_SOIL        FALSE   # read soil parameters from ARC/INFO ASCII grids\n')
 #fp.write('SOIL            ../DATA/VIC/INPUT/tmp.txt\n')
 fp.write('SOIL            /tigress/nwanders/Scripts/VIC/1.0deg/global_soils_combined_smoothed_3hourly_calib_depth2_0.7.txt\n')
 fp.write('VEGPARAM        /tigress/nwanders/Scripts/VIC/1.0deg/global_lai.txt\n')
 fp.write('VEGLIB          /tigress/nwanders/Scripts/VIC/veglib.dat\n')
 fp.write('GLOBAL_LAI      TRUE      # true if veg param file has monthly LAI\n')
 fp.write('RESULT_DIR      /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/output/\n')
 fp.write('INIT_STATE /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/states/state_%d%.2d%.2d\n' % (idate.year,idate.month, idate.day))
 fp.write('STATENAME /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/states/state\n')
 fp.write('STATEYEAR %d \n' % fdate.year)
 fp.write('STATEMONTH %d \n' % fdate.month)
 fp.write('STATEDAY %d \n' % fdate.day) 
 
 #Close the file
 fp.close()

dims = {}
dims['minlat'] = -89.5
dims['minlon'] = 0.5
dims['nlat'] = 180 
dims['nlon'] = 360
dims['res'] = 1.000
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)

for year in range(1978,1981):
  forcing_file = '/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/forcing/forcing_'+str(year)+'0101'
  fp = open(forcing_file,'wb')
  d = 0
  day = datetime.datetime(year, 1,1) + datetime.timedelta(days=d)
  while year == day.year:
    i = day - datetime.datetime(1948, 1,1)
    m = (year - 1948)*12 + day.month - 1
    precNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/prec_PGF.nc","prec", DOY=i.days)
    #tempNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tas_PGF.nc","tas", DOY=i.days)-273.15
    windNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/wind_PGF.nc","wind", DOY=day.month)
    tmaxNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmax_PGF.nc","tmax", DOY=i.days)-273.15
    tminNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmin_PGF.nc","tmin", DOY=i.days)-273.15
    prec = np.array(precNC, dtype=np.float32)
    tmax = np.array(tmaxNC, dtype=np.float32)
    tmin = np.array(tminNC, dtype=np.float32) 
    wind = np.array(windNC, dtype=np.float32)
    #Append to the outgoing file
    prec.tofile(fp)
    tmax.tofile(fp)
    tmin.tofile(fp)
    wind.tofile(fp)
    d += 1
    day = datetime.datetime(year, 1,1) + datetime.timedelta(days=d)
  #Close the outgoing file
  fp.close()
  
  VIC_global = '/tigress/nwanders/Scripts/VIC/VIC_4.0.5_image_mode/VIC_dev.exe'
  
  Prepare_VIC_Global_Parameter_File(datetime.datetime(year, 1,1),datetime.datetime(year+1, 1,1),dims,0)
  
  print time.strftime("%H:%M:%S")
  os.system(VIC_global + ' -g /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/Global_Parameter.txt >& /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/VIC_log.txt')
  print time.strftime("%H:%M:%S")

for year in range(1981,2015):
  for month in range(1,13):
    startDate = datetime.datetime(year, month,1)
    forcing_file = '/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/forcing/forcing_%d%.2d01' %(startDate.year, startDate.month)
    fp = open(forcing_file,'wb')
    d = 0
    day = datetime.datetime(year, month,1) + datetime.timedelta(days=d)
    while month == day.month:
      i = day - datetime.datetime(1948, 1,1)
      precNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/prec_PGF.nc","prec", DOY=i.days)
      #tempNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tas_PGF.nc","tas", DOY=i.days)-273.15
      windNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/wind_PGF.nc","wind", DOY=day.month)
      tmaxNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmax_PGF.nc","tmax", DOY=i.days)-273.15
      tminNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmin_PGF.nc","tmin", DOY=i.days)-273.15
      prec = np.array(precNC, dtype=np.float32)
      tmax = np.array(tmaxNC, dtype=np.float32)
      tmin = np.array(tminNC, dtype=np.float32)
      wind = np.array(windNC, dtype=np.float32)
      #Append to the outgoing file
      prec.tofile(fp)
      tmax.tofile(fp)
      tmin.tofile(fp)
      wind.tofile(fp)
      d += 1
      day = datetime.datetime(year, month,1) + datetime.timedelta(days=d)
    #Close the outgoing file
    fp.close()
    
    VIC_global = '/tigress/nwanders/Scripts/VIC/VIC_4.0.5_image_mode/VIC_dev.exe'
    
    if month != 12:
      Prepare_VIC_Global_Parameter_File(datetime.datetime(year, month,1),datetime.datetime(year, month+1,1),dims,0)
    else:
      Prepare_VIC_Global_Parameter_File(datetime.datetime(year, month,1),datetime.datetime(year+1, 1,1),dims,0)
    print time.strftime("%H:%M:%S")
    os.system(VIC_global + ' -g /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/Global_Parameter.txt >& /tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/test/VIC_log.txt')
    print time.strftime("%H:%M:%S")

