import datetime
import os
import grads
import numpy as np
import fileinput
import netCDF4 as nc
import subprocess
import time
import random
import sys
import cPickle as pickle
import shutil
from calendar import monthrange 
#import matplotlib as mpl
#mpl.use('GTKAgg')
#import matplotlib.pyplot as plt
grads_exe = '/tigress/nwanders/Programs/opengrads-2.1.a2.oga.1.princeton/opengrads'
ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)

startTime = datetime.datetime(int(sys.argv[1][0:4]), int(sys.argv[1][4:6]), int(sys.argv[1][6:8]))
endTime = datetime.datetime(int(sys.argv[2][0:4]), int(sys.argv[2][4:6]), int(sys.argv[2][6:8]))
model = sys.argv[3]
refForcing = sys.argv[4]

print startTime
print endTime

MV = 1e20
smallNumber = 1E-39

def readNCforcing(ncFile,varName, dateInput, latPoint = None, lonPoint = None, endDay = None, useDoy = None, LatitudeLongitude = False, specificFillValue = None):

    # Get netCDF file and variable name:
    f = nc.Dataset(ncFile)
    print "New: ", ncFile

    #print ncFile
    #f = nc.Dataset(ncFile)  
    varName = str(varName)

    if LatitudeLongitude == True:
        try:
            f.variables['lat'] = f.variables['latitude']
            f.variables['lon'] = f.variables['longitude']
        except:
            pass
    orgDate = datetime.datetime(1901,1,1)
    orgDateEnd = orgDate
    date = dateInput
    if useDoy == "Yes":
        idx = dateInput - 1
    elif endDay != "None":
        if isinstance(date, str) == True and isinstance(endDay, str) == True:
            startDay = datetime.datetime.strptime(str(date),'%Y-%m-%d')
            lastDay = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
        dateDif = datetime.datetime(startDay.year,startDay.month,startDay.day) - orgDate
        deltaDays = datetime.datetime(lastDay.year,lastDay.month,lastDay.day) - orgDateEnd + datetime.timedelta(days=1)
        # time index (in the netCDF file)
        nctime = f.variables['time']  # A netCDF time variable object.
        if model == "FLOR" or model == "CCSM":
            #print int(np.where(nctime[:] == int(dateDif.days)+0.5)[0])
            #print deltaDays.days
            #print nctime[:]
            #print int(np.where(nctime[:] == int(deltaDays.days)-0.5)[0])
            print np.minimum(int(deltaDays.days)-0.5, np.max(nctime[:]))
            idx = range(int(np.where(nctime[:] == int(dateDif.days)+0.5)[0]), int(np.where(nctime[:] == np.minimum(int(deltaDays.days)-0.5, np.max(nctime[:])))[0])+1)
        else:
            idx = range(int(np.where(nctime[:] == int(dateDif.days))[0]), int(np.where(nctime[:] == int(deltaDays.days)-1)[0])+1)
    else:
        if isinstance(date, str) == True:
          date = datetime.datetime.strptime(str(date),'%Y-%m-%d')
        dateDif = datetime.datetime(date.year,date.month,date.day) - orgDate
        # time index (in the netCDF file)
        nctime = f.variables['time']  # A netCDF time variable object.
        idx = int(np.where(nctime[:] == int(dateDif.days))[0])

    outputData = f.variables[varName][idx,:,:]       # still original data

    f.close()

    return(outputData)


def createNetCDF(ncFileName, varName, varUnits, latitudes, longitudes,\
                                      longName = None, loop=False):
    
    rootgrp= nc.Dataset(ncFileName,'w', format= 'NETCDF4')
    
    #-create dimensions - time is unlimited, others are fixed
    rootgrp.createDimension('time',None)
    rootgrp.createDimension('lat',len(latitudes))
    rootgrp.createDimension('lon',len(longitudes))
    
    date_time= rootgrp.createVariable('time','f4',('time',))
    date_time.standard_name= 'time'
    date_time.long_name= 'Days since 1901-01-01'
    
    date_time.units= 'Days since 1901-01-01' 
    date_time.calendar= 'standard'
    
    lat= rootgrp.createVariable('lat','f4',('lat',))
    lat.long_name= 'latitude'
    lat.units= 'degrees_north'
    lat.standard_name = 'latitude'
    
    lon= rootgrp.createVariable('lon','f4',('lon',))
    lon.standard_name= 'longitude'
    lon.long_name= 'longitude'
    lon.units= 'degrees_east'
    
    lat[:]= latitudes
    lon[:]= longitudes
    
    if loop:
        for i in range(len(varName)):
            shortVarName = varName[i]
            longVarName  = varName[i]
            if longName != None: longVarName = longName
            var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=MV,zlib=True)
            var.standard_name = varName[i]
            var.long_name = longVarName
            var.units = varUnits[i]
    else:    
        shortVarName = varName
        longVarName  = varName
        if longName != None: longVarName = longName
        var= rootgrp.createVariable(shortVarName,'f4',('time','lat','lon',) ,fill_value=MV,zlib=True)
        var.standard_name = varName
        var.long_name = longVarName
        var.units = varUnits
    rootgrp.sync()
    rootgrp.close()

def readGrads(gradsfile,gradsVarName, gradsTime, lon=[0.5, 359.5], lat=[-89.5, 89.5]):
  ga = grads.GrADS(Bin=grads_exe,Window=False,Echo=False)
  
  ga("open " + gradsfile)
  
  if gradsTime != None:
    ga("set t " + gradsTime)
        
  ga("set lon "+ str(lon[0]) +" "+str(lon[1]))
  ga("set lat "+ str(lat[0]) +" "+str(lat[1]))
  longitudes = ga.coords().lon
  latitudes = ga.coords().lat
  time = ga.coords().denv.tyme[0]
  
  data=ga.exp(gradsVarName)
  data.longitudes = longitudes
  data.latitudes = latitudes
  data.time = time
  del(ga)
  return(data)

def data2NetCDF(ncFile,varName,varField,timeStamp,posCnt = None):
  #-write data to netCDF
  rootgrp= nc.Dataset(ncFile,'a')    
  
  shortVarName= varName        
  
  date_time= rootgrp.variables['time']
  if posCnt == None: posCnt = len(date_time)
  
  date_time[posCnt]= nc.date2num(timeStamp,date_time.units,date_time.calendar)
  rootgrp.variables[shortVarName][posCnt,:,:]= (varField)
  
  rootgrp.sync()
  rootgrp.close()


def datetime2gradstime(date):
 
 #Convert datetime to grads time
 str = date.strftime('%HZ%d%b%Y')
 
 return str

def gradstime2datetime(str):

 #Convert grads time to datetime
 date = datetime.datetime.strptime(str,'%HZ%d%b%Y')

 return date

def monthYear(time):
    return('%.2d%d' %(time.month, time.year))

def yearMonth(time):
    return('%d%.2d' %(time.year, time.month))

def yearMonthDay(time):
    return('%d%.2d%.2d' %(time.year, time.month, time.day))

def findMonthEnd(time, model):
    day = time.day
    month = time.month
    year = time.year
    out = datetime.datetime(year, month, day)
    if model == "CanCM3" or model == "CanCM4":
        if day == 29:
            out = datetime.datetime(year, month, 28)
    return(out)

def ncFileName(model, varName, startTime, ensNumber = 1, dirLoc = '/tigress/nwanders/Scripts/Seasonal/'):
    endTime = datetime.datetime(startTime.year+1, startTime.month, startTime.day) - datetime.timedelta(days = 1)
    if model == "FLOR":
        modelName = "_day_GFDL-FLORB01_FLORB01-P1-ECDA-v3.1-"
        ncFile = dirLoc+model+'/'+varName+modelName+monthYear(startTime)+"_r"+str(ensNumber)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM3":
        modelName = "_day_CanCM3_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"_r"+str(ensNumber)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM4":
        modelName = "_day_CanCM4_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"_r"+str(ensNumber)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    return(ncFile)   

def readNC(ncFile,varName, DOY=1):
    
    # Get netCDF file and variable name:
    f = nc.Dataset(ncFile)
    print "New: ", ncFile
    
    #print ncFile
    #f = nc.Dataset(ncFile)  
    varName = str(varName)
    if DOY == "all":
      outputData = f.variables[varName][:,:,:]       # still original data    
    else:
      outputData = f.variables[varName][DOY,:,:]       # still original data
    
    f.close()
    
    return(outputData)

def readNCMatching(ncFile,varName, DOY=1):
    
    # Get netCDF file and variable name:
    f = nc.Dataset(ncFile)
    print "New: ", ncFile
    
    #print ncFile
    #f = nc.Dataset(ncFile)  
    varName = str(varName)
    if DOY == "all":
      outputData = f.variables[varName][:,:,:,:]       # still original data    
    else:
      outputData = f.variables[varName][DOY,:,:,:]       # still original data
    
    f.close()
    
    return(outputData)


def matchCDF(data, orgDataCDF, refDataCDF, var="prec"):
    nx, ny, orgMax = orgDataCDF.shape
    refMax = refDataCDF.shape[2]
    if var == "prec":
        nonZero = orgDataCDF > 0.0
        nonZeroOrg = nonZero.argmax(2)-1
        
        nonZero = refDataCDF > 0.0
        nonZeroRef = nonZero.argmax(2)-1
        
    optDif = np.zeros((nx,ny))+9e9
    lowDif = np.zeros((nx,ny))
    highDif = np.zeros((nx,ny))
    optVal = np.zeros((nx,ny))+9e10
    out = np.zeros((nx,ny))
    
    for p in range(orgMax):
        absDif = np.abs(orgDataCDF[:,:,p] - data)
        improveVal = optDif > absDif
        optVal[improveVal] = p
        optDif[improveVal] = absDif[improveVal]
        lowDif[improveVal] = np.abs(orgDataCDF[:,:,np.max([p-1, 0])] - data)[improveVal]
        highDif[improveVal] = np.abs(orgDataCDF[:,:,np.min([p+1, orgMax-1])] - data)[improveVal]
    
    upInt = highDif+1e-05 < lowDif
    lowInt = highDif > lowDif+1e-05
    out = optVal
    out[upInt] = optVal[upInt] + optDif[upInt]/(highDif[upInt] + optDif[upInt])
    out[lowInt] = optVal[lowInt] - optDif[lowInt]/np.maximum(lowDif[lowInt] + optDif[lowInt], 1e-10)
    out = out/(orgMax-1.)
    out[out < 0.0] = 0.0
    out[out > 100.0] = 100.0
    
    transData= np.zeros((nx,ny))
    out = out * (refMax-1.)
    for p in range(refMax):
        selVal = np.floor(out) == p
        maxVal = refDataCDF[:,:,np.min([p+1,refMax-1])][selVal]
        minVal = refDataCDF[:,:,np.max([p,0])][selVal]
        lowData = minVal + (out[selVal] - p) * (maxVal - minVal)
        highData = minVal + ((p+1) - out[selVal]) * (maxVal - minVal)
        transData[selVal] = (lowData + highData)/2
        if p == 0 and var == "prec":
            randomRainChance = np.maximum(nonZeroOrg - nonZeroRef,0.0)/np.maximum(nonZeroOrg, 1e-10)
            randomRainChance[randomRainChance >= 0.99999] = 2
            randomRain = (np.random.random((nx,ny)) - (1-randomRainChance)) / np.maximum(randomRainChance, 1e-10)
            randomRain[randomRain < 0] = 0.0
            randomRain[randomRain > 1.0] = 0.0
            rainPercentile = (nonZeroRef-0.001) + np.ceil((nonZeroOrg - nonZeroRef) * randomRain)
            out[selVal] = rainPercentile[selVal]
    
    return(transData)

def Prepare_VIC_Global_Parameter_File(idate,fdate,dims, model, refForcing):
 try:
  os.mkdir('/tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/%d-%.2d-%.2d/' %(idate.year, idate.month, idate.day))
 except:
  foo = 0
 try:
  os.mkdir('/tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/%d-%.2d-%.2d/resultRAW' %(idate.year, idate.month, idate.day))
 except:
  foo = 0
 file = '/tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/controlFiles/%d%.2d%.2d.txt' %(idate.year, idate.month, idate.day)
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
 fp.write('FORCING1        /tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/forcing/forcing_\n')
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
 fp.write('RESULT_DIR      /tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/%d-%.2d-%.2d/resultRAW/\n' %(idate.year, idate.month, idate.day))
 fp.write('INIT_STATE /tigress/nwanders/Scripts/hydroSeasonal/'+refForcing+'/VIC/test/states/state_%d%.2d%.2d\n' % (idate.year,idate.month, idate.day))
 fp.write('STATENAME /tigress/nwanders/Scripts/hydroSeasonal/'+model+'/VIC/'+refForcing+'/states/state\n')
 fp.write('STATEYEAR %d \n' % fdate.year)
 fp.write('STATEMONTH %d \n' % fdate.month)
 fp.write('STATEDAY %d \n' % fdate.day) 
 
 #Close the file
 fp.close()
 return(file)

dims = {}
dims['minlat'] = -89.5
dims['minlon'] = 0.5
dims['nlat'] = 180 
dims['nlon'] = 360
dims['res'] = 1.000
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)

totEns = 10
year = startTime.year
randomYear = np.random.choice(np.concatenate((np.arange(1981,year),np.arange(year+1,2013))), totEns, replace=False)

for ens in range(totEns):
  forecastDate = startTime
  year = startTime.year
  month = startTime.month
  forcing_file = '/tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/forcing/forcing_%d%.2d01' %(year, month)
  fp = open(forcing_file,'wb')
  d = 0
  day = forecastDate + datetime.timedelta(days=d)
  randomDate = datetime.datetime(randomYear[ens],month, 1)
  randomDay = randomDate + datetime.timedelta(days=d)
  DOY = day - datetime.datetime(day.year, 1, 1)
  while endTime >= day:
    if day.day == 1:
      windNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/wind_clim_PGF.nc","wind", DOY=day.month-1)
    tmaxNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmax_diff_PGF.nc","tmax", DOY=DOY.days)
    tminNC = readNC("/tigress/nwanders/Scripts/Seasonal/refData/tmin_diff_PGF.nc","tmin", DOY=DOY.days)
    try:
      precNC = readNCforcing("/tigress/nwanders/Scripts/Seasonal/refData/prec_PGF_PCR.nc4","prec", dateInput=randomDay, endDay = "None")
      tempNC = readNCforcing("/tigress/nwanders/Scripts/Seasonal/refData/tas_PGF_PCR.nc4","tas", dateInput=randomDay, endDay = "None")
      tempMaxNC = tempNC + tmaxNC
      tempMinNC = tempNC - tminNC
      print randomDay
      print day
    except:
      print "Leap year"
    prec = np.zeros((180,360), dtype=np.float32)
    prec[:,0:180] = precNC[::-1,180:360]*1000.
    prec[:,180:360] = precNC[::-1,0:180]*1000.
    tmax = np.zeros((180,360), dtype=np.float32)
    tmax[:,0:180] = tempMaxNC[::-1,180:360]
    tmax[:,180:360] = tempMaxNC[::-1,0:180]
    tmin = np.zeros((180,360), dtype=np.float32)
    tmin[:,0:180] = tempMinNC[::-1,180:360]
    tmin[:,180:360] = tempMinNC[::-1,0:180]
    wind = np.array(windNC, dtype=np.float32)
    #Append to the outgoing file
    prec.tofile(fp)
    tmax.tofile(fp)
    tmin.tofile(fp)
    wind.tofile(fp)
    d += 1
    day = forecastDate + datetime.timedelta(days=d)
    randomDay = randomDate + datetime.timedelta(days=d)
    DOY = day - datetime.datetime(day.year, 1, 1)
  #Close the outgoing file
  fp.close()
  
  VIC_global = '/tigress/nwanders/Scripts/VIC/VIC_4.0.5_image_mode/VIC_dev.exe'
  
  settingsFile = Prepare_VIC_Global_Parameter_File(forecastDate,datetime.datetime(forecastDate.year, forecastDate.month, forecastDate.day) + datetime.timedelta(days = d),dims,model, refForcing)
  
  print time.strftime("%H:%M:%S")
  os.system(VIC_global + ' -g '+settingsFile+' >& /tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/logFiles/VIC_'+model+'%d%.2d_%d.txt' %(forecastDate.year, forecastDate.month, ens+1))
  print time.strftime("%H:%M:%S")
  
  varNames = ["prec","evap", "runoff", "baseflow", "sm1", "sm2", "sm3", "surf_temp", "swq", "snow_depth"]
  fileName = '/tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/%d-%.2d-%.2d/resultRAW/output_grid_%d%.2d%.2d00.ctl' %(forecastDate.year, forecastDate.month, forecastDate.day, forecastDate.year, forecastDate.month, forecastDate.day)
  try:
    os.mkdir('/tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/%d-%.2d-%.2d/%d' %(forecastDate.year, forecastDate.month, forecastDate.day,ens+1))
  except:
    foo = 0
  ncFile = '/tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/%d-%.2d-%.2d/%d/output_%d%.2d%.2d.nc' %(forecastDate.year, forecastDate.month, forecastDate.day,ens+1,forecastDate.year, forecastDate.month, forecastDate.day)
  data = readGrads(fileName, "prec", str(1), lon=[-179.5, 179.5])
  createNetCDF(ncFile, varNames, ["mm","mm","mm","mm","mm","mm","mm","C", "mm","cm"], latitudes=data.latitudes, longitudes=data.longitudes, loop=True)
  
  for i in range(d):
    print i
    for var in varNames:
      data = readGrads(fileName, var, str(i+1), lon=[-179.5, 179.5])
      data2NetCDF(ncFile, var, data, data.time, posCnt = i)
  print time.strftime("%H:%M:%S")
  os.remove(forcing_file)
  shutil.rmtree('/tigress/nwanders/Scripts/hydroSeasonal/ESP/VIC/'+refForcing+'/%d-%.2d-%.2d/resultRAW' %(forecastDate.year, forecastDate.month, forecastDate.day))
