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

MV = 1e20
smallNumber = 1E-39

month = int(sys.argv[1])

def returnForecast(dateInput, varName, doy):
  outPut = np.zeros((3, 180, 360))
  modTel = 0
  for model in ["CanCM3", "CanCM4", "FLOR"]:
    if model == "CanCM3":
      dirLoc = "/tigress/nwanders/Scripts/Seasonal/CanCM3/"
      ensNr = 10
      if varName == "tas":
        factor = 1.
      else:
        varName = "prlr"
        factor = 1000.
    if model == "CanCM4":
      dirLoc = "/tigress/nwanders/Scripts/Seasonal/CanCM4/"
      ensNr = 10
      if varName == "tas":
        factor = 1.
      else:
        varName = "prlr"
        factor = 1000.
    if model == "FLOR":
      dirLoc = "/tigress/nwanders/Scripts/Seasonal/FLOR/"
      ensNr = 12
      if varName == "tas":
        factor = 1.
      else:
        varName = "pr"
        factor = 1000.
    if model == "CCSM":
      dirLoc = "/tigress/nwanders/Scripts/Seasonal/CCSM/"
      ensNr = 10
      if varName == "tas":
        factor = 1.
      else:
        varName = "prec"
        factor = 1.
    out = np.zeros((ensNr, 180, 360))
    for ens in range(ensNr):
      out[ens,:,:] = readNC(ncFileName(model, varName, dateInput, ensNumber = ens+1),varName, DOY=doy) * factor
    outPut[modTel,:,:] = np.mean(out,axis=0)
    modTel += 1
  return(outPut)


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
    if model == "CCSM":
        modelName = "_day_CCSM4_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"01_r"+str(ensNumber)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    return(ncFile)

def ncFileNameEns(model, varName, startTime, endTime, ens, dirLoc = '/tigress/nwanders/Scripts/Seasonal/'):
    if model == "FLOR":
        modelName = "_day_GFDL-FLORB01_FLORB01-P1-ECDA-v3.1-"
        ncFile = dirLoc+model+'/'+varName+modelName+monthYear(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM3":
        modelName = "_day_CanCM3_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM4":
        modelName = "_day_CanCM4_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CCSM":
        modelName = "_day_CCSM4_"
        ncFile = dirLoc+model+'/'+varName+modelName+yearMonth(startTime)+"01_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    return(ncFile)


def checkForcingFiles(model, precName, tempName, startTime, endTime):
    totEns = 10
    if model == "FLOR":
        totEns = 12
    ens = 0
    fileExist = True
    while ens < totEns and fileExist:
        print ncFileNameEns(model, precName, startTime, endTime, ens+1)
        fileExist = os.path.isfile(ncFileNameEns(model, precName, startTime, endTime, ens+1)) \
           and os.path.isfile(ncFileNameEns(model, tempName, startTime, endTime, ens+1))
        ens += 1
    return fileExist

def readNC(ncFile,varName, DOY=1):
    
    # Get netCDF file and variable name:
    f = nc.Dataset(ncFile)    
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
    
    #print ncFile
    #f = nc.Dataset(ncFile)  
    varName = str(varName)
    if DOY == "all":
      outputData = f.variables[varName][:,:,:,:]       # still original data    
    else:
      outputData = f.variables[varName][DOY,:,:,:]       # still original data
    
    f.close()
    
    return(outputData)

def calcEnsVar(covMatrix, coef):
  totalVar = 0
  for i in range(covMatrix.shape[0]):
    for j in range(covMatrix.shape[0]):
      totalVar += coef[i]*coef[j]*abs(covMatrix[i,j])
  return(totalVar)


dims = {}
dims['minlat'] = -89.5
dims['minlon'] = 0.5
dims['nlat'] = 180 
dims['nlon'] = 360
dims['res'] = 1.000
dims['maxlat'] = dims['minlat'] + dims['res']*(dims['nlat']-1)
dims['maxlon'] = dims['minlon'] + dims['res']*(dims['nlon']-1)

varNames = ["prec","tas", "prec_var", "tas_var"]

modelS = ["CanCM3", "CanCM4", "FLOR"]

for year in range(1981,2013):
    d = 0
    forecastDate = datetime.datetime(year, month, 1)
    endTime = datetime.datetime(year+1, month, 1) - datetime.timedelta(days=d)
    day = forecastDate + datetime.timedelta(days=d)
    subTel = 0
    ncFile = '/tigress/nwanders/Scripts/Seasonal/WeightedEqual/%d%.2d%.2d_forecasts_CanCM3_CanCM4_FLOR.nc' %(forecastDate.year, forecastDate.month, forecastDate.day)
    createNetCDF(ncFile, varNames, ["mm","C", "mm","C"], latitudes=np.arange(89.5,-90,-1), longitudes=np.arange(-179.5,180,1), loop=True)
    weightFilePrec = '/tigress/nwanders/Scripts/Seasonal/resultsNetCDF/NNLS_weights_month_%d_var_prec.nc4' %(month)
    weightFileTemp = '/tigress/nwanders/Scripts/Seasonal/resultsNetCDF/NNLS_weights_month_%d_var_tas.nc4' %(month)
    while endTime > day:
      print day
      if day.day == 1 or day.day == 16:
        precNC = np.zeros((16,3,180,360))
        tempNC = np.zeros((16,3, 180,360))
        dayCount = 0
        weightsPrec = np.zeros((3,180,360))+1./3.
        refPrec = np.zeros((3,180,360))
        refPrec[0,:,:] = readNC(weightFilePrec, "CanCM3ref", subTel)
        refPrec[1,:,:] = readNC(weightFilePrec, "CanCM4ref", subTel)
        refPrec[2,:,:] = readNC(weightFilePrec, "FLORref", subTel)
        refLevelPrec = readNC(weightFilePrec, "PGFref", subTel)
        weightsTemp = np.zeros((3,180,360))+1./3.
        refTemp = np.zeros((3,180,360))
        refTemp[0,:,:] = readNC(weightFileTemp, "CanCM3ref", subTel)
        refTemp[1,:,:] = readNC(weightFileTemp, "CanCM4ref", subTel)
        refTemp[2,:,:] = readNC(weightFileTemp, "FLORref", subTel)
        refLevelTemp = readNC(weightFileTemp, "PGFref", subTel)
        subTel += 1
      try:
        precNC[dayCount,:,:,:] = returnForecast(forecastDate, "prec", d)
        tempNC[dayCount,:,:,:] = returnForecast(forecastDate, "tas", d)
      except:
        print "Leap year"
      if  (day + datetime.timedelta(days=1)).day == 1 or (day + datetime.timedelta(days=1)).day == 16:
        precEns = np.mean(precNC[0:(dayCount+1),:,:,:], axis=0)
        tempEns = np.mean(tempNC[0:(dayCount+1),:,:,:], axis=0)
        precAnomaly = np.sum((precEns - refPrec)*weightsPrec, axis=0)
        tempAnomaly = np.sum((tempEns - refTemp)*weightsTemp, axis=0)
        precNew = np.minimum(np.maximum(np.sum(((precNC[0:(dayCount+1),:,:,:]/np.sum(precNC[0:(dayCount+1),:,:,:], axis=0))*(weightsPrec/np.sum(weightsPrec, axis=0))) * (precAnomaly + refLevelPrec), axis=1) * (dayCount+1),0.0),100.)
        tempNew = np.maximum(np.sum((tempNC[0:(dayCount+1),:,:,:]-np.mean(tempNC[0:(dayCount+1),:,:,:], axis=0))*(weightsTemp/np.sum(weightsTemp, axis=0)), axis=1) + tempAnomaly + refLevelTemp, np.min(refLevelTemp))
        varPrec = np.zeros((180,360))
        varTemp = np.zeros((180,360))
        for y in range(180):
          for x in range(360):
            varPrec[y,x] = calcEnsVar(np.cov(precNC[0:(dayCount+1),:,y,x], rowvar=0), weightsPrec[:,y,x])
            varTemp[y,x] = calcEnsVar(np.cov(tempNC[0:(dayCount+1),:,y,x], rowvar=0), weightsTemp[:,y,x])        
        for i in range(dayCount+1):
          print forecastDate+datetime.timedelta(days=d-dayCount+i)
          print d-dayCount+i
          out = precNew[i,:,:]
          out[np.isnan(out)] = refLevelPrec[np.isnan(out)]
          data2NetCDF(ncFile, varNames[0], precNew[i,:,:], forecastDate+datetime.timedelta(days=d-dayCount+i), posCnt = d-dayCount+i)
          data2NetCDF(ncFile, varNames[2], varPrec, forecastDate+datetime.timedelta(days=d-dayCount+i), posCnt = d-dayCount+i)
          out = tempNew[i,:,:]
          out[np.isnan(out)] = refLevelTemp[np.isnan(out)]
          data2NetCDF(ncFile, varNames[1], tempNew[i,:,:], forecastDate+datetime.timedelta(days=d-dayCount+i), posCnt = d-dayCount+i)
          data2NetCDF(ncFile, varNames[3], varTemp, forecastDate+datetime.timedelta(days=d-dayCount+i), posCnt = d-dayCount+i)
      print time.strftime("%H:%M:%S")
      d += 1
      dayCount += 1
      day = forecastDate + datetime.timedelta(days=d)
      #Close the outgoing file

