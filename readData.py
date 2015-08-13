import netCDF4 as nc
import datetime

import numpy as np
import numpy.random as random

# Global variables:
MV = 1e20
smallNumber = 1E-39

# file cache to minimize/reduce opening/closing files.  

def readNC(ncFile,varName, dateInput, latPoint = None, lonPoint = None, endDay = None, useDoy = None, LatitudeLongitude = False, specificFillValue = None, model = "NMME"):
    
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
    elif endDay != None:
        if isinstance(date, str) == True and isinstance(endDay, str) == True:
            startDay = datetime.datetime.strptime(str(date),'%Y-%m-%d')
            lastDay = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
        dateDif = datetime.datetime(startDay.year,startDay.month,startDay.day) - orgDate
        deltaDays = datetime.datetime(lastDay.year,lastDay.month,lastDay.day) - orgDateEnd + datetime.timedelta(days=1)
        # time index (in the netCDF file)
        nctime = f.variables['time']  # A netCDF time variable object.
        print nctime[:]
        print int(dateDif.days)
        print int(deltaDays.days)-1
        print np.where(nctime[:] >= int(dateDif.days))[0][0]
        print np.where(nctime[:] <= int(deltaDays.days)-1)[0][-1]
        idx = range(int(np.where(nctime[:] >= int(dateDif.days))[0][0]), int(np.where(nctime[:] <= int(deltaDays.days)-1)[0][-1])+1)
    else:
        if isinstance(date, str) == True:
            date = datetime.datetime.strptime(str(date),'%Y-%m-%d') 
        dateDif = datetime.datetime(date.year,date.month,date.day) - orgDate
        # time index (in the netCDF file)
        nctime = f.variables['time']  # A netCDF time variable object.
        idx = int(np.where(nctime[:] == int(dateDif.days))[0])
    
    outputData = f.variables[varName][idx,:,:]       # still original data
    outputData[outputData == MV] = np.nan
    f.close()
    
    return(outputData)

def createNetCDF(ncFileName, varName, varUnits, latitudes, longitudes,\
                                      longName = None, loop=False):
    
    rootgrp= nc.Dataset(ncFileName,'w', format='NETCDF4')
    
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


def leapCount(start):
    tempStart = datetime.datetime.strptime(str(start),'%Y-%m-%d')
    if tempStart.month <= 2:
        leapDif = np.floor(((tempStart.year-1) - 1852)/4.)
    else:
        leapDif = np.floor((tempStart.year - 1852)/4.)
    return(leapDif)

def findMonthEnd(year, month, day, model):
    if month + 1 <= 12:
        firstDay = datetime.datetime(year, month+1,1)
    else:
        firstDay = datetime.datetime(year+1, month-11,1)
    out =firstDay-datetime.timedelta (days = 1)
    if day > 25: day=31
    if out.day > day:
        out = datetime.datetime(year, month, day)
    if model == "CanCM3" or model == "CanCM4":
        if out.day == 29:
            out = datetime.datetime(year, month, 28)
    return(out)

def lagToDateStr(date, lag, model):
    try:
        startDate = datetime.datetime.strptime(date,'%Y-%m-%d')
    except:
        startDate = date
    y = startDate.year
    m = startDate.month
    d = startDate.day
    if m+lag > 12:
        tempEndDate = findMonthEnd(y+1,m-12+lag,d, model)
        zero = ""
        if len(str(m-11+lag)) < 2: zero = "0"
        zero2 = ""
        if len(str(startDate.day)) < 2: zero2 = "0"
        newDate = str(tempEndDate.year) +"-"+ zero + str(tempEndDate.month) + "-" + zero2+str(tempEndDate.day)
    else:
        tempEndDate = findMonthEnd(y,m+lag,d, model)        
        zero = ""
        if len(str(m+lag)) < 2: zero = "0"
        zero2 = ""
        if len(str(startDate.day)) < 2: zero2 = "0"
        newDate = str(tempEndDate.year) +"-"+ zero + str(tempEndDate.month) + "-" + zero2+str(tempEndDate.day)
    return(newDate)


def lagToDateTime(date, lag, model):
    try:
        startDate = datetime.datetime.strptime(date,'%Y-%m-%d')
    except:
        startDate = date
    y = startDate.year
    m = startDate.month
    d = startDate.day
    if m+lag > 12:
        tempEndDate = findMonthEnd(y+1,m-12+lag,d,model)
    else:
        tempEndDate = findMonthEnd(y,m+lag,d,model)
    return(tempEndDate)

def returnSeasonalForecast(dateInput, endDay, model, varName, lag, month = 0, ensNr = 1, dirLoc=""):
    if varName == "discharge":
      varFile = "discharge_seasoAvg_output.nc"
    model = "PGF"
    deltaDay = lagToDateTime(endDay, lag, model).day - lagToDateTime(dateInput, lag, model).day + 1
    deltaYear = lagToDateTime(endDay, lag, model).year - lagToDateTime(dateInput, lag, model).year + 1
    data = np.zeros((deltaYear,360,720))
    print data.shape
    start = datetime.datetime.strptime(str(dateInput),'%Y-%m-%d')
    end = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
    lastEntry = 0
    m = start.month
    for y in range(start.year, end.year+1):
        tempStartDate = datetime.datetime.strptime(str(str(y)+"-"+str(m)+"-01"),'%Y-%m-%d')
        zero = ""
        if len(str(m)) < 2: zero = "0"
        zeroDay = ""
        if len(str(start.day)) < 2: zeroDay="0"
        startDate = str(tempStartDate.year)+"-"+zero+str(tempStartDate.month)+"-"+zeroDay+str(start.day)
        tempEnd = lagToDateTime(findMonthEnd(y, m, 31, model),11, model)
        if tempStartDate.year >= start.year and tempStartDate < (end - datetime.timedelta (days = 1)):
            zero = ""
            if len(str(m)) < 2: zero = "0"
            zero2 = ""
            if len(str(tempEnd.month)) < 2: zero2 = "0"
            tempEndDate = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            #tempEndDate = lagToDateTime(str(y)+"-"+zero+str(m)+"-"+str(end.day), lag)
            zero = ""
            if len(str(tempEndDate.month)) < 2: zero = "0"
            startDateTime = lagToDateTime(startDate, lag, model)
            endDateTime = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            addYear = 0
            print endDateTime.month
            print startDateTime.month
            if endDateTime.year < startDateTime.year: addYear = 1
            endDateTime = lagToDateTime(findMonthEnd(y+addYear,end.month,end.day, model), lag, model)
            print startDate
            print findMonthEnd(y,end.month,end.day, model)
            print lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            print lag
            print endDateTime
            print lagToDateTime(startDate, lag, model)
            if endDateTime.month < startDateTime.month and endDateTime.year == startDateTime.year: addYear = 1
            endDate = lagToDateStr(findMonthEnd(y+addYear, end.month, end.day, model), lag, model)
            deltaDay = (datetime.datetime.strptime(endDate,'%Y-%m-%d')-datetime.datetime.strptime(lagToDateStr(startDate, lag, model),'%Y-%m-%d')).days + 1
            for ens in range(ensNr):
                zero = ""
                if len(str(m)) < 2: zero = "0"
                ncFile = dirLoc+str(y)+"-"+zero+str(m)+"-01/netcdf/"+str(ens+1)+"/"+varFile
                print ncFile
                print lagToDateStr(startDate, lag, model)
                print endDate
                if ens == 0:
                    temp = readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model)
                    tempData = np.zeros((ensNr, 360,720))
                    tempData[ens,:,:] = aggregateTime(temp, var=varName)
                else:
                    tempData[ens,:,:] = aggregateTime(readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model), var=varName)
            data[lastEntry,:,:] = ensembleMean(tempData)
            lastEntry += 1
    return(data)



def aggregateTime(data, timeDimension = 0, var="prec"):
    if var == "tas":
        outPut = data.mean(axis=timeDimension)
    else:
        outPut = data.sum(axis=timeDimension)
    return(outPut)

def ensembleMean(data, ensDimension = 0):
    outPut = data.mean(axis=ensDimension)
    return(outPut)

def aggregateSpace(data, extent = 0):
    if extent != 0:
        outPut = data
        dy = data.shape[1]
        dx = data.shape[2]
        largerData = np.zeros(np.add(data.shape, (0,extent*2,extent*2)))
        yLen = range(extent,(dy+extent))
        xLen = range(extent,(dx+extent))
        largerData[:,extent:(dy+extent), extent:(dx+extent)] = data
        largerData[:,yLen,0:extent] = data[:,:,dx-extent:dx]
        largerData[:,yLen,dx:(dx+extent)] = data[:,:,0:extent]
        largerData[:,0:extent,xLen] = data[:,dy-extent:dy,:]
        largerData[:,dy:(dy+extent),xLen] = data[:,0:extent,:]
        for x in range(-extent, extent+1):
            for y in range(-extent, extent+1):
                if x == 0 and y == 0:
                    pass
                else:
                    outPut += largerData[:,extent+y:dy+extent+y, extent+x:dx+extent+x]
        return(outPut/(extent**2))
    else:
        return(data)

def readForcing(ncFile, varName, dateInput, endDay, lag=0, model="PGF"):
    if varName == "discharge":
      varFile = "discharge_dailyTot_output.nc"
    deltaDay = lagToDateTime(endDay, lag, model).day - lagToDateTime(dateInput, lag, model).day + 1
    deltaYear = lagToDateTime(endDay, lag, model).year - lagToDateTime(dateInput, lag, model).year + 1
    data = np.zeros((deltaYear,360,720))
    print data.shape
    start = datetime.datetime.strptime(str(dateInput),'%Y-%m-%d')
    end = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
    lastEntry = 0
    m = start.month
    d = start.day
    for y in range(start.year, end.year+1):
        tempStartDate = datetime.datetime.strptime(str(str(y)+"-"+str(m)+"-"+str(d)),'%Y-%m-%d')
        zero = ""
        if len(str(m)) < 2: zero = "0"
        zeroDay = ""
        if len(str(start.day)) < 2: zeroDay="0"
        startDate = str(tempStartDate.year)+"-"+zero+str(tempStartDate.month)+"-"+zeroDay+str(start.day)
        tempEnd = datetime.datetime.strptime(str(str(y+1)+"-"+str(m)+"-01"),'%Y-%m-%d') - datetime.timedelta (days = 1)
        print tempStartDate
        if tempStartDate >= start and tempStartDate < (end - datetime.timedelta (days = 1)):
            startDateTime = lagToDateTime(startDate, lag, model)
            endDateTime = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            addYear = 0
            if endDateTime.year < startDateTime.year: addYear = 1
            endDateTime = lagToDateTime(findMonthEnd(y+addYear,end.month,end.day, model), lag, model)
            if endDateTime.month < startDateTime.month and endDateTime.year == startDateTime.year: addYear = 1
            endDate = lagToDateStr(findMonthEnd(y+addYear, end.month, end.day, model), lag, model)
            print lagToDateStr(startDate, lag, model)
            print endDate
            data[lastEntry,:,:] = aggregateTime(readNC(ncFile, varName, lagToDateStr(startDate, lag, model), endDay=endDate, model="PGF"), var=varName)
            lastEntry += 1
    return(data)


def readRandomForcing(ncFile, varName, dateInput, endDay, lag=0, model="PGF", ensNr = 10):
    deltaDay = lagToDateTime(endDay, lag, model).day - lagToDateTime(dateInput, lag, model).day + 1
    deltaYear = lagToDateTime(endDay, lag, model).year - lagToDateTime(dateInput, lag, model).year + 1
    data = np.zeros((deltaYear,180,360))
    print data.shape
    start = datetime.datetime.strptime(str(dateInput),'%Y-%m-%d')
    end = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
    lastEntry = 0
    m = start.month
    d = start.day
    yearS = range(start.year, end.year+1)
    for y in range(start.year, end.year+1):
        realY = y
        tempYears = np.delete(yearS, realY-yearS[0])
        randomYears = random.choice(tempYears, ensNr, replace=False)
        tempData = np.zeros((ensNr, 180,360))
        for ens in range(ensNr):
            y = randomYears[ens]
            tempStartDate = datetime.datetime.strptime(str(str(y)+"-"+str(m)+"-01"),'%Y-%m-%d')
            zero = ""
            if len(str(m)) < 2: zero = "0"
            zeroDay = ""
            if len(str(start.day)) < 2: zeroDay="0"
            startDate = str(tempStartDate.year)+"-"+zero+str(tempStartDate.month)+"-"+zeroDay+str(start.day)
            tempEnd = datetime.datetime.strptime(str(str(y+1)+"-"+str(m)+"-01"),'%Y-%m-%d') - datetime.timedelta (days = 1)
            if tempStartDate >= start and tempStartDate < (end - datetime.timedelta (days = 1)):
                startDateTime = lagToDateTime(startDate, lag, model)
                endDateTime = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
                addYear = 0
                if endDateTime.year < startDateTime.year: addYear = 1
                endDateTime = lagToDateTime(findMonthEnd(y+addYear,end.month,end.day, model), lag, model)
                if endDateTime.month < startDateTime.month and endDateTime.year == startDateTime.year: addYear = 1
                endDate = lagToDateStr(findMonthEnd(y+addYear, end.month, end.day, model), lag, model)
                print lagToDateStr(startDate, lag, model)
                tempData[ens,:,:] = aggregateTime(readNC(ncFile, varName, lagToDateStr(startDate, lag, model), endDay=endDate, model="PGF"), var=varName)
        data[lastEntry,:,:] = ensembleMean(tempData)
        lastEntry += 1
    return(data)

def returnSeasonalEnsembleForecast(dateInput, endDay, model, varName, lag, month = 0, ensNr = 1, dirLoc=""):
    deltaDay = lagToDateTime(endDay, lag, model).day - lagToDateTime(dateInput, lag, model).day + 1
    deltaYear = lagToDateTime(endDay, lag, model).year - lagToDateTime(dateInput, lag, model).year + 1
    data = np.zeros((deltaYear,ensNr,180,360))
    print data.shape
    start = datetime.datetime.strptime(str(dateInput),'%Y-%m-%d')
    end = datetime.datetime.strptime(str(endDay),'%Y-%m-%d')
    lastEntry = 0
    m = start.month
    for y in range(start.year, end.year+1):
        tempStartDate = datetime.datetime.strptime(str(str(y)+"-"+str(m)+"-01"),'%Y-%m-%d')
        zero = ""
        if len(str(m)) < 2: zero = "0"
        zeroDay = ""
        if len(str(start.day)) < 2: zeroDay="0"
        startDate = str(tempStartDate.year)+"-"+zero+str(tempStartDate.month)+"-"+zeroDay+str(start.day)
        tempEnd = lagToDateTime(findMonthEnd(y, m, 31, model),11, model)
        if tempStartDate.year >= start.year and tempStartDate < (end - datetime.timedelta (days = 1)):
            zero = ""
            if len(str(m)) < 2: zero = "0"
            zero2 = ""
            if len(str(tempEnd.month)) < 2: zero2 = "0"
            tempEndDate = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            #tempEndDate = lagToDateTime(str(y)+"-"+zero+str(m)+"-"+str(end.day), lag)
            zero = ""
            if len(str(tempEndDate.month)) < 2: zero = "0"
            startDateTime = lagToDateTime(startDate, lag, model)
            endDateTime = lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            addYear = 0
            print endDateTime.month
            print startDateTime.month
            if endDateTime.year < startDateTime.year: addYear = 1
            endDateTime = lagToDateTime(findMonthEnd(y+addYear,end.month,end.day, model), lag, model)
            print startDate
            print findMonthEnd(y,end.month,end.day, model)
            print lagToDateTime(findMonthEnd(y,end.month,end.day, model), lag, model)
            print lag
            print endDateTime
            print lagToDateTime(startDate, lag, model)
            if endDateTime.month < startDateTime.month and endDateTime.year == startDateTime.year: addYear = 1
            endDate = lagToDateStr(findMonthEnd(y+addYear, end.month, end.day, model), lag, model)
            deltaDay = (datetime.datetime.strptime(endDate,'%Y-%m-%d')-datetime.datetime.strptime(lagToDateStr(startDate, lag, model),'%Y-%m-%d')).days + 1
            for ens in range(ensNr):
                if model == "FLOR":
                    zero = ""
                    if len(str(m)) < 2: zero = "0"
                    ncFile = dirLoc+varName+"_day_GFDL-FLORB01_FLORB01-P1-ECDA-v3.1-"+zero+str(m)+str(y)+"_r"+str(ens+1)+"i1p1_"+str(y)+zero+str(m)+"01-"+str(tempEnd.year)+zero2+str(tempEnd.month)+str(tempEnd.day)+".nc4"
                    print ncFile
                    print lagToDateStr(startDate, lag, model)
                    print endDate
                    if ens == 0:
                        temp = readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model)
                        tempData = np.zeros((ensNr, 180,360))
                        tempData[ens,:,:] = aggregateTime(temp, var=varName)
                    else:
                        tempData[ens,:,:] = aggregateTime(readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model), var=varName)
                if model == "CanCM3":
                    zero = ""
                    if len(str(m)) < 2: zero = "0"
                    ncFile = dirLoc+varName+"_day_"+model+"_"+str(y)+zero+str(m)+"_r"+str(ens+1)+"i1p1_"+str(y)+zero+str(m)+"01-"+str(tempEnd.year)+zero2+str(tempEnd.month)+str(tempEnd.day)+".nc4"
                    print ncFile
                    print lagToDateStr(startDate, lag, model)
                    print endDate
                    if ens == 0:
                        temp = readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model)
                        tempData = np.zeros((ensNr, 180,360))
                        tempData[ens,:,:] = aggregateTime(temp, var=varName)
                    else:
                        tempData[ens,:,:] = aggregateTime(readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model), var=varName)
                if model == "CanCM4":
                    zero = ""
                    if len(str(m)) < 2: zero = "0"
                    ncFile = dirLoc+varName+"_day_"+model+"_"+str(y)+zero+str(m)+"_r"+str(ens+1)+"i1p1_"+str(y)+zero+str(m)+"01-"+str(tempEnd.year)+zero2+str(tempEnd.month)+str(tempEnd.day)+".nc4"
                    print ncFile
                    print lagToDateStr(startDate, lag, model)
                    print endDate
                    if ens == 0:
                        temp = readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model)
                        tempData = np.zeros((ensNr, 180,360))
                        tempData[ens,:,:] = aggregateTime(temp, var=varName)
                    else:
                        tempData[ens,:,:] = aggregateTime(readNC(ncFile,varName, lagToDateStr(startDate, lag, model), endDay = endDate, model=model), var=varName)
            data[lastEntry,:,:,:] = tempData
            lastEntry += 1
    return(data)

def nashSutcliffe(obs, mod):
  MSE = np.sum((obs-mod)**2)
  MSEclim = np.sum((obs - np.mean(obs))**2)
  NSE = 1-MSE/MSEclim
  return np.maximum(NSE,-999.)

def RMSE(obs, mod):
  out = np.mean((obs-mod)**2)**0.5
  return out


