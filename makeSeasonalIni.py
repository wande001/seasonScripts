import datetime
import os.path

def timeToStr(tempDate):
    zero = ""
    if len(str(tempDate.month)) < 2: zero = "0"
    zero2 = ""
    if len(str(tempDate.day)) < 2: zero2 = "0"
    newDate = str(tempDate.year) +"-"+ zero + str(tempDate.month) + "-" + zero2+str(tempDate.day)
    return(newDate)

def findMonthEnd(time, model):
    day = time.day
    month = time.month
    year = time.year
    out = datetime.datetime(year, month, day)
    if model == "CanCM3" or model == "CanCM4":
        if day == 29:
            out = datetime.datetime(year, month, 28)
    return(out)


def monthYear(time):
    month = time.month
    year = time.year
    zero = ""
    if month < 10: zero = "0"
    return(zero+str(month)+str(year))

def yearMonth(time):
    month = time.month
    year = time.year
    zero = ""
    if month < 10: zero = "0"
    return(str(year)+zero+str(month))

def yearMonthDay(time):
    day = time.day
    month = time.month
    year = time.year
    zero = ""
    if month < 10: zero = "0"
    zero2 = ""
    if day < 10: zero2 = "0"
    return(str(year)+zero+str(month)+zero2+str(day))


def ncFileName(model, varName, startTime, endTime):
    if model == "FLOR":
        modelName = "_day_GFDL-FLORB01_FLORB01-P1-ECDA-v3.1-"
        ncFile = varName+modelName+monthYear(startTime)+"_r1i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM3":
        modelName = "_day_CanCM3_"
        ncFile = varName+modelName+yearMonth(startTime)+"_r1i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM4":
        modelName = "_day_CanCM4_"
        ncFile = varName+modelName+yearMonth(startTime)+"_r1i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CCSM":
        modelName = "_day_CCSM4_"
        ncFile = varName+modelName+yearMonth(startTime)+"01_r1i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    return(ncFile)

def ncFileNameEns(model, varName, startTime, endTime, ens):
    if model == "FLOR":
        modelName = "_day_GFDL-FLORB01_FLORB01-P1-ECDA-v3.1-"
        ncFile = varName+modelName+monthYear(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM3":
        modelName = "_day_CanCM3_"
        ncFile = varName+modelName+yearMonth(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CanCM4":
        modelName = "_day_CanCM4_"
        ncFile = varName+modelName+yearMonth(startTime)+"_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    if model == "CCSM":
        modelName = "_day_CCSM4_"
        ncFile = varName+modelName+yearMonth(startTime)+"01_r"+str(ens)+"i1p1_"+yearMonthDay(startTime)+"-"+yearMonthDay(findMonthEnd(endTime, model))+".nc4"
    return(ncFile)


def checkForcingFiles(path, model, precName, tempName, startTime, endTime):
    totEns = 10
    if model == "FLOR":
        totEns = 12
    ens = 0
    fileExist = True
    while ens < totEns and fileExist:
        fileExist = os.path.isfile(path+ncFileNameEns(model, precName, startTime, endTime, ens+1)) \
           and os.path.isfile(path+ncFileNameEns(model, tempName, startTime, endTime, ens+1))
        ens += 1
    return fileExist
    

beginYear = 1981
endYear = 2011
forcingInput = "/tigress/nwanders/Scripts/Seasonal/"
modelS = ["CCSM"] #["CanCM3","CanCM4", "FLOR"]
precVarNameS = ["prec"] #["prlr","prlr","pr"]
tempVarName = "tas"
refInput = "/tigress/nwanders/Scripts/Seasonal/refData/"
refModelS = ["PGF"] #["PGF", "CFS"]
precRefVarName = "prec"
tempRefVarName = "tas"
precCorFactor = 1.0
pctlInput = "/tigress/nwanders/Scripts/Seasonal/resultsNetCDF/"

master = open("/tigress/nwanders/Scripts/hydroSeasonal/jobs/master.sh", "w")

for m in range(len(modelS)):
  model = modelS[m]
  precVarName = precVarNameS[m]
  for r in range(len(refModelS)):
    print r
    refModel = refModelS[r]
    for year in range(beginYear,endYear+1):
      for month in range(1,13):
        day = 1
        stateTime = timeToStr(datetime.datetime(year, month, day) - datetime.timedelta(days=1))
        startTime = datetime.datetime(year, month, day)
        endTime = datetime.datetime(year+1, month, day) - datetime.timedelta(days=1)
        if checkForcingFiles(forcingInput+model+"/", model, precVarName, tempVarName, startTime, endTime):

          outputDir = "/tigress/nwanders/Scripts/hydroSeasonal/"+model+"/PCRGLOBWB/"+refModel+"/"+timeToStr(startTime)
          stateDir = "/tigress/nwanders/Scripts/hydroSeasonal/"+refModel+"/states/"

          con = open("/tigress/nwanders/Scripts/hydroSeasonal/setup_Original.ini")
          iniFile = con.readlines()
          con.close()

          for i in range(len(iniFile)):
            mapEnd = iniFile[i].find('31.map')
            mapPath = iniFile[i].find('/states/')
            startPath = iniFile[i].find('/home/')  
            if iniFile[i][0:9] == 'outputDir':
              iniFile[i] = 'outputDir = '+outputDir+'\n'
            if iniFile[i][0:5] == 'model':
              iniFile[i] = 'model = '+model+'\n'
            if iniFile[i][0:9] == 'startTime':
              iniFile[i] = 'startTime = '+timeToStr(startTime)+'\n'
            if iniFile[i][0:7] == 'endTime':
              iniFile[i] = 'endTime = '+timeToStr(endTime)+'\n'
            if mapEnd > -1:
              iniFile[i] = iniFile[i][0:mapEnd-8] + stateTime + iniFile[i][mapEnd+2:]
            if iniFile[i][0:15] == 'precipitationNC':
              iniFile[i] = 'precipitationNC = '+forcingInput+model+"/"+ncFileName(model, precVarName, startTime, endTime)+'\n'
            if iniFile[i][0:13] == 'temperatureNC':
              iniFile[i] = 'temperatureNC = '+forcingInput+model+"/"+ncFileName(model, tempVarName, startTime, endTime)+'\n'
            if iniFile[i][0:20] == 'precipitationVarName':
              iniFile[i] = 'precipitationVarName = '+precVarName+'\n'
            if iniFile[i][0:18] == 'temperatureVarName':
              iniFile[i] = 'temperatureVarName = '+tempVarName+'\n'
            if iniFile[i][0:30] == 'precipitationCorrectionFactor':
              iniFile[i] = 'precipitationCorrectionFactor = '+precCorFactor+'\n'
            if iniFile[i][0:21] == 'precipitationInputCDF':
              iniFile[i] = 'precipitationInputCDF = '+pctlInput+model+"_"+precVarName+'_pctl.nc4'+'\n'
            if iniFile[i][0:19] == 'temperatureInputCDF':
              iniFile[i] = 'temperatureInputCDF = '+pctlInput+model+"_"+tempVarName+'_pctl.nc4'+'\n'
            if iniFile[i][0:25] == 'precipitationReferenceCDF':
              iniFile[i] = 'precipitationReferenceCDF = '+pctlInput+refModel+"_"+precRefVarName+'_pctl.nc4'+'\n'
            if iniFile[i][0:23] == 'temperatureReferenceCDF':
              iniFile[i] = 'temperatureReferenceCDF = '+pctlInput+refModel+"_"+tempRefVarName+'_pctl.nc4'+'\n'
            if iniFile[i][0:9] == 'nrSamples' and model == "FLOR":
              iniFile[i] = 'nrSamples = 12'+'\n'
            if mapPath > -1:
              iniFile[i] = iniFile[i][0:startPath] + stateDir + iniFile[i][mapPath+8:]


          out = open("/tigress/nwanders/Scripts/hydroSeasonal/config/setup_"+model+"_"+refModel+"_"+timeToStr(startTime)+".ini", "w")
          out.writelines(iniFile)
          out.close()
    
          job = open("/tigress/nwanders/Scripts/hydroSeasonal/jobs/"+model+"_"+refModel+"_"+timeToStr(startTime)+".sh", "w")
          job.writelines("#!/bin/bash\n")
          job.writelines("#SBATCH -n 1   # node count\n")
          job.writelines("#SBATCH -t 23:59:59\n")
          job.writelines("#SBATCH --mail-type=fail\n")
          job.writelines("#SBATCH --mail-user=nwanders@princeton.edu\n")
          job.writelines("cd /tigress/nwanders/Scripts/hydroSeasonal/\n")
          job.writelines("python /tigress/nwanders/Scripts/PCR-GLOBWB/model/seasonal_runner.py /tigress/nwanders/Scripts/hydroSeasonal/config/setup_"+model+"_"+refModel+"_"+timeToStr(startTime)+".ini")
          job.close()
    
          master.writelines("sbatch /tigress/nwanders/Scripts/hydroSeasonal/jobs/"+model+"_"+refModel+"_"+timeToStr(startTime)+".sh\n")

master.close()
