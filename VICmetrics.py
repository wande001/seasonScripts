from readData import *
#from plotMatrix import *
from scipy.stats.stats import spearmanr
import sys

lag = int(sys.argv[5])
step = int(sys.argv[2])
end = int(sys.argv[3])
month = int(sys.argv[4])
tempScale = int(sys.argv[1])
model = sys.argv[6]
varName = sys.argv[8]
ref = sys.argv[7]

varOutPutName = varName

if model == "CanCM3":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM3/VIC/"+ref+"/"
    ensNr = 10

if model == "CanCM4":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM4/VIC/"+ref+"/"
    ensNr = 10

if model == "FLOR":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/FLOR/VIC/"+ref+"/"
    ensNr = 12

if model == "EPS":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/EPS/VIC/"+ref+"/"
    ensNr = 10

if model == "CCSM":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CCSM/VIC/"+ref+"/"
    ensNr = 10

if ref == "PGF":
    ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/output_1972_2014.nc"


ncOutputFile = "/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_"+model+"_"+ref+"_"+varOutPutName+"_lag_"+str(lag)+".nc4"

startDays = np.tile(["01","16"],24)
endDays = np.tile(["15","31","15","28","15","31","15","30","15","31","15","30","15","31","15","31","15","30","15","31","15","30","15","31"],2)
inputMonth = np.tile(np.repeat(["01","02","03","04","05","06","07","08","09","10","11","12"],2),2)
inputYear = np.repeat(["2010","2011"],24)
varNames = ["correlation","signif", "bias","NSE", "RMSE","CRPS", "NCRPS"]
varUnits = ["-","-","m3/s","-","m3/s","m3/s","-"]
MV = -999.
createNetCDF(ncOutputFile, varNames, varUnits, np.arange(-89.5,90,1.), np.arange(-179.5,180,1.), loop=True)
posCount = 0

for event in range(0,end,step):
    dateInput = "1981-"+inputMonth[event]+"-"+startDays[event]
    endDay = inputYear[event+month-1]+"-"+inputMonth[event+month-1]+"-"+endDays[event+month-1]
    print dateInput
    print endDay
    
    NMME = returnVICForecast(dateInput, endDay, model, varName, lag, dirLoc = dirLoc, ensNr = ensNr)
    
    print lagToDateStr(dateInput, lag, model)
    print lagToDateStr(endDay, lag, model)
    dataPGF = readForcingVIC(ncRef, varName, dateInput, endDay=endDay, lag=lag, model="PGF")
    
    for space in range(1):
        print space
        spaceNMME = ensembleMean(NMME, ensDimension = 1)
        spacePGF = aggregateSpace(dataPGF, extent=space)
        
        corMap = np.zeros((180,360))-999.
        signMap = np.zeros((180,360))-999.
        biasMap = np.zeros((180,360))-999.
        NSmap = np.zeros((180,360))-999.
        RMSEmap = np.zeros((180,360))-999.
        CRPSmap = np.zeros((180,360))-999.
        NCRPSmap = np.zeros((180,360))-999.

        for i in range(180):
          print i
          for j in range(360):
            if spacePGF[1,i,j] < 1e+10:
              try:
                out = spearmanr(spacePGF[:,i,j], spaceNMME[:,i,j])
              except:
                out = np.ones(2)
                out[0] = -999.
              try:
                nsOut = nashSutcliffe(spacePGF[:,i,j], spaceNMME[:,i,j])
              except:
                nsOut = -999.
              try:
                rmseOut = RMSE(spacePGF[:,i,j], spaceNMME[:,i,j])
              except:
                rmseOut = -999.
              try:
                crpsOut = crps(dataPGF[:,i,j], NMME[:,:,i,j])
              except:
                crpsOut = -999.
              try:
                ncrpsOut = crps(normalize(dataPGF[:,i,j]), normalize(NMME[:,:,i,j]))
              except:
                ncrpsOut = -999.
              corMap[i,j] = out[0]
              signMap[i,j] = out[1]
              biasMap[i,j] = np.mean(spaceNMME[:,i,j]) - np.mean(spacePGF[:,i,j])
              NSmap[i,j] = nsOut
              RMSEmap[i,j] = rmseOut
              CRPSmap[i,j] = crpsOut
              NCRPSmap[i,j] = ncrpsOut
        
        data2NetCDF(ncOutputFile, "correlation", corMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "signif", signMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "bias", biasMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "NSE", NSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "RMSE", RMSEmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "CRPS", CRPSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "NCRPS", NCRPSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
    posCount += 1
    filecache = None

