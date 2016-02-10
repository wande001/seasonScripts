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

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/majorCatchments.nc"
f = nc.Dataset(ncFile)
catchments = f.variables["Band1"][::-1,:]
f.close()
multiCatch = np.zeros((30,360,720))
for i in range(30):
  multiCatch[i,:,:] = catchments


if varName == "discharge":
  varFile = "discharge_dailyTot_output.nc"
  varOutPutName = "discharge"

if varName == "groundwater":
  varFile = "storGroundwater"
  varOutPutName = "groundwater"

if varName == "soilMoistureLow":
  varFile = "satDegLow"
  varOutPutName = "soilMoistureLow"

if varName == "soilMoistureUp":
  varFile = "satDegUpp"
  varOutPutName = "soilMoistureUp"

if varName == "evap":
  varName_1 = "reference_potential_evaporation"
  varFile_1 = "referencePotET"
  varName_2 = "total_evaporation"
  varFile_2 = "totalEvaporation"
  varOutPutName = "evapDeficit"


if model == "CanCM3":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM3/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "CanCM4":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM4/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "FLOR":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/FLOR/PCRGLOBWB/"+ref+"/"
    ensNr = 12

if model == "EPS":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/ESP/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "CCSM":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CCSM/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "Weighted":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/Weighted/PCRGLOBWB/"+ref+"/"
    ensNr = 9

if model == "WeightedEqual":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/WeightedEqual/PCRGLOBWB/"+ref+"/"
    ensNr = 9

if ref == "PGF":
    if varName == "discharge":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/discharge_seasoAvg_output.nc"
    if varName == "groundwater":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/storGroundwater_seasoAvg_output.nc"
    if varName == "soilMoistureLow":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/satDegLow_seasoAvg_output.nc"
    if varName == "soilMoistureUp":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/satDegUpp_seasoAvg_output.nc"
    if varName == "evap":
        ncRef1 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/referencePotET_seasoAvg_output.nc"
        ncRef2 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/totalEvaporation_seasoAvg_output.nc"

if ref == "CFS":
    if varName == "discharge":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/CFS/netcdf/discharge_seasoAvg_output.nc"
    if varName == "groundwater":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/CFS/netcdf/storGroundwater_seasoAvg_output.nc"
    if varName == "soilMoistureLow":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/CFS/netcdf/satDegLow_seasoAvg_output.nc"
    if varName == "soilMoistureUp":
        ncRef = "/tigress/nwanders/Scripts/hydroSeasonal/CFS/netcdf/satDegUpp_seasoAvg_output.nc"


ncOutputFile = "/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_"+model+"_"+ref+"_"+varOutPutName+"_lag_"+str(lag)+"Catchment.nc4"

startDays = np.tile(["01","16"],24)
endDays = np.tile(["15","31","15","28","15","31","15","30","15","31","15","30","15","31","15","31","15","30","15","31","15","30","15","31"],2)
inputMonth = np.tile(np.repeat(["01","02","03","04","05","06","07","08","09","10","11","12"],2),2)
inputYear = np.repeat(["2010","2011"],24)
varNames = ["correlation","signif", "bias","NSE", "RMSE", "CRPS", "NCRPS"]
varUnits = ["-","-","m3/s","-","m3/s", "m3/s", "-"]
MV = -999.
createNetCDF(ncOutputFile, varNames, varUnits, np.arange(89.75,-90,-0.5), np.arange(-179.75,180,0.5), loop=True)
posCount = 0

for event in range(0,end,step):
    dateInput = "1981-"+inputMonth[event]+"-"+startDays[event]
    endDay = inputYear[event+month-1]+"-"+inputMonth[event+month-1]+"-"+endDays[event+month-1]
    print dateInput
    print endDay
    
    if varName != "evap":
      NMME = returnSeasonalForecast(dateInput, endDay, model, varName, lag, dirLoc = dirLoc, ensNr = ensNr)
    else:
      potEvap = returnSeasonalForecast(dateInput, endDay, model, varName_1, lag, dirLoc = dirLoc, ensNr = ensNr)
      evap = returnSeasonalForecast(dateInput, endDay, model, varName_2, lag, dirLoc = dirLoc, ensNr = ensNr)
      NMME = potEvap - evap
    print lagToDateStr(dateInput, lag, model)
    print lagToDateStr(endDay, lag, model)
    if varName != "evap":
      dataPGF = readForcing(ncRef1, varName, dateInput, endDay=endDay, lag=lag, model="PGF")
    else:
      potEvapPGF = readForcing(ncRef1, varName_1, dateInput, endDay=endDay, lag=lag, model="PGF")
      evapPGF = readForcing(ncRef2, varName_2, dateInput, endDay=endDay, lag=lag, model="PGF")
      dataPGF = potEvapPGF - evapPGF
    
    for space in range(1):
        print space
        ensNMME = ensembleMean(NMME, ensDimension = 1)
        spaceNMME = ensNMME - np.mean(ensNMME, axis=0)
        spacePGF = dataPGF - np.mean(dataPGF, axis=0)
        
        corMap = np.zeros((360,720))-999.
        signMap = np.zeros((360,720))-999.
        biasMap = np.zeros((360,720))-999.
        NSmap = np.zeros((360,720))-999.
        RMSEmap = np.zeros((360,720))-999.
        CRPSmap = np.zeros((360,720))-999.
        NCRPSmap = np.zeros((360,720))-999.

        for c in range(1,(np.max(catchments)+1)):
          print c
          if len(np.where(catchments == c)[0]) > 10.:
            maskNMME = np.ma.masked_where(multiCatch != c, spaceNMME).mean(axis=1).mean(axis=1)
            maskPGF = np.ma.masked_where(multiCatch != c, spacePGF).mean(axis=1).mean(axis=1)
            try:
              out = spearmanr(maskPGF, maskNMME)
            except:
              out = np.ones(2)
              out[0] = -999.
            try:
              nsOut = nashSutcliffe(maskPGF, maskNMME)
            except:
              nsOut = -999.
            try:
              rmseOut = RMSE(maskPGF, maskNMME)
            except:
              rmseOut = -999.
            try:
              crpsOut = crps(maskPGF, maskNMME)
            except:
              crpsOut = -999.
            try:
              ncrpsOut = crps(normalize(maskPGF, maskNMME))
            except:
              ncrpsOut = -999.
            corMap[catchments==c] = out[0]
            signMap[catchments==c] = out[1]
            biasMap[catchments==c] = np.mean(maskPGF -  maskNMME)
            NSmap[catchments==c] = nsOut
            RMSEmap[catchments==c] = rmseOut
            CRPSmap[catchments==c] = crpsOut
            NCRPSmap[catchments==c] = ncrpsOut
        
        data2NetCDF(ncOutputFile, "correlation", corMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "signif", signMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "bias", biasMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "NSE", NSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "RMSE", RMSEmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "CRPS", CRPSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile, "NCRPS", NCRPSmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
    posCount += 1
    filecache = None

