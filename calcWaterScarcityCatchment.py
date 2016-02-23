from readData import *
#from plotMatrix import *
from scipy.stats.stats import spearmanr
import sys
import pandas as pd


#lag = int(sys.argv[5])
#step = int(sys.argv[2])
#end = int(sys.argv[3])
#month = int(sys.argv[4])
#tempScale = int(sys.argv[1])
#model = sys.argv[6]
#varName = sys.argv[8]
#ref = sys.argv[7]

lag = int(0)
step = int(1)
end = int(1)
month = int(1)
tempScale = int(0)
model = "ESP"
varName = "scarcity"
ref = "PGF"

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/majorCatchments.nc"
f = nc.Dataset(ncFile)
catchments = f.variables["Band1"][::-1,:]
f.close()
catchFlat = np.array(catchments).reshape(360*720)
multiCatch = np.zeros((30,360,720))
for i in range(30):
  multiCatch[i,:,:] = catchments

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/hydrobasins.nc"
f = nc.Dataset(ncFile)
catchments = f.variables["basin"][:,:]
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

if varName == "scarcity":
  varName_1 = "reference_potential_evaporation"
  varFile_1 = "referencePotET"
  varName_2 = "total_evaporation"
  varFile_2 = "totalEvaporation"
  varName_3 = "discharge"
  varFile_3 = "discharge"
  varOutPutName = "waterScarcity"
  ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/routing/cellarea30min.nc"
  f = nc.Dataset(ncFile)
  area = f.variables["area"][::-1,:]
  f.close()
  multiarea = np.zeros((30,360,720))
  for i in range(30):
    multiarea[i,:,:] = area
  ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/efficiency/efficiency.nc"
  f = nc.Dataset(ncFile)
  efficiency = f.variables["efficiency"][::-1,:]
  f.close()
  multiEfficiency = np.zeros((30,360,720))
  for i in range(30):
    multiEfficiency[i,:,:] = efficiency


if model == "CanCM3":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM3/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "CanCM4":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/CanCM4/PCRGLOBWB/"+ref+"/"
    ensNr = 10

if model == "FLOR":
    dirLoc = "/tigress/nwanders/Scripts/hydroSeasonal/FLOR/PCRGLOBWB/"+ref+"/"
    ensNr = 12

if model == "ESP":
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
    if varName == "scarcity":
        ncRef3 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/discharge_seasoAvg_output.nc"
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


ncOutputFile1 = "/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_"+model+"_"+ref+"_waterDemand_lag_"+str(lag)+"Catchment.nc4"
ncOutputFile2 = "/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_"+model+"_"+ref+"_waterScarcity_lag_"+str(lag)+"Catchment.nc4"

startDays = np.tile(["01","16"],24)
endDays = np.tile(["15","31","15","28","15","31","15","30","15","31","15","30","15","31","15","31","15","30","15","31","15","30","15","31"],2)
inputMonth = np.tile(np.repeat(["01","02","03","04","05","06","07","08","09","10","11","12"],2),2)
inputYear = np.repeat(["2010","2011"],24)
varNames = ["correlation","signif", "bias","RMSE"]
varUnits = ["-","-","m3/s","m3/s"]
MV = 1.e+20
createNetCDF(ncOutputFile1, varNames, varUnits, np.arange(89.75,-90,-0.5), np.arange(-179.75,180,0.5), loop=True)
createNetCDF(ncOutputFile2, varNames, varUnits, np.arange(89.75,-90,-0.5), np.arange(-179.75,180,0.5), loop=True)
posCount = 0

for event in range(0,end,step):
    dateInput = "1981-"+inputMonth[event]+"-"+startDays[event]
    endDay = inputYear[event+month-1]+"-"+inputMonth[event+month-1]+"-"+endDays[event+month-1]
    print dateInput
    print endDay    
    if varName != "evap" and varName != "scarcity":
      NMME = returnSeasonalForecast(dateInput, endDay, model, varName, lag, dirLoc = dirLoc, ensNr = ensNr)
    elif varName == "evap":
      potEvap = returnSeasonalForecast(dateInput, endDay, model, varName_1, lag, dirLoc = dirLoc, ensNr = ensNr)
      evap = returnSeasonalForecast(dateInput, endDay, model, varName_2, lag, dirLoc = dirLoc, ensNr = ensNr)
      NMME = potEvap - evap
    else:
      dateInputStatic = "1981-"+inputMonth[event]+"-01"
      domestic = readForcingstatistic("/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/domesticWaterDemand30ArcMin.nc", "domesticGrossDemand", dateInputStatic, endDay=endDay, lag=lag, model="PGF")
      industry = np.zeros((30,360,720))
      ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/industryWaterDemand30ArcMin.nc"
      f = nc.Dataset(ncFile)
      for i in range(21,51):
        industry[i-21,:,:] = f.variables["industryGrossDemand"][i,:,:]
      f.close()
      irriArea = np.zeros((30,360,720))
      ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/irrigated_areas/irrigationArea30ArcMin.nc"
      f = nc.Dataset(ncFile)
      for i in range(21,51):
        irriArea[i-21,:,:] = f.variables["irrigationArea"][i,:,:]
      f.close()
      potEvap = returnSeasonalForecast(dateInput, endDay, model, varName_1, lag, dirLoc = dirLoc, ensNr = ensNr) 
      evap = returnSeasonalForecast(dateInput, endDay, model, varName_2, lag, dirLoc = dirLoc, ensNr = ensNr)
      discharge = ensembleMean(returnSeasonalForecast(dateInput, endDay, model, varName_3, lag, dirLoc = dirLoc, ensNr = ensNr), ensDimension = 1)*86400.
      totalDemand = ensembleMean(potEvap-evap, ensDimension = 1)*irriArea*multiEfficiency + (domestic + industry)*multiarea
      del(potEvap)
      del(evap)
    print lagToDateStr(dateInput, lag, model)
    print lagToDateStr(endDay, lag, model)
    if varName != "evap" and varName != "scarcity":
      dataPGF = readForcing(ncRef1, varName, dateInput, endDay=endDay, lag=lag, model="PGF")
    elif varName == "evap":
      potEvapPGF = readForcing(ncRef1, varName_1, dateInput, endDay=endDay, lag=lag, model="PGF")
      evapPGF = readForcing(ncRef2, varName_2, dateInput, endDay=endDay, lag=lag, model="PGF")
    else:
      potEvapPGF = readForcing(ncRef1, varName_1, dateInput, endDay=endDay, lag=lag, model="PGF")
      evapPGF = readForcing(ncRef2, varName_2, dateInput, endDay=endDay, lag=lag, model="PGF")
      dischargePGF = readForcing(ncRef3, varName_3, dateInput, endDay=endDay, lag=lag, model="PGF")*86400.
      totalDemandPGF = (potEvapPGF-evapPGF)*irriArea*multiEfficiency + (domestic + industry)*multiarea
      del(potEvapPGF)
      del(evapPGF)
    for space in range(1):
        print space
        spaceNMME = totalDemand
        spacePGF = totalDemandPGF
        corMap = np.zeros((360,720))-np.nan
        signMap = np.zeros((360,720))-np.nan
        biasMap = np.zeros((360,720))-np.nan
        RMSEmap = np.zeros((360,720))-np.nan
        corMapScar = np.zeros((360,720))-np.nan
        signMapScar = np.zeros((360,720))-np.nan
        biasMapScar = np.zeros((360,720))-np.nan
        RMSEmapScar = np.zeros((360,720))-np.nan
        annualCatchDemand = np.zeros((30,np.max(catchments)+1))
        annualCatchDemandPGF = np.zeros((30,np.max(catchments)+1))
        annualCatchScarDemand = np.zeros((30,np.max(catchments)+1))
        annualCatchScarPGF = np.zeros((30,np.max(catchments)+1))
        for year in range(30):
	  df = pd.DataFrame({ 'catchment' : catchFlat.reshape((360*720)), 'demand' : totalDemand[year,:,:].reshape(360*720), 'demandPGF' : totalDemandPGF[year,:,:].reshape(360*720)})
          result = df.groupby(['catchment']).sum()
          locs = result.index.values[1:(np.max(catchments)+1)]
          annualCatchDemand[year,locs] = np.array(result["demand"])[1:np.max(catchments)+2]
          annualCatchDemandPGF[year,locs] = np.array(result["demandPGF"])[1:np.max(catchments)+2]
	  df = pd.DataFrame({ 'catchment' : catchFlat.reshape((360*720)), 'discharge' : discharge[year,:,:].reshape(360*720), 'dischargePGF' : dischargePGF[year,:,:].reshape(360*720)})
          result = df.groupby(['catchment']).max()
          annualCatchScarDemand[year,locs] = np.array(result["discharge"])[1:np.max(catchments)+2]
          annualCatchScarPGF[year,locs] = np.array(result["dischargePGF"])[1:np.max(catchments)+2]
        for c in locs:
          print c
          try:
            out = spearmanr(annualCatchDemandPGF[:,c], annualCatchDemand[:,c])
            outScar = spearmanr(annualCatchDemandPGF[:,c]/annualCatchScarPGF[:,c], annualCatchDemand[:,c]/annualCatchScarDemand[:,c])
          except:
            out = np.ones(2)
            out[0] = np.nan
            outScar = np.ones(2)
            outScar[0] = np.nan
          try:
            rmseOut = RMSE(annualCatchDemandPGF[:,c], annualCatchDemand[:,c])
            rmseOutScar = RMSE(annualCatchDemandPGF[:,c]/annualCatchScarPGF[:,c], annualCatchDemand[:,c]/annualCatchScarDemand[:,c])
          except:
            rmseOut = np.nan
            rmseOutScar = np.nan
          resultsMask = catchments==c
          corMap[resultsMask] = out[0]
          signMap[resultsMask] = out[1]
          biasMap[resultsMask] = np.mean(annualCatchDemandPGF[:,c] - annualCatchDemand[:,c])
          RMSEmap[resultsMask] = rmseOut        
          corMapScar[resultsMask] = outScar[0]
          signMapScar[resultsMask] = outScar[1]
          biasMapScar[resultsMask] = np.mean(annualCatchDemandPGF[:,c]/annualCatchScarPGF[:,c] - annualCatchDemand[:,c]/annualCatchScarDemand[:,c])
          RMSEmapScar[resultsMask] = rmseOutScar        
        resultsMask = catchments==0
        corMap[resultsMask] = MV
        signMap[resultsMask] = MV
        biasMap[resultsMask] = MV
        RMSEmap[resultsMask] = MV
        corMapScar[resultsMask] = MV
        signMapScar[resultsMask] = MV
        biasMapScar[resultsMask] = MV
        RMSEmapScar[resultsMask] = MV
        data2NetCDF(ncOutputFile1, "correlation", corMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile1, "signif", signMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile1, "bias", biasMap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile1, "RMSE", RMSEmap, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile2, "correlation", corMapScar, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile2, "signif", signMapScar, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile2, "bias", biasMapScar, lagToDateTime(dateInput, 0, model), posCnt = posCount)
        data2NetCDF(ncOutputFile2, "RMSE", RMSEmapScar, lagToDateTime(dateInput, 0, model), posCnt = posCount)        
    posCount += 1
    filecache = None
