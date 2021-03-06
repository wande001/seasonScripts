from readData import *
import sys
import pandas as pd

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/hydrobasins.nc"
f = nc.Dataset(ncFile)
catchments = f.variables["basin"][:,:]
catchFlat = np.array(catchments).reshape(360*720)
f.close()

ncRef3 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/discharge_seasoAvg_output.nc"
ncRef1 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/referencePotET_seasoAvg_output.nc"
ncRef2 = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/totalEvaporation_seasoAvg_output.nc"

ncOutputFile = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/basinLosses_seasoAvg_output.nc"

startDays = np.tile(["01","16"],24)
endDays = np.tile(["15","31","15","28","15","31","15","30","15","31","15","30","15","31","15","31","15","30","15","31","15","30","15","31"],2)
inputMonth = np.tile(np.repeat(["01","02","03","04","05","06","07","08","09","10","11","12"],2),2)
inputYear = np.repeat(["2010","2011"],24)
varNames = ["hydroPower", "hydroPowerEnvi", "hydroPowerReduction", "irrigationReduction", "irrigationDemand"]
varUnits = ["W","W", "W", "m3/s", "m3/s"]
MV = 1.e+20
createNetCDF(ncOutputFile, varNames, varUnits, np.arange(89.75,-90,-0.5), np.arange(-179.75,180,0.5), loop=True)
posCount = 0

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/routing/cellarea30min.nc"
f = nc.Dataset(ncFile)
area = f.variables["area"][::-1,:]
f.close()
ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/efficiency/efficiency.nc"
f = nc.Dataset(ncFile)
efficiency = f.variables["efficiency"][::-1,:]
f.close()
ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/routing/downstream_elevation.nc"
f = nc.Dataset(ncFile)
potHeight = f.variables["elev"][::-1,:]
f.close()

f1 = nc.Dataset(ncRef1)
f2 = nc.Dataset(ncRef2)
f3 = nc.Dataset(ncRef3)

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/irrigated_areas/irrigationArea30ArcMin.nc"
f = nc.Dataset(ncFile)
irriTime = f.variables["time"][:]

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/domesticWaterDemand30ArcMin.nc"
domestic = nc.Dataset(ncFile)
domTime = domestic.variables["time"][:]

ncFile = "/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/landSurface/waterDemand/industryWaterDemand30ArcMin.nc"
industry = nc.Dataset(ncFile)
indTime = industry.variables["time"][:]

ncFile = "/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/discharge_envirFlow_output.nc"
environment = nc.Dataset(ncFile)

for d in range(816):
  print d
  day = datetime.datetime(1901,1,1)+datetime.timedelta(days = int(f1.variables["time"][d]))
  potPGF = f1.variables["reference_potential_evaporation"][d,:,:]
  actPGF = f2.variables["total_evaporation"][d,:,:]
  dischargePGF = f3.variables["discharge"][d,:,:] * 86400.
  powerPGF = f3.variables["discharge"][d,:,:] * potHeight * 1000. * 9.81
  
  irriSel = np.where(irriTime < int(f1.variables["time"][d]))[0][-1]
  irriArea = f.variables["irrigationArea"][irriSel,:,:] * 100. * 100.
  domSel = np.where(domTime < int(f1.variables["time"][d]))[0][-1]
  domDemand = domestic.variables["domesticGrossDemand"][domSel,:,:]
  indSel = np.where(indTime < int(f1.variables["time"][d]))[0][-1]
  indDemand = industry.variables["industryGrossDemand"][indSel,:,:]
  
  envDemand = environment.variables["discharge"][day.month-1,:,:] * 86400.
  powerEnvi = environment.variables["discharge"][day.month-1,:,:] * potHeight * 1000. * 9.81
  
  Demand = (potPGF-actPGF)*irriArea/efficiency + (domDemand + indDemand) * area
  
  outMap = np.zeros((360,720))
  envMap = np.zeros((360,720))
  hydroMap = np.zeros((360,720))
  hydroRedMap = np.zeros((360,720))
  hydroEnviMap = np.zeros((360,720))
  catchDemand = np.zeros((np.max(catchments)+1))
  catchDischarge = np.zeros((np.max(catchments)+1))
  catchEnvironment = np.zeros((np.max(catchments)+1))
  hydroPower = np.zeros((np.max(catchments)+1))
  hydroEnvi = np.zeros((np.max(catchments)+1))
  df = pd.DataFrame({ 'catchment' : catchFlat.reshape((360*720)), 'discharge' : dischargePGF.reshape(360*720), 'Power' : powerPGF.reshape(360*720), 'PowerEnvi' : powerEnvi.reshape(360*720), 'Demand' : Demand.reshape(360*720), 'environment' : dischargePGF.reshape(360*720)-envDemand.reshape(360*720)})
  result = df.groupby(['catchment']).sum()
  locs = result.index.values[1:(np.max(catchments)+1)]
  catchDemand[locs] = np.array(result["Demand"])[1:np.max(catchments)+2]
  hydroPower[locs] = np.array(result["Power"])[1:np.max(catchments)+2]
  hydroEnvi[locs] = np.array(result["PowerEnvi"])[1:np.max(catchments)+2]
  result = df.groupby(['catchment']).max()
  catchDischarge[locs] = np.array(result["discharge"])[1:np.max(catchments)+2]
  catchEnvironment[locs] = np.array(result["environment"])[1:np.max(catchments)+2]
  
  for c in locs:
    if catchDischarge[c] > 0.:
      resultsMask = catchments==c
      try:
        outMap[resultsMask] = catchDemand[c]
      except:
        outMap[resultsMask] = 0.0
      try:
        envMap[resultsMask] = np.maximum(catchDemand[c] - np.minimum(catchDemand[c],catchEnvironment[c]), 0.0)
      except:
        envMap[resultsMask] = 0.0
      try:
        hydroRedMap[resultsMask] = (1.0 - np.maximum(np.minimum(catchDemand[c]/catchEnvironment[c], 1.0), 0.0)) * (hydroPower[c]-hydroEnvi[c]) + hydroEnvi[c]
      except:
        hydroRedMap[resultsMask] = 0.0
      try:
        hydroMap[resultsMask] = hydroPower[c]
        hydroEnviMap[resultMask] = hydroEnvi[c]
      except:
        hydroMap[resultsMask] = 0.0
        hydroEnviMap[resultsMask] = 0.0
  resultsMask = catchments==0
  outMap[resultsMask] = MV
  envMap[resultsMask] = MV
  hydroMap[resultsMask] = MV
  hydroRedMap[resultsMask] = MV
  hydroEnviMap[resultsMask] = MV
  data2NetCDF(ncOutputFile, "irrigationDemand", outMap, day,posCnt = posCount)
  data2NetCDF(ncOutputFile, "irrigationReduction", envMap, day,posCnt = posCount)
  data2NetCDF(ncOutputFile, "hydroPower", hydroMap, day,posCnt = posCount)
  data2NetCDF(ncOutputFile, "hydroPowerEnvi", hydroEnviMap, day,posCnt = posCount)
  data2NetCDF(ncOutputFile, "hydroPowerReduction", hydroRedMap, day,posCnt = posCount)
  posCount += 1
  filecache = None

industry.close()
domestic.close()
environment.close()
f.close()
f1.close()
f2.close()
f3.close()
