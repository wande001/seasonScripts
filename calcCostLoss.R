require(ncdf4)
require(ncdf)

modelS = c("CanCM3", "CanCM4", "FLOR", "EPS")
varnameS = c("discharge", "storGroundwater", "satDegUpp","satDegLow")
varncS = c("discharge", "groundwater_storage", "upper_soil_saturation_degree","lower_soil_saturation_degree")
refS = c("PGF")
totEns = c(10,10,12,10,9,9)
meanS = matrix(NA, 8,24)
Hit = list()
False = list()
Miss = list()
for(m in 1:length(modelS)){
  Hit[[m]] = array(0, c(720,360,24))
  False[[m]] = array(0, c(720,360,24))
  Miss[[m]] = array(0, c(720,360,24))
}
for(v in 1:length(varnameS)){
  var = varnameS[v]
  for(ref in refS){
    thres = array(NA, c(720, 360, 24))
    months = rep(seq(1,24,1),32)
    fileName = paste(var,"_seasoAvg_output.nc",sep="")
    NC = nc_open(paste("../",ref,"/netcdf/",fileName, sep=""))
    dataRef = ncvar_get(NC, varncS[v])
    nc_close(NC)
    for(x in 1:720){
      print(x)
      for(y in 1:360){
        if(is.na(dataRef[x,y,1]) == FALSE){
          for(month in 1:24){
            thres[x,y,month] = as.numeric(quantile(dataRef[x,y,months==month], probs = 0.1, na.rm=T))
          }
        }
      }
    }
    quantile(thres, probs=seq(0,1,0.1), na.rm=T)
    for(year in 1981:2010){
      for(month in 1:12){
        threshold = thres[,,c((month*2-1):24, 1:(month*2-2))[1:24]]
        if(nchar(month) < 2){
           month = paste("0",as.character(month),sep="")
	      }
        startDay = (as.numeric(year)-1980)*24+(as.numeric(month)-1)*2 + 1
        floodObs = dataRef[,,startDay:(startDay+23)] < threshold
        for(m in 1:length(modelS)){
          nEns = totEns[m]
          model = modelS[m]
          print(paste(model, month, year))
          ensMean = array(0, c(720,360, 24))
          for(n in c(1:nEns)){
            NC = nc_open(paste("../",model,"/PCRGLOBWB/",ref,"/",year,"-",month,"-01/netcdf/",n,"/",fileName, sep=""))
            data = ncvar_get(NC, varncS[v])
            ensMean = ensMean + data/nEns
            nc_close(NC)
          }
          floodMod = ensMean < threshold
          Hit[[m]][which(floodObs & floodMod)] = Hit[[m]][which(floodObs & floodMod)] + 1
          False[[m]][which(floodObs == FALSE & floodMod)] = False[[m]][which(floodObs ==FALSE & floodMod)] + 1
          Miss[[m]][which(floodObs & floodMod == FALSE)] = Miss[[m]][which(floodObs & floodMod == FALSE)] + 1
        }
      }
    }
  }

  HitTemp = Hit
  FalseTemp = False
  MissTemp = Miss

  for(m in 1:4){
    Hit[[m]][is.na(thres)] == NA
    False[[m]][is.na(thres)] == NA
    Miss[[m]][is.na(thres)] == NA
  }

  ratio = list()
  for(m in 1:3){
    ratio[[m]] = (Miss[[4]] - Miss[[m]])/(Hit[[m]] + False[[m]] - Hit[[4]] - False[[4]])
    ratio[[m]][is.finite(ratio[[m]]) == FALSE] = 1.0
    ratio[[m]][ratio[[m]] < 0] = 0.0
  }

  dimLon <- dim.def.ncdf( "Longitude", "-", seq(-179.75,179.75,0.5))
  dimLat <- dim.def.ncdf( "Latitude", "-", seq(89.75,-89.75,-0.5))
  dimTime = dim.def.ncdf( "Time", "-", c(1:24))

  mv <- NA  	# missing value to use
  Lon <- var.def.ncdf( "Longitude", "-", list(dimLon), mv )
  Lat <- var.def.ncdf( "Latitude", "-", list(dimLat), mv)
  Time <- var.def.ncdf( "Time", "-", list(dimTime), mv)
  CanCM3Hits <- var.def.ncdf( "Hits_CanCM3", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM4Hits <- var.def.ncdf( "Hits_CanCM4", "-", list(dimLon, dimLat, dimTime), mv )
  FLORHits <- var.def.ncdf( "Hits_FLOR", "-", list(dimLon, dimLat, dimTime), mv )
  ESPHits <- var.def.ncdf( "Hits_ESP", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM3False <- var.def.ncdf( "False_CanCM3", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM4False <- var.def.ncdf( "False_CanCM4", "-", list(dimLon, dimLat, dimTime), mv )
  FLORFalse <- var.def.ncdf( "False_FLOR", "-", list(dimLon, dimLat, dimTime), mv )
  ESPFalse <- var.def.ncdf( "False_ESP", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM3Miss <- var.def.ncdf( "Miss_CanCM3", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM4Miss <- var.def.ncdf( "Miss_CanCM4", "-", list(dimLon, dimLat, dimTime), mv )
  FLORMiss <- var.def.ncdf( "Miss_FLOR", "-", list(dimLon, dimLat, dimTime), mv )
  ESPMiss <- var.def.ncdf( "Miss_ESP", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM3ratio <- var.def.ncdf( "ratio_CanCM3", "-", list(dimLon, dimLat, dimTime), mv )
  CanCM4ratio <- var.def.ncdf( "ratio_CanCM4", "-", list(dimLon, dimLat, dimTime), mv )
  FLORratio <- var.def.ncdf( "ratio_FLOR", "-", list(dimLon, dimLat, dimTime), mv )


  nc <- create.ncdf(paste("ForecastValueDrought_",var,".nc", sep=""), list(CanCM3Hits, CanCM4Hits, FLORHits, ESPHits, CanCM3False, CanCM4False, FLORFalse, ESPFalse, CanCM3Miss, CanCM4Miss, FLORMiss, ESPMiss, CanCM3ratio, CanCM4ratio, FLORratio))

  put.var.ncdf(nc, CanCM3Hits, Hit[[1]])
  put.var.ncdf(nc, CanCM4Hits, Hit[[2]])
  put.var.ncdf(nc, FLORHits, Hit[[3]])
  put.var.ncdf(nc, ESPHits, Hit[[4]])
  put.var.ncdf(nc, CanCM3False, False[[1]])
  put.var.ncdf(nc, CanCM4False, False[[2]])
  put.var.ncdf(nc, FLORFalse, False[[3]])
  put.var.ncdf(nc, ESPFalse, False[[4]])
  put.var.ncdf(nc, CanCM3Miss, Miss[[1]])
  put.var.ncdf(nc, CanCM4Miss, Miss[[2]])
  put.var.ncdf(nc, FLORMiss, Miss[[3]])
  put.var.ncdf(nc, ESPMiss, Miss[[4]])
  put.var.ncdf(nc, CanCM3ratio, ratio[[1]])
  put.var.ncdf(nc, CanCM4ratio, ratio[[2]])
  put.var.ncdf(nc, FLORratio, ratio[[3]])

  close.ncdf(nc)
}


