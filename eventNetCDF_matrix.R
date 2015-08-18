require(ncdf)

modelS = c("CanCM4")
forcingS = c("PGF")
varNameS = c("discharge", "groundwater", "soilMoistureLow","soilMoistureUp")
lim = 0.05
Rlim = 0.0

for(model in modelS){
for(forcing in forcingS){
for(varName in varNameS){

NC = open.ncdf(paste("../resultsNetCDF/",model,"_",forcing,"_",varName,"_tempScale_",1,"_lag_",0,".nc4", sep=""))

Lon = get.var.ncdf(NC, "lon")
Lat = get.var.ncdf(NC, "lat")
T = get.var.ncdf(NC, "time")

close.ncdf(NC)

dimX <- dim.def.ncdf( "lon", "degrees_north", Lon)
dimY <- dim.def.ncdf( "lat", "degrees_east", Lat)
dimT <- dim.def.ncdf( "time", "Days since 1901-01-01", T, unlim=TRUE)

mv <- -9999

PPM = list()
tel = 0
for(i in 0:11){
  for(t in 0:(12-i)){
    tel = tel + 1
    PPM[[tel]] <- var.def.ncdf(paste("PPM",i,t,sep="_"), "-", list(dimX, dimY, dimT), mv,prec="double")
  }
}

nc <- create.ncdf(paste("../resultsNetCDF/",model,"_",forcing,"_",varName,"_PPM_matrix_highLim.nc4",sep=""), PPM)

CCevents = array(0,c(720, 360, 12))
count = 0
varTel = 0
for(lag in 0:11){
  for(temp in 0:(12-lag)){
    varTel = varTel + 1
    NC = open.ncdf(paste("../resultsNetCDF/",model,"_",forcing,"_",varName,"_tempScale_",temp,"_lag_",lag,".nc4", sep=""))
    print(paste("../resultsNetCDF/",model,"_",forcing,"_",varName,"_tempScale_",temp,"_lag_",lag,".nc4", sep=""))
    for(spat in 0:8){
      if(temp ==0){
        for(time in 1:24){
          R = get.var.ncdf(NC, "correlation", start=c(1,1,time), count=c(720,360,1))
          sign = get.var.ncdf(NC, "signif", start=c(1,1,time), count=c(720,360,1))
          out = matrix(0, 720, 360)
          out[sign< lim & R > Rlim] = 1
          count = count + 1
          CCevents[,,ceiling(time/2)] = CCevents[,,ceiling(time/2)] + out
        }
      }
      else{
        for(time in 1:12){
          R = get.var.ncdf(NC, "correlation", start=c(1,1,time), count=c(720,360,1))
          sign = get.var.ncdf(NC, "signif", start=c(1,1,time), count=c(720,360,1))
          out = matrix(0, 720, 360)
          out[sign< lim & R > Rlim] = 1
          count = count + 1
          CCevents[,,time] = CCevents[,,time] + out
        }
      }
    }
    close.ncdf(NC)
    put.var.ncdf(nc, PPM[[varTel]], CCevents/(count/12))
    print(count)
    CCevents = array(0,c(720, 360, 12))
    count = 0
  }
}
rm(CCevents)
rm(R)
rm(sign)
rm(out)
}}}

