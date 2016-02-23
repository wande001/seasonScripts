require(rgdal)
require(ncdf4)
require(fields)
require(akima)

makeMatrix <- function(rows, cols, size=5, legendSize = 2, labelSize = 1){
  A = matrix(rows*cols+rows+1, rows*size, cols*size+legendSize+labelSize)
  A[,1] = rep(c(1:rows), each=size)
  for(i in c(1:(rows*cols))){
    row = floor((i-1)/cols)
    col = (i-1) - row * cols
    A[(row*size+1):(row*size+size),(col*size+1+labelSize):(col*size+size+labelSize)] = i + rows
  }
  return(A)
}

makeMatrixDis <- function(rows, cols, size=5, legendSize = 2, labelSize = 1){
  A = matrix(rows*cols+rows+2, (rows+1)*size, cols*size+legendSize+labelSize)
  A[,1] = c(rep(c(1:rows), each=size), rep(rows+rows*cols+1, size))
  for(i in c(1:(rows*cols))){
    row = floor((i-1)/cols)
    col = (i-1) - row * cols
    A[(row*size+1):(row*size+size),(col*size+1+labelSize):(col*size+size+labelSize)] = i + rows
  }
  row = floor((i)/cols)
  A[(row*size+1):(row*size+size),(1+labelSize):(cols*size+labelSize)] = i + rows + 1
  return(A)
}

makeMatrixDisLarge <- function(rows, cols, size=5, legendSize = 2, labelSize = 1, extra=2){
  A = matrix(rows*cols+rows+extra+1, (rows+1)*size, cols*size+legendSize+labelSize)
  A[,1] = c(rep(c(1:rows), each=size), rep(rows+rows*cols+1, size))
  for(i in c(1:(rows*cols))){
    row = floor((i-1)/cols)
    col = (i-1) - row * cols
    A[(row*size+1):(row*size+size),(col*size+1+labelSize):(col*size+size+labelSize)] = i + rows
  }
  row = floor((i)/cols)
  A[(row*size+1):(row*size+size),(1+labelSize):(cols/2*size+labelSize)] = i + rows + 1
  A[(row*size+1):(row*size+size),(cols/2*size+labelSize):(cols*size+labelSize)] = i + rows + 2
  return(A)
}



modelS = c("CanCM3", "CanCM4", "FLOR","CCSM", "Weighted", "WeightedEqual", "EPS", "Weigthed", "ACF")
modelNames = c("CanCM3", "CanCM4", "FLOR","CCSM", "Weighted", "WeightedEqual", "ESP", "Weigthed", "ACF")
varnameRef = c("discharge", "storGroundwater", "satDegUpp","satDegLow")
varnameS = c("discharge", "groundwater", "soilMoistureUp","soilMoistureLow")
varncS = c("discharge", "groundwater_storage", "upper_soil_saturation_degree","lower_soil_saturation_degree")
refS = c("PGF")

area = matrix(readGDAL("/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/routing/catchmentarea_30min.map")@data$band1, 720, 360)
NC = nc_open("/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/majorCatchments.nc")
catchments = ncvar_get(NC, "Band1")
nc_close(NC)

areaSel = list()
classes = 20
for(i in 1:classes){
  areaSel[[i]] = which(log(area) >= quantile(log(area), probs=((i-1)/classes), na.rm=T) & log(area) <= quantile(log(area), probs=(i/classes), na.rm=T))
}

interpolation = list()
latInter = list()
forecastMonth = c(1:24,1:24)


pdf("interpolatedCorrelations.pdf", width=7, height=6)

for(v in 1){
  var = varnameS[v]
  ref= refS
  for(m in c(1:3,7)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:classes){
        tel = tel + 1
        vals[tel] = median(dataCOR[areaSel[[i]]], na.rm=T)
        leads[tel] = lag
        sizes[tel] = median(log(area[areaSel[[i]]]), na.rm=T)
      }
      for(l in 1:360){
        latTel = latTel + 1
        latVals[latTel] = median(dataCOR[,l,], na.rm=T)
        lats[latTel] = 90.25 - l/2
        latLeads[latTel] = lag
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  sel = which(is.na(latVals) == FALSE)
  latInter[[m]] <- interp(latLeads[sel], lats[sel], latVals[sel], yo=seq(min(lats[sel]), max(lats[sel]), length=40))
  }


layout(makeMatrix(2,2))
par(mar=c(0,0,0,0))
for(var in c(1:2)){
  par(mar=c(2.0,0.2,2.0,0))
  if(var == 2){
    par(mar=c(2,0.2,2.0,0)) 
  }
  plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
  text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
  text(0.22, 6.0, "Lead time", srt=90)
}
for(m in c(1:2)){
  par(mar=c(2.0,0.2,2.0,0))
  toPlot = interpolation[[m]]
  toPlot$z[toPlot$z < 0] = 0.0
  image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
  contour(toPlot, add=TRUE)
}
for(m in c(3,7)){
  par(mar=c(2,0.2,2.0,0)) 
  toPlot = interpolation[[m]]
  toPlot$z[toPlot$z < 0] = 0.0
  image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
  axis(1, at = seq(19,26,1), labels=signif(exp(seq(19,26,1))/(1000*1000), digits=3))
  contour(toPlot, add=TRUE)
}
par(mar=c(1,2.7,1,3.5))
colLen = 100
cols = tim.colors
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))



layout(makeMatrix(2,2))
par(mar=c(0,0,0,0))
for(var in c(1:2)){
  par(mar=c(2.0,0.2,2.0,0))
  if(var == 2){
    par(mar=c(2,0.2,2.0,0)) 
  }
  plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
  text(0.75, seq(0.1,11.75, length=15), seq(-55,85,10))
  text(0.22, 6.0, "Lead time", srt=90)
}
for(m in c(1:2)){
  par(mar=c(2.0,0.2,2.0,0))
  toPlot = latInter[[m]]
  toPlot$z[toPlot$z < 0] = 0.0
  image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
  contour(toPlot, add=TRUE)
}
for(m in c(3,7)){
  par(mar=c(2,0.2,2.0,0)) 
  toPlot = latInter[[m]]
  toPlot$z[toPlot$z < 0] = 0.0
  image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
  axis(1, at = seq(0,11,1), labels=seq(0,11,1))
  contour(toPlot, add=TRUE)
}
par(mar=c(1,2.7,1,3.5))
colLen = 100
cols = tim.colors
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
dev.off()

pdf("timeCorrelationsNH.pdf", width=7, height=6)

for(v in 1:4){
  var = varnameS[v]
  ref= refS
  for(m in c(1:3,7)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,1:180,i], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrix(2,2))
  par(mar=c(0,0,0,0))
  for(var in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(var == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    contour(toPlot, add=TRUE)
  }
  for(m in c(3,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE)
  }
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
dev.off()

pdf("timeCorrelationsSH.pdf", width=7, height=6)

for(v in 1:4){
  var = varnameS[v]
  ref= refS
  for(m in c(1:3,7)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,181:360,i], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrix(2,2))
  par(mar=c(0,0,0,0))
  for(var in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(var == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    contour(toPlot, add=TRUE)
  }
  for(m in c(3,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE)
  }
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
dev.off()

pdf("timeCorrelationsTropics.pdf", width=7, height=6)

for(v in 1:4){
  var = varnameS[v]
  ref= refS
  for(m in c(1:3,7)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,140:220,i], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrix(2,2))
  par(mar=c(0,0,0,0))
  for(var in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(var == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    contour(toPlot, add=TRUE)
  }
  for(m in c(3,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE)
  }
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
dev.off()



pdf("timeCorrelationsVICNH.pdf", width=7, height=6)

for(v in 1){
  var = varnameS[v]
  ref= refS
  for(m in c(1:3,7)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,1:180,i], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrix(2,2))
  par(mar=c(0,0,0,0))
  for(var in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(var == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    contour(toPlot, add=TRUE)
  }
  for(m in c(3,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE)
  }
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
dev.off()


pdf("timeCorrelationsMajorRivers.pdf", width=10, height=7)
catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray", "Colorado")
catchTel = 0
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828, 7047)){
catchTel = catchTel + 1
sel = which(catchments == catch)

for(v in 1){
  var = varnameS[v]
  ref= refS
  for(m in c(1:4,6:8,9)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      if(m < 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      }
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",ref,"_",var,"_lag_",lag,"_Coef.nc4", sep=""))
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      }
      if(m == 9){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/output_1981_2012_ACF_lag_",lag,".nc", sep=""))
        dataCOR = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,24))
      }
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
        leads[tel] = lag+((i+1)/2)-floor(i/2)
        sizes[tel] = ceiling(i/2)
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrixDis(2,4))
  par(mar=c(0,0,0,0))
  for(i in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(i == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:3,4)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  for(m in c(6,8,7,9)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/netcdf/",var,"_seasoAvg_annualCycle.nc", sep=""))
  data = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,24))
  nc_close(NC)
  if(var == "discharge"){
    selOutlet = which(catchments == catch)[which.max(data[,360:1,1][which(catchments == catch)])]
    multiply = 1
  }
  else{
    selOutlet = which(catchments == catch)
    multiply = area[,360:1][selOutlet]
  }
  season = array()
  for(i in 1:24){
    season[i] = sum(data[,360:1,i][selOutlet]*multiply)
  }
  par(mar=c(4,4,2.2,0)) 
  plot(1,1, type="n", xlim=c(0,13),ylim=range(season), axes=FALSE, xaxs="i", xlab="Month", ylab=var, main=catchNames[catchTel])
  lines(seq(0,12,0.5), c(season, season[1]))
  axis(1, at=seq(0.5,11.5,1), labels=c(1:12))
  axis(2)  
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
}
dev.off()

pdf("timeCorrelationsMajorRiversVIC.pdf", width=8, height=7)
catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray", "Colorado")
catchTel = 0
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828, 7047)){
catchTel = catchTel + 1
sel = which(catchments == catch)

for(v in 1){
  var = varnameS[v]
  ref= refS
  for(m in c(1:4,6,8,7,9)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      if(m < 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      }
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_",ref,"_",var,"_lag_",lag,"_Coef.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      }
      if(m == 9){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/output_1981_2012_ACF_lag_",lag,".nc", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,time))
        dataCOR[dataCOR < -1] = NA
        dataCOR[dataCOR > 1] = NA        
      }
      nc_close(NC)
      for(i in 1:time){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrixDis(2,4))
  par(mar=c(0,0,0,0))
  for(i in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(i == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:4)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  for(m in c(6,8,7,9)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/VIC_",var,"_seasoAvg_annualCycle.nc", sep=""))
  data = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,24))
  nc_close(NC)
  if(var == "discharge"){
    selOutlet = which(catchments == catch)[which.max(data[,360:1,1][which(catchments == catch)])]
    multiply = 1
  }
  else{
    selOutlet = which(catchments == catch)
    multiply = area[,360:1][selOutlet]
  }
  season = array()
  for(i in 1:24){
    season[i] = sum(data[,360:1,i][selOutlet]*multiply)
  }
  par(mar=c(4,4,2.2,0)) 
  plot(1,1, type="n", xlim=c(0,13),ylim=range(season), axes=FALSE, xaxs="i", xlab="Month", ylab=var, main=catchNames[catchTel])
  lines(seq(0,12,0.5), c(season, season[1]))
  axis(1, at=seq(0.5,11.5,1), labels=c(1:12))
  axis(2)  
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
}
dev.off()


pdf("timeCorrelationsIrriMajorRivers.pdf", width=7, height=6)
catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray")
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828)){
sel = which(catchments == catch)

for(v in 1){
  var = "evapDeficit"
  ref= refS
  for(m in c(1:3,6:8)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrixDis(2,3))
  par(mar=c(0,0,0,0))
  for(var in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(var == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:3)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE)
  }
  for(m in c(6,8,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }
    axis(1, at = seq(0,23,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE)
  }
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
}
dev.off()


pdf("timeCorrelationsIrriMajorRivers.pdf", width=8, height=7)
catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray", "Colorado")
catchTel = 0
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828, 7047)){
catchTel = catchTel + 1
sel = which(catchments == catch)[1]

for(v in 1){
  var = "evapDeficit"
  ref= refS
  for(m in c(1:3,6:8)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      time = 0
      if(m != 6 & m!= 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,"Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }
      if(m == 6){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,"_Coef_Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }      
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",ref,"_",var,"_lag_",lag,"_Coef_Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }      
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      nc_close(NC)
      for(i in 1:time){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
        lagDif = 0
        if(i/2 == ceiling(i/2)){
          lagDif = 0.5
        }
        leads[tel] = lag + lagDif
        sizes[tel] = ceiling(forecastMonth[i+lag*2]/2)
        if(ceiling(forecastMonth[i+lag*2]/2) == 1){
          tel = tel + 1
          vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
          leads[tel] = lag + lagDif
          sizes[tel] = 13
        }
        if(ceiling(forecastMonth[i+lag*2]/2) == 12){
          tel = tel + 1
          vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
          leads[tel] = lag + lagDif
          sizes[tel] = 0
        }
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }


  layout(makeMatrixDis(2,3))
  par(mar=c(0,0,0,0))
  for(v in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(v == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:3)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m], xlim=c(0.5,12.5), xaxs="i")
    for(l in 1:24){
      abline(a = -l+12-0.5, b = 1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1))
  }
  for(m in c(6,8,7)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m], xlim=c(0.5,12.5), xaxs="i")
    for(l in 1:24){
      abline(a = -l+12-0.1, b = 1, col="grey", lty=2)
    }
    axis(1, at = seq(1,12,1), labels=seq(1,12,1))
    contour(toPlot, add=TRUE, zlim=c(0,1))
  }
  NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/netcdf/",var,"_seasoAvg_annualCycle.nc", sep=""))
  data = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,24))
  nc_close(NC)
  if(var == "discharge"){
    selOutlet = which(catchments == catch)[which.max(data[,360:1,1][which(catchments == catch)])]
    multiply = 1
  }
  else{
    selOutlet = which(catchments == catch)
    multiply = area[,360:1][selOutlet]
    divide = sum(area[,360:1][selOutlet])
  }
  season = array()
  for(i in 1:24){
    season[i] = sum(data[,360:1,i][selOutlet]*multiply)/divide*1000
  }
  par(mar=c(4,4,2.0,0)) 
  plot(1,1, type="n", xlim=c(0,13),ylim=range(season), axes=FALSE, xaxs="i", xlab="Month", ylab=paste(var,"mm/d"), main=catchNames[catchTel])
  lines(seq(0,12,0.5), c(season, season[1]))
  axis(1, at=seq(0.5,11.5,1), labels=c(1:12))
  axis(2)
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}
}
dev.off()

pdf("timeCorrelationsScarcityMajorRivers.pdf", width=7, height=6)
catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray", "Colorado")
catchTel = 0
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828, 7047)){
catchTel = catchTel + 1
sel = which(catchments == catch)[1]

scenTel = 0
scen = array()
for(var in c("waterDemand")){
  ref= refS
  for(m in c(1:3,7)){
    scenTel = scenTel + 1
    scen[scenTel] = m
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      time = 0
      if(m != 6 & m!= 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,"Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }
      if(m == 6){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,"_Coef_Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }      
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",ref,"_",var,"_lag_",lag,"_Coef_Catchment.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))
      }      
      dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      nc_close(NC)
      for(i in 1:time){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][sel], na.rm=T)
        leads[tel] = lag
        sizes[tel] = forecastMonth[i+lag*2]
      }
    }
  interpolation[[scenTel]] <- interp(sizes, leads, vals)
  }
}

  layout(makeMatrixDis(2,2))
  par(mar=c(0,0,0,0))
  for(v in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    if(v == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:4)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=paste(modelNames[scen[m]], "Demand"))
    for(l in 1:24){
      abline(a = l, b = 1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1))
  }
#   for(m in c(5:8)){
#     par(mar=c(2,0.2,2.0,0)) 
#     toPlot = interpolation[[m]]
#     toPlot$z[toPlot$z < 0] = 0.0
#     image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=paste(modelNames[scen[m]], "Demand"))
#     for(l in 1:24){
#       abline(a = l, b = -1, col="grey", lty=2)
#     }
#     axis(1, at = seq(0,11,1), labels=seq(0,11,1))
#     contour(toPlot, add=TRUE, zlim=c(0,1))
#   }
  NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/netcdf/","evapDeficit","_seasoAvg_annualCycle.nc", sep=""))
  data = ncvar_get(NC, "evapDeficit", start=c(1,1,1), count=c(720,360,24))
  nc_close(NC)
  if(var == "discharge"){
    selOutlet = which(catchments == catch)[which.max(data[,360:1,1][which(catchments == catch)])]
    multiply = 1
  }
  else{
    selOutlet = which(catchments == catch)
    multiply = area[,360:1][selOutlet]
    divide = sum(area[,360:1][selOutlet])
  }
  season = array()
  for(i in 1:24){
    season[i] = sum(data[,360:1,i][selOutlet]*multiply)/divide*1000
  }
  par(mar=c(4,4,2.0,0)) 
  plot(1,1, type="n", xlim=c(0,13),ylim=range(season), axes=FALSE, xaxs="i", xlab="Month", ylab=paste(var,"mm/d"), main=catchNames[catchTel])
  lines(seq(0,12,0.5), c(season, season[1]))
  axis(1, at=seq(0.5,11.5,1), labels=c(1:12))
  axis(2)
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))
}
dev.off()



################# With Obs ##################

pdf("timeCorrelationsObsMajorRivers.pdf", width=10, height=12)

startMonths = c(12,rep(c(1:12), each=2))
startDays = c(16,rep(c(1,16), times=12))
endDays = c(32,rep(c(15,32), times=12))
lagSel = rep(c(1:24),times=2)

catchNames = c("Amazon", "Mississippi", "Bramaputhra", "Niger", "Darling", "Colorado", "Danube", "Rhine", "Orinoco")
catchTel = 0

for(catch in c(25, 36, 12, 6,53, 44,64,58,19)){
  Corrected = read.table("/tigress/nwanders/Scripts/hydroSeasonal/GRDC_DailyDischarge/Metadata_cor.csv", fill=T, header=T, sep=";")
  sel = which(Corrected[,1] == catch)
  file = Corrected[sel,2]
  if(length(sel) > 0){
    print(sel)
    tel = tel + 1
    Qdata = read.table(paste("/tigress/nwanders/Scripts/hydroSeasonal/GRDC_DailyDischarge",file,sep="/"), skip=40, header=T, sep=";")
    if(Corrected[sel,8] == "OK"  & Corrected[sel,9] == "OK"){
      Lat = Corrected[sel,4]
      Lon = Corrected[sel,3]
    }
    if(Corrected[sel,8] == "NOTOK"  & Corrected[sel,9] == "OK"){
      Lat = Corrected[sel,12]
      Lon = Corrected[sel,11]
    }
  }
  years = as.numeric(strftime(as.POSIXlt(as.Date(Qdata[,1], origin="1901-01-01")), format="%Y"))
  months = as.numeric(strftime(as.POSIXlt(as.Date(Qdata[,1], origin="1901-01-01")), format="%m"))
  days = as.numeric(strftime(as.POSIXlt(as.Date(Qdata[,1], origin="1901-01-01")), format="%d"))
  Qobs = Qdata[,4]
  Qobs[Qobs < 0] = NA
  Qcycle = array()
  vals = array()
  leads = array()
  sizes = array()
  tel = 0
  for(i in 1:24){
    sel = which(months == startMonths[i] & days == startDays[i])
    for(lag in 0:11){
      tel = tel + 1
      j = lagSel[i+lag*2+1]
      selLag = which(months == startMonths[j] & days == startDays[j])
      if(sel[1] > selLag[1]){
        selLag = selLag[2:length(sel)]
      }
      else{
        selLag = selLag[1:(length(sel)-1)]
      }
      meanQobs = array()
      meanQlagObs = array()
      for(t in 1:length(selLag)){
        meanQobs[t] = mean(Qobs[sel[t]:(sel[t]+14)], na.rm=T)
        meanQlagObs[t] = mean(Qobs[selLag[t]:(selLag[t]+14)], na.rm=T)
      }
      vals[tel] = cor(meanQobs, meanQlagObs, use="complete.obs")
      leads[tel] = lag+((i+1)/2)-floor(i/2)
      sizes[tel] = ceiling(i/2)
      if(lag == 0){
        Qcycle[i] = mean(meanQlagObs, na.rm=T)
      }
    }
  }
  obsACF <- interp(sizes, leads, vals)

  catchTel = catchTel + 1

  var = varnameS[1]
  ref= refS
  for(m in c(1:4,6:8,9)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      if(m < 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      }
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/PCRGLOBWB_",ref,"_",var,"_lag_",lag,"_Coef.nc4", sep=""))
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,24))
      }
      if(m == 9){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/PCRGLOBWB/netcdf/1/output_1981_2012_ACF_lag_",lag,".nc", sep=""))
        dataCOR = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,24))
      }
      nc_close(NC)
      for(i in 1:24){
        tel = tel + 1
        vals[tel] = dataCOR[,360:1,i][(Lon*2+360.5), (Lat*2+180.5)]
        leads[tel] = lag+((i+1)/2)-floor(i/2)
        sizes[tel] = ceiling(i/2)
      }
    }
  interpolation[[m]] <- interp(sizes, leads, vals)
  }
  
  for(m in c(1:4,6,8,7,9)){
    leads = array()
    sizes = array()
    vals = array()
    latVals = array()
    lats = array()
    latLeads = array()
    latTel = 0
    tel = 0
    for(lag in 0:11){
      print(lag)
      model = modelS[m]
      if(m < 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_",model,"_",ref,"_",var,"_lag_",lag,".nc4", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      }
      if(m == 8){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/resultsNetCDF/VIC_",ref,"_",var,"_lag_",lag,"_Coef.nc4", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, "correlation", start=c(1,1,1), count=c(720,360,time))
      }
      if(m == 9){
        NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/PCRGLOBWB/netcdf/1/output_1981_2012_ACF_lag_",lag,".nc", sep=""))
        time= length(ncvar_get(NC, "time"))        
        dataCOR = ncvar_get(NC, var, start=c(1,1,1), count=c(720,360,time))
        dataCOR[dataCOR < -1] = NA
        dataCOR[dataCOR > 1] = NA        
      }
      nc_close(NC)
      for(i in 1:time){
        tel = tel + 1
        vals[tel] = median(dataCOR[,360:1,i][(Lon*2+360.5), (Lat*2+180.5)], na.rm=T)
        leads[tel] = lag
        sizes[tel] = i/2
      }
    } 
  interpolation[[m+9]] <- interp(sizes, leads, vals)
  }


  A = makeMatrixDisLarge(4,4)
  layout(A)
  par(mar=c(0,0,0,0))
  for(i in c(1:4)){
    par(mar=c(2.0,0.2,2.0,0))
    if(i == 2){
      par(mar=c(2,0.2,2.0,0)) 
    }
    plot(1,1,type="n", xlim=c(0,1), ylim=c(-0.15,12.2), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="")
    text(0.75, seq(0.1,11.75, length=12), seq(0,11,1))
    text(0.22, 6.0, "Lead time", srt=90)
  }
  for(m in c(1:3,4)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  for(m in c(6,8,7,9)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  for(m in c(1:3,4)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = interpolation[[m+9]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }    
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  for(m in c(6,8,7,9)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = interpolation[[m+9]]
    toPlot$z[toPlot$z < 0] = 0.0
    image(toPlot, axes=FALSE, zlim=c(0,1), col=tim.colors(200), main=modelNames[m])
    for(l in 1:24){
      abline(a = l, b = -1, col="grey", lty=2)
    }
    axis(1, at = seq(0,11,1), labels=seq(0,11,1))
    contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  }
  NC = nc_open(paste("/tigress/nwanders/Scripts/hydroSeasonal/PGF/VIC/","discharge","_seasoAvg_annualCycle.nc", sep=""))
  dataVIC = ncvar_get(NC, var, start=c((Lon*2+360.5), 361-(Lat*2+180.5),1), count=c(1,1,24))
  nc_close(NC)
  NC = nc_open("/tigress/nwanders/Scripts/hydroSeasonal/PGF/netcdf/discharge_seasoAvg_annualCycle.nc")
  dataPCR = ncvar_get(NC, var, start=c((Lon*2+360.5), 361-(Lat*2+180.5),1), count=c(1,1,24))
  nc_close(NC)
  par(mar=c(4,4,2.0,0)) 
  plot(1,1, type="n", xlim=c(0,13),ylim=range(c(Qcycle, dataVIC, dataPCR), na.rm=T), axes=FALSE, xaxs="i", xlab="Month", ylab=paste(var,"m3/s"), main=catchNames[catchTel])
  lines(seq(0,12,0.5), c(Qcycle, Qcycle[1]), col=1, lwd=2)
  lines(seq(0,12,0.5), c(dataVIC, dataVIC[1]), col=2)
  lines(seq(0,12,0.5), c(dataPCR, dataPCR[1]), col=3)
  legend("bottomleft", c("Obs", "VIC", "PCR"), col=c(1,2,3), lwd=c(2,1,1))
  axis(1, at=seq(0.5,11.5,1), labels=c(1:12))
  axis(2)

  par(mar=c(2,4,2,2))
  toPlot = obsACF
  toPlot$z[toPlot$z <0] = 0
  image(toPlot, zlim=c(0,1), col=tim.colors(200), main=catchNames[catchTel], ylab="Lead Time", xlab="", axes=FALSE)
  for(l in 1:24){
    abline(a = l, b = -1, col="grey", lty=2)
  }
  axis(1, at = seq(0,11,1), labels=seq(0,11,1))
  contour(toPlot, add=TRUE, zlim=c(0,1), nlevels=10)
  
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(0,1,0.1), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))

}

dev.off()





