require(rgdal)
require(ncdf4)
require(fields)
require(akima)

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

area = matrix(readGDAL("/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/routing/catchmentarea_30min.map")@data$band1, 720, 360)
NC = nc_open("/tigress/nwanders/Scripts/PCR-GLOBWB/input30min/global/majorCatchments.nc")
catchments = ncvar_get(NC, "Band1")
nc_close(NC)

catchNames = c("Amazon", "Mississippi", "Ganges", "Indus", "Nile", "Niger", "Murray", "Colorado")
catchTel = 0
selCoords = matrix(NA, 8,2)
for(catch in c(9131, 7225, 7663, 7593, 7064, 8942,10828, 7047)){
  catchTel = catchTel + 1
  selCoords[catchTel,] = which.max(catchments == catch, arr.ind=T)[1,]
}

lonCell = selCoords[[8,1]
latCell = selCoords[[8,1]
lonCell = 249
latCell = 184
modelS = c("CanCM3", "CanCM4", "FLOR","ESP")
varnameS = c("discharge") #, "storGroundwater", "satDegUpp","satDegLow")
varncS = c()c("discharge") #, "groundwater_storage", "upper_soil_saturation_degree","lower_soil_saturation_degree")
refS = c("PGF")
totEns = c(10,10,12,9,9,10,9,9)
meanS = array(NA, c(4,372,24))
cycleMat = matrix(NA, 372, 24)
  var = "discharge"
  ref = refS
  NC = nc_open(paste("../",ref,"/netcdf/",var,"_seasoAvg_annualCycle.nc", sep=""))
  dataCycle = ncvar_get(NC, varncS[v], start=c(lonCell,latCell,1), count=c(1,1,24))
  nc_close(NC)
  startDay = (as.numeric(1981)-1979)*24+(as.numeric(1)-1)*2 + 1
  fileName1 = paste("discharge","_seasoAvg_output.nc",sep="")
  NC = nc_open(paste("../",ref,"/PCRGLOBWB/netcdf/1/",fileName1, sep=""))
  dataRef1 = ncvar_get(NC, "discharge", start=c(lonCell,latCell,startDay), count=c(1,1,744))
  nc_close(NC)
  #fileName2 = paste("totalEvaporation","_seasoAvg_output.nc",sep="")
  #NC = nc_open(paste("../",ref,"/PCRGLOBWB/netcdf/1/",fileName2, sep=""))
  #dataRef2 = ncvar_get(NC, "total_evaporation", start=c(lonCell,latCell,startDay), count=c(1,1,744))
  #nc_close(NC)
  dataRef = dataRef1 #- dataRef2
  forecastTel = 0
  for(year in 1981:2012){
    for(month in c("01","02","03","04","05","06","07","08","09","10","11","12")){
      forecastTel = forecastTel + 1
      cycleMat[forecastTel,] = dataCycle[c(((as.numeric(month)-1)*2 + 1):24,1:((as.numeric(month)-1)*2))][1:24]
      startDay = (as.numeric(year)-1980)*24+(as.numeric(month)-1)*2 + 1
      for(m in 1:length(modelS)){
        print(paste(m, forecastTel))
        nEns = totEns[m]
        model = modelS[m]
        ensMean1 = array(0, 24)
        for(n in c(1:nEns)){
          NC = nc_open(paste("../",model,"/PCRGLOBWB/",ref,"/",year,"-",month,"-01/netcdf/",n,"/",fileName1, sep=""))
          data = ncvar_get(NC, "discharge", start=c(lonCell,latCell,1), count=c(1,1,24))
          ensMean1 = ensMean1 + data/nEns
          nc_close(NC)
        }
        #ensMean2 = array(0, 24)
        #for(n in c(1:nEns)){
        #  NC = nc_open(paste("../",model,"/PCRGLOBWB/",ref,"/",year,"-",month,"-01/netcdf/",n,"/",fileName2, sep=""))
        #  data = ncvar_get(NC, "total_evaporation", start=c(lonCell,latCell,1), count=c(1,1,24))
        #  ensMean2 = ensMean2 + data/nEns
        #  nc_close(NC)
        #}
        meanS[m,forecastTel,] = ensMean1 #-ensMean2
      }
    }
  }
}
pdf("AmazonDischargeForecast.pdf", width=7, height=8)
  toPlot = dataRef - rep(dataCycle,times=31)
  zlims = c(-max(abs(toPlot), na.rm=T), max(abs(toPlot), na.rm=T))
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
  for(m in c(1:2)){
    par(mar=c(2.0,0.2,2.0,0))
    toPlot = meanS[m,,] - cycleMat
    image(toPlot, axes=FALSE, zlim=zlims, col=tim.colors(200), main=modelS[m])
  }
  for(m in c(3:4)){
    par(mar=c(2,0.2,2.0,0)) 
    toPlot = meanS[m,,] - cycleMat
    image(toPlot, axes=FALSE, zlim=zlims, col=tim.colors(200), main=modelS[m])
    axis(1, at = seq(0,1,length=7), labels=seq(1980,2010,5))
  }
  par(mar=c(4,4,2.0,0)) 
  toPlot = dataRef - rep(dataCycle,times=31)
  image(as.matrix(toPlot,1, 744), zlim=zlims, col=tim.colors(200), axes=FALSE)
  axis(1, at = seq(0,1,length=7), labels=seq(1980,2010,5))
  par(mar=c(1,2.7,1,3.5))
  colLen = 100
  cols = tim.colors
  plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="R")
  axis(4,seq(0,1.0,0.1), seq(zlims[1],zlims[2],length=11), las=1)
  symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols(colLen), bg=cols(colLen))
dev.off()
