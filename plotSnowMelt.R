#####R##########
require(ncdf)
require(fields)

mapFlip <- function(map){
  temp = matrix(NA, dim(map)[1], dim(map)[2])
  half = dim(map)[1]/2
  end = dim(map)[1]
  temp[1:half,]=map[(half+1):end,]
  temp[(half+1):end,]=map[1:half,]
  return(temp)
}

continentNC = open.ncdf("continents.nc")
continent = get.var.ncdf(continentNC, "con")
close.ncdf(continentNC)
continent[continent < 1] = NA
continent = mapFlip(continent)

makeMatrix <- function(rows, cols, size=5, legendSize = 2, labelSize = 1){
  A = matrix(rows*cols+3, rows*size, cols*size+legendSize+labelSize)
  A[1:size,1:labelSize] = 1
  A[(size+1):(2*size),1:labelSize] = 2
  for(i in c(1:(rows*cols))){
    row = floor((i-1)/cols)
    col = (i-1) - row * cols
    A[(row*size+1):(row*size+size),(col*size+1):(col*size+size)+labelSize] = i + 2
  }
  return(A)
}

NC = open.ncdf("OffsetSnowSeason.nc")
onsetPGF = get.var.ncdf(NC, "startPGF")
offsetCanCM3 = get.var.ncdf(NC, "offsetCanCM3")
offsetCanCM4 = get.var.ncdf(NC, "offsetCanCM4")
offsetFLOR = get.var.ncdf(NC, "offsetFLOR")
offsetCCSM = get.var.ncdf(NC, "offsetCCSM")
close.ncdf(NC)

minPGF = matrix(999, 360,180)
maxPGF = matrix(0, 360,180)
snowYears = matrix(0, 360,180)
for(i in 1:32){
  temp = onsetPGF[,,i]
  temp[temp == 0] = NA
  sel = which(temp < minPGF & is.na(temp) == FALSE)
  minPGF[sel] = onsetPGF[,,i][sel]
  sel = which(temp > maxPGF & is.na(temp) == FALSE)
  maxPGF[sel] = onsetPGF[,,i][sel]
  sel = which(is.na(temp) == FALSE)
  snowYears[sel] = snowYears[sel] + 1
}

rangePGF = maxPGF - minPGF
rangePGF[rangePGF == 0.0] = NA

Brier = list()
Brier[[1]] = pmax(rowMeans(offsetCanCM3, dims=2),0.0)/rangePGF
Brier[[2]] = pmax(rowMeans(offsetCanCM4, dims=2),0.0)/rangePGF
Brier[[3]] = pmax(rowMeans(offsetFLOR, dims=2),0.0)/rangePGF
Brier[[4]] = pmax(rowMeans(offsetCCSM, dims=2),0.0)/rangePGF

pdf("MeltOnsetPGF.pdf", width=6.5, height=5)
A = matrix(1, 2, 9)
A[1,9] = 2
A[2,] = 3
A[2,9] = 4
layout(A)
colLen = 100
cols=colorRampPalette(c("yellow","Green","Blue", "Purple", "Black"))(100)
par(mar=c(0,0,2,0))
toPlot = rowMeans(onsetPGF, dims=2)
toPlot[toPlot > 200] = 200
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main="Average onset of snow melt")
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
world(add=T)
par(mar=c(1,1.8,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Julian day")
axis(4, seq(0,1.0,0.1),seq(0,200,20), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)

par(mar=c(0,0,2,0))
toPlot = rangePGF
toPlot[toPlot > 200] = 200
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main="Range in snow melt onset")
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
world(add=T)
par(mar=c(1,1.8,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Days")
axis(4, seq(0,1.0,0.1),seq(0,200,20), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)

dev.off()

pdf("BrierSnowMelt.pdf", width=6.5, height=6)
A = makeMatrix(3,1,labelSize=0, size=15)
A = A - 2
layout(A)
colLen = 100
cols=colorRampPalette(c("Blue","Yellow", "red"))(100)
par(mar=c(0,0,2,0))
models = c("CanCM3", "CanCM4", "FLOR", "CCSM")
for(model in 1:3){
  toPlot = Brier[[model]]
  toPlot[toPlot > 1] = 1
  image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main=models[model])
  image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(0,1), ylim=c(25,87), xlim=c(-180,180), add=T)
  image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(0,1), ylim=c(25,87), xlim=c(-180,180), add=T)
  world(add=T)
}
par(mar=c(1,1.6,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Brier")
axis(4, seq(0,1.0,0.1),seq(0,1,0.1), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)
dev.off()

#########Season Start################

NC = open.ncdf("startSnowSeason.nc")
onsetPGF = get.var.ncdf(NC, "startPGF")
offsetCanCM3 = get.var.ncdf(NC, "offsetCanCM3")
offsetCanCM4 = get.var.ncdf(NC, "offsetCanCM4")
#offsetFLOR = get.var.ncdf(NC, "offsetFLOR")
#offsetCCSM = get.var.ncdf(NC, "offsetCCSM")
close.ncdf(NC)

minPGF = matrix(999, 360,180)
maxPGF = matrix(0, 360,180)
snowYears = matrix(0, 360,180)
for(i in 1:32){
  temp = onsetPGF[,,i]
  temp[temp == 0] = NA
  sel = which(temp < minPGF & is.na(temp) == FALSE)
  minPGF[sel] = onsetPGF[,,i][sel]
  sel = which(temp > maxPGF & is.na(temp) == FALSE)
  maxPGF[sel] = onsetPGF[,,i][sel]
  sel = which(is.na(temp) == FALSE)
  snowYears[sel] = snowYears[sel] + 1
}

rangePGF = maxPGF - minPGF
rangePGF[rangePGF == 0.0] = NA

Brier = list()
Brier[[1]] = pmax(rowMeans(offsetCanCM3, dims=2),0.0)/rangePGF
Brier[[2]] = pmax(rowMeans(offsetCanCM4, dims=2),0.0)/rangePGF
Brier[[3]] = pmax(rowMeans(offsetFLOR, dims=2),0.0)/rangePGF
Brier[[4]] = pmax(rowMeans(offsetCCSM, dims=2),0.0)/rangePGF

pdf("AccuOnsetPGF.pdf", width=6.5, height=5)
A = matrix(1, 2, 9)
A[1,9] = 2
A[2,] = 3
A[2,9] = 4
layout(A)
colLen = 100
cols=colorRampPalette(c("yellow","Green","Blue", "Purple", "Black"))(100)
par(mar=c(0,0,2,0))
toPlot = onsetPGF
toPlot[toPlot == 0] = NA
toPlot[snowYears < 15] = NA
toPlot = rowMeans(toPlot, dims=2, na.rm=T)
toPlot[toPlot > 200] = 200
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main="Average onset of snow accumulation")
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
world(add=T)
par(mar=c(1,1.8,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Julian day")
axis(4, seq(0,1.0,0.1),seq(0,200,20), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)

par(mar=c(0,0,2,0))
toPlot = rangePGF
toPlot[toPlot > 200] = 200
toPlot[snowYears < 15] = NA
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main="Range in snow accumulation onset")
image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(1,200), ylim=c(20,87), xlim=c(-180,180), add=T)
world(add=T)
par(mar=c(1,1.8,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Days")
axis(4, seq(0,1.0,0.1),seq(0,200,20), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)

dev.off()

pdf("BrierSnowAccu.pdf", width=6.5, height=6)
A = makeMatrix(3,1,labelSize=0, size=15)
A = A - 2
layout(A)
colLen = 100
cols=colorRampPalette(c("Blue","Yellow", "red"))(100)
par(mar=c(0,0,2,0))
models = c("CanCM3", "CanCM4", "FLOR", "CCSM")
for(model in 1:3){
  toPlot = Brier[[model]]
  toPlot[toPlot > 1] = 1
  toPlot[snowYears < 25] = NA
  image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), continent, ylim=c(25,87), xlim=c(-170,190), zlim=c(0,10), col="grey", axes=FALSE, xlab="", ylab="", main=models[model])
  image(seq(-179.5,179.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(0,1), ylim=c(25,87), xlim=c(-180,180), add=T)
  image(seq(180.5,539.5,1), seq(-89.5, 89.5,1), toPlot, col=cols, zlim=c(0,1), ylim=c(25,87), xlim=c(-180,180), add=T)
  world(add=T)
}
par(mar=c(1,1.6,2,2.3))
plot(1,1,type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, ylab="", xlab="", main="Brier")
axis(4, seq(0,1.0,0.1),seq(0,1,0.1), las=1, tick=FALSE, line=-0.7)
symbols(rep(0.5,colLen), seq(0,1,length=colLen), rectangles=matrix(rep(c(1,1/colLen),2*colLen),colLen,2, byrow=T), inches=F, add=T, fg=cols, bg=cols)
dev.off()
