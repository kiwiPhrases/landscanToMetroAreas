require(raster)
require(igraph)

setwd('C:/Users/SpiffyApple/Documents/USC/Courses/Fall2018/JorgeWinterWork/LandScan Global 2015/LandScan Global 2015/lspop2015')

## load the raster data for Peru
layer = raster("Peru landscan 2.tif")

## set key parameters
coreSize = 50e3
settlementSize = 5e3
fringeMin = 500
fringeMax = 1500
coreMin = 1500

## -------------------- Smoothing ------------------------- ##
## maka a weight matrix to use in focal()
W = matrix(1,7,7)

smthLyr = stack(focal(layer,w=W, fun='mean')) ## apply smoothing function to raster

#plot the result
par(mfrow=c(1,2))
spplot(smthLyr)
spplot(layer)

## --------------------- find cores ----------------------- ##
coreClumps = clump(smthLyr>=coreMin, directions=4) # find clumps of core-density cells

coreClumpSums = zonal(smthLyr, coreClumps, fun='sum', na.rm=TRUE) # find total pop in each clump
colnames(coreClumpSums) = c("id", "v")
coreClumpsVals = subs(coreClumps,data.frame(coreClumpSums)) # sub in pop sums into zone nums
coreAreas = coreClumpsVals>=coreSize # select min cores
plot(coreAreas)## plot results
paste("Number of core area clumps:",nrow(coreClumpSums[coreClumpSums[,2]>50e3,]))

# dump the core data to a data frame and csv file
coreCentroids = t(sapply(unique(coreClumps[coreClumpsVals>coreSize]), function(x){colMeans(xyFromCell(coreClumps, Which(coreClumps==x,cell=TRUE)))}))
coreDensity = sapply(unique(coreClumps[coreClumpsVals>coreSize]), function(x){length(Which(coreClumps==x,cell=TRUE))})
coreAreaDF = cbind(coreCentroids, coreDensity, coreClumpSums[,2][coreClumpSums[,2]>50e3])
colnames(coreAreaDF)[3:4] = c("sqrkm", "pop")
coreAreaDF = data.frame(coreAreaDF)
coreAreaDF[,'density'] = coreAreaDF$pop/coreAreaDF$sqrkm
write.csv(coreAreaDF,'coreArea.csv')

## ------------------- find fringe areas ------------------ ##
fringePotentials = (smthLyr >= fringeMin) & (smthLyr < fringeMax) # set fringe cell densities


# first, find cells wit fringe density near core
# then, find cells with fringe density adjacent to fringe cells
# iterate to expand the fringe cells
fringeAreaCells = adjacent(fringePotentials,Which(coreAreas == 1, cells = TRUE),
                           target = Which(fringePotentials == 1, cells = TRUE), directions=4, include=FALSE)
listLength = 20
fringeAreaList = list(length=listLength)
fringeAreaList[[1]] = fringeAreaCells[,2]

for(i in 2:listLength){
  print(i)
  fringeAreaList[[i]] = adjacent(fringePotentials,fringeAreaList[[i-1]],pair=FALSE,
                             target = Which(fringePotentials == 1, cells = TRUE), directions=4, include=FALSE)}
paste("Number of fringe cells:", length(unique(unlist(fringeAreaList))))

fringeAreas =  smthLyr ## initiate an empty layer, populate with 0s
fringeAreas[] = 0
fringeAreas[unique(unlist(fringeAreaList))] = 1 ## set 1 at fringe cells
plot(fringeAreas)

# associate fringe areas to Cores by:
# find fringe clumps
# compute dist btw core clump centroids and fringe clump centroids
# for each fringe clump, find closest clump
# associate closest core clump with each fringe clump

fringeClumps = clump(fringeAreas[[1]], directions=4)
## find center of each fringe:
fringeCentroids = t(sapply(c(1:maxValue(fringeClumps)),
                           function(x){colMeans(xyFromCell(fringeClumps, Which(fringeClumps==x,cell=TRUE)))}))

# plot this stuff -> can't believe i have to do this
plot(fringeCentroids,col='blue',pch=17)
points(t(coreCentroids),col='red',pch=16)

# wow R, I manually compute distance matrix. thankfully, this is rather small
Xs = sapply(fringeCentroids[,'x'],function(z){(z-coreCentroids[,'x'])^2})
Ys = sapply(fringeCentroids[,'y'],function(z){(z-coreCentroids[,'y'])^2})
distMat = sqrt(Xs+Ys)

fringeClumpCoreNums = apply(distMat,2, which.min)
fringeClumpCores = subs(fringeClumps,data.frame(id=c(1:maxValue(fringeClumps)), value=fringeClumpCoreNums))


# get pop and size of each fringe (as summed by corresponding core)
fringeSize = zonal(fringeAreas,fringeClumpCores, fun='sum', na.rm=TRUE) #i think?
#fringeSize = sapply(c(1:maxValue(fringeClumpCores)),function(x)length(Which(fringeClumpCores==x,cells=TRUE)))
fringePop = zonal(smthLyr,fringeClumpCores, fun='sum', na.rm=TRUE)
fringeDF = data.frame(size= fringeSize, pop = fringePop[,'value'], coreCentroids)
colnames(fringeDF)[3:4] = c('coreX', 'coreY')
write.csv(fringeDF,"fringeAreas.csv")

##-------------------- find settlements ------------------- ##
# find clusters of settlement-density cells
settlementClumps = clump((smthLyr >= fringeMin) & (smthLyr < fringeMax), directions=4)
paste("Number of potentials settlements:",maxValue(settlementClumps))
settlementClumpSums = zonal(smthLyr, settlementClumps, fun='sum', na.rm=TRUE) #compute sums of clusters
colnames(settlementClumpSums) = c("id", "v")
settlementClumpsVals = subs(settlementClumps,data.frame(settlementClumpSums)) #replace zone vals with pop sums
settlementPotentials = settlementClumpsVals >= settlementSize
settlementCells = adjacent(settlementPotentials,Which((coreAreas ==1) | (fringeAreas==1), cells = TRUE),
                           target = Which(settlementPotentials==1, cells = TRUE), directions=4, include=FALSE)

settlementAreas = settlementPotentials
settlementAreas[settlementCells[,2]] = 0 ## annull the settlement potentials cells adjacent to fringes and cores

newClumpSet = clump(settlementAreas,directions=4)
paste("Final settlement count:",maxValue(newClumpSet))

settlementCentroids = sapply(c(1:maxValue(newClumpSet)), function(x){colMeans(xyFromCell(newClumpSet, Which(newClumpSet==x,cell=TRUE)))})
settlementSize = sapply(c(1:maxValue(newClumpSet)), function(x){length(Which(newClumpSet==x,cell=TRUE))})
settlementPop = zonal(smthLyr, newClumpSet, fun='sum', na.rm=TRUE)

settlementMat =cbind(t(settlementCentroids), settlementSize, settlementPop)
settlementDF = data.frame(settlementMat)
colnames(settlementDF)[5]='pop'
write.csv(settlementDF, 'settlementAreas.csv')

## -------------------- plot the results ----------------- ##
coreAreas[coreAreas == 0] = NA
fringeAreas[fringeAreas == 0] = NA
settlementAreas[settlementAreas == 0] = NA

fringeAreaplt = fringeAreas*2
settlementAreaplt = settlementAreas*3

areas =cover(cover(coreAreas, fringeAreaplt[[1]]),settlementAreaplt[[1]])

par(mfrow=c(1,1))
plot(smthLyr)
plot(areas, col = colorRampPalette(c('red','purple','blue'))(10))
writeRaster(areas, 'peruAreas.tif')

# write results to disc:

writeRaster(coreAreas, 'peruCores.tif')
writeRaster(fringeAreas, 'peruFringe.tif')
writeRaster(settlementAreas, 'peruSettlements.tif')

#####---------------------- test area ----------------------############
# make data
set.seed(2)
testVals = matrix(runif(100,1,100),nrow=10,ncol=10)
testVals[testVals<25] = 0
testVals[(testVals>=25) & (testVals<50)] = 25
testVals[(testVals>=50) & (testVals<=75)] = 50
testVals[(testVals>75)] = 100

test = raster(nrow=10,ncol=10)
values(test) = testVals

# find cores:
coreClumps = clump(test>=90, directions=4)

coreClumpSums = zonal(test, coreClumps, fun='sum', na.rm=TRUE)
colnames(coreClumpSums) = c("id", "v")
coreClumpsVals = subs(coreClumps,data.frame(coreClumpSums))
coreAreas = coreClumpsVals>200

# find fringe areas:
plot(fringePotentials)
plot(coreAreas)

fringePotentials = (test == 50)
fringeAreaCells = adjacent(fringePotentials,Which(coreAreas == 1, cells = TRUE),
                           target = Which(fringePotentials == 1, cells = TRUE), directions=4, include=FALSE)
length(fringeAreaCells[,1])

fringeAreas =  test
fringeAreas[] = 0
fringeAreas[fringeAreaCells[,2]] = 1


#find settlements:
settlementClumps = clump(test==25, directions=4)
settlementClumpSums = zonal(test, settlementClumps, fun='sum', na.rm=TRUE)
colnames(settlementClumpSums) = c("id", "v")
settlementClumpsVals = subs(settlementClumps,data.frame(settlementClumpSums))
settlementPotentials = settlementClumpsVals >= 75
settlementCells = adjacent(settlementPotentials,Which((coreAreas ==1) | (fringeAreas==1), cells = TRUE),
                           target = Which(settlementPotentials==1, cells = TRUE), directions=4, include=FALSE)

settlementAreas = settlementPotentials
settlementAreas[settlementCells[,2]] = 0 ## annull the settlement potentials cells adjacent to fringes and cores
plot(settlementAreas)

# plot the results:
coreAreas[coreAreas == 0] = NA
fringeAreas[fringeAreas == 0] = NA
settlementAreas[settlementAreas == 0] = NA

fringeAreaplt = fringeAreas*2
settlementAreaplt = settlementAreas*3

areas =cover(cover(coreAreas, fringeAreaplt),settlementAreaplt)

par(mfrow=c(1,2))
plot(test,col = colorRampPalette(c('white','blue','purple','red'))(10))
plot(areas, col = colorRampPalette(c('red','purple','blue'))(10))


