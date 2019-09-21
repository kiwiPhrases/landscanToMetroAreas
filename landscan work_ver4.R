require(raster)
require(igraph)
require("rgdal")
require(rgeos)
home_path = 'C:/Users/SpiffyApple/Documents/USC/Courses/Fall2018/JorgeWinterWork/LandScan Global 2015/LandScan Global 2015/lspop2015'
setwd(home_path)

## convenience function
reName = function(fileName,strInsert){
  old = strsplit(fileName,split = '\\.')
  return(paste(old[[1]][1],"_",strInsert,".",old[[1]][2], sep=""))
}

## load the raster data for Peru
layer = raster("Peru landscan 2.tif")

## load Peru border
peru = readOGR(dsn="./peruBorder")

## iterate over key parameters
#matDists = c(5,7,9)
matDists = c(7)
#fringeMin = c(250,500)
fringeMin = c(250)
for(ij in matDists){
  for(fringeMin in fringeMin){
    print(paste('smoothing distance:',ij,'fringe thresh:' ,fringeMin))

    # reset home path
    home_path = 'C:/Users/SpiffyApple/Documents/USC/Courses/Fall2018/JorgeWinterWork/LandScan Global 2015/LandScan Global 2015/lspop2015'
    setwd(home_path)
## -------------------- Smoothing ------------------------- ##
# set smoothing distance
#ij = 5

# maka a weight matrix to use in focal()
W = matrix(1,ij,ij)

smthLyr = stack(focal(layer,w=W, fun='mean')) ## apply smoothing function to raster
print("smoothing layer")
smth = smthLyr[[1]]
smth[smth == 0] = NA
smth[smth < 11] = NA ## 3rd quartile

## mask(cut out) ocean parts
smth = mask(smth, peru)

## save it as a polygon because arcGIS
#tmpVal = rasterToPolygons(smth, dissolve=FALSE)
#writeOGR(obj = tmpVal, dsn = "PeruLandScan",layer="peru",driver="ESRI Shapefile")

#plot the result
#par(mfrow=c(1,2))
#spplot(smthLyr)
#spplot(layer)

## set key parameters
coreSize = 50e3
settlementSize = 10e3
#fringeMin = 500
fringeMax = 1500
coreMin = 1500

##--------------- Shift File Destinations ---------------- ##
home_path = 'C:/Users/SpiffyApple/Documents/USC/Courses/Fall2018/JorgeWinterWork/LandScan Global 2015/LandScan Global 2015/lspop2015'
newFolder = paste("fringe",fringeMin, "smth",ij, sep='')
dir.create(newFolder)
new_path = paste(home_path, newFolder)
print(paste("Setting new path to", new_path))
setwd(paste(home_path, newFolder, sep='/'))

#for(fringeMin in seq(250,750,by=50)){
print("Finding core,fringe, and rural areas")
## --------------------- find cores ----------------------- ##
coreClumps = clump(smthLyr>=coreMin, directions=4) # find clumps of core-density cells

coreClumpSums = zonal(smthLyr, coreClumps, fun='sum', na.rm=TRUE) # find total pop in each clump
colnames(coreClumpSums) = c("id", "v")
coreClumpsVals = subs(coreClumps,data.frame(coreClumpSums)) # sub in pop sums into zone nums
coreAreas = coreClumpsVals>=coreSize # select min cores
#plot(coreAreas)## plot results
paste("Number of core area clumps:",nrow(coreClumpSums[coreClumpSums[,2]>coreSize,]))

# dump the core data to a data frame and csv file
coreCentroids = t(sapply(unique(coreClumps[coreClumpsVals>coreSize]), function(x){colMeans(xyFromCell(coreClumps, Which(coreClumps==x,cell=TRUE)))}))
coreDensity = sapply(unique(coreClumps[coreClumpsVals>coreSize]), function(x){length(Which(coreClumps==x,cell=TRUE))})
coreAreaDF = cbind(coreCentroids, coreDensity, coreClumpSums[coreClumpSums[,2]>coreSize,])
colnames(coreAreaDF)[3:5] = c("sqrkm", "zone","pop")
coreAreaDF = data.frame(coreAreaDF)
#coreAreaDF[,'density'] = coreAreaDF$pop/coreAreaDF$sqrkm
coreAreaDF[,'areaType'] = 'core'

#coreAreaDFsp = SpatialPointsDataFrame(coreAreaDF[,c('x', 'y')], coreAreaDF[,3:4], proj4string = crs(coreAreas))
#coreAreaDF[,'zone'] = NA ## add this to make colnames consistent across data framers
#write.csv(coreAreaDF,'coreArea.csv')

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
  #print(i)
  fringeAreaList[[i]] = adjacent(fringePotentials,fringeAreaList[[i-1]],pair=FALSE,
                             target = Which(fringePotentials == 1, cells = TRUE), directions=4, include=FALSE)}
paste("Number of fringe cells:", length(unique(unlist(fringeAreaList))))

fringeAreas =  smthLyr ## initiate an empty layer, populate with 0s
fringeAreas[] = NA
fringeAreas[unique(unlist(fringeAreaList))] = 1


## associate fringes to cores
coreFringes = merge(fringeAreas,coreAreas) ## merge fringes to cores
coreFringeClumps = clump(coreFringes==1,directions=8)

coreFringeClumpPop = zonal(smthLyr,coreFringeClumps,fun='sum', na.rm=TRUE) ## compute pop
coreFringeClumpSize = zonal(coreFringes==1,coreFringeClumps,fun='sum',na.rm=TRUE)
all(coreClumpSums[coreClumpSums[,2]>=coreSize,2] <= coreFringeClumpPop[,2]) #this should be true

coreFringeDF = cbind(coreClumpSums[coreClumpSums[,2]>=coreSize,], coreFringeClumpPop[,2],coreFringeClumpSize[,2])
colnames(coreFringeDF) = c('zone','corePop','coreFringePop','coreFringeSize')

## make the fringe data frame
fringeDF = data.frame(zone= coreFringeDF[,'zone'],
                 pop = coreFringeDF[,'coreFringePop']-coreFringeDF[,'corePop'],
                 sqrkm = coreFringeDF[,'coreFringeSize'] - coreAreaDF[,'sqrkm'],
                 coreAreaDF[,c('x','y')])
fringeDF[,'areaType'] = 'fringe'

## make the coreFringe data frame
coreFringeDF = data.frame(coreFringeDF[,c('zone','coreFringePop','coreFringeSize')])
coreFringeDF['areaType'] = 'core&fringe'
coreFringeDF = cbind(coreFringeDF,coreAreaDF[,c('x','y')])
colnames(coreFringeDF)[2:3] = c('pop','sqrkm')
#write.csv(fringeDF,"fringeAreas.csv")

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
settlementSizes = sapply(c(1:maxValue(newClumpSet)), function(x){length(Which(newClumpSet==x,cell=TRUE))})
settlementPop = zonal(smthLyr, newClumpSet, fun='sum', na.rm=TRUE)

settlementMat =cbind(t(settlementCentroids), settlementSizes, settlementPop)
settlementDF = data.frame(settlementMat)
colnames(settlementDF)[c(3,5)]=c('sqrkm','pop')
settlementDF[,'areaType'] = 'settlement'
settlementDF[,'zone'] = settlementDF$zone + max(coreFringeDF$zone)
#write.csv(settlementDF, 'settlementAreas.csv')

##-------------------- dump the three datasets ----------- ##
print("Saving found areas")
colOrder = c("zone" ,  "sqrkm" ,   "pop" ,     "x"  ,      "y"   ,     "areaType")
write.csv(rbind(coreAreaDF[,colOrder], fringeDF[,colOrder],coreFringeDF[,colOrder],
                settlementDF[,colOrder]),reName("peruAreas.csv",''))

## ------------------- Rasters into polygons ------------- ##
print("Saving areas as polygons: dissolved")
### this makes it easier to work with in GIS programs
coreFringesPoly = rasterToPolygons(coreFringes, dissolve=TRUE)
coreFringesPoly = coreFringesPoly[coreFringesPoly$layer>0,]
coreFringesPoly = disaggregate(coreFringesPoly)
#coreFringesPoly@data = data.frame(coreFringesPoly@data, coreFringeDF)
writeOGR(obj = coreFringesPoly, dsn = "coreFringe",layer="coreFringe",driver="ESRI Shapefile")

settlementPoly = rasterToPolygons(settlementAreas, dissolve=TRUE)
settlementPoly = settlementPoly[settlementPoly$layer>0,]
settlementPoly = disaggregate(settlementPoly)
#settlementPoly@data = data.frame(settlementPoly@data, settlementDF)
writeOGR(obj = settlementPoly, dsn = "settlements",layer="settlements",driver="ESRI Shapefile")

corePoly = rasterToPolygons(coreAreas,dissolve=TRUE)
corePoly = corePoly[corePoly$layer>0,]
corePoly = disaggregate(corePoly)
#corePoly@data = data.frame(corePoly@data, coreAreaDF)
writeOGR(obj = corePoly, dsn = "core",layer="core",driver="ESRI Shapefile")

##------------------- Layers for Pugo/Roca measure --------##
## sub in landscan data to the cores and dump to raster
print("Saving areas as polygons: not dissolved")
coreAreasVal = coreAreas
coreAreasVal[coreAreasVal != 1] = NA
coreAreasVal[coreAreasVal==1] = smthLyr[coreAreasVal==1]

tmpVal = rasterToPolygons(coreAreasVal, dissolve=FALSE)
writeOGR(obj = tmpVal, dsn = "coreLandScan",layer="coreAreas",driver="ESRI Shapefile",overwrite_layer=TRUE, delete_dsn=TRUE)

## repeat for core fringes
coreFringeVals = coreFringes
coreFringeVals[coreFringeVals != 1] = NA
coreFringeVals[coreFringeVals==1] = smthLyr[coreFringeVals==1]

tmpVal = rasterToPolygons(coreFringeVals, dissolve=FALSE)
writeOGR(obj = tmpVal, dsn = "coreFringeLandScan",layer="coreFringe",driver="ESRI Shapefile",overwrite_layer=TRUE, delete_dsn=TRUE)

## repeat for settlements
settlementAreaVals = settlementAreas
settlementAreaVals[settlementAreaVals != 1] = NA
settlementAreaVals[settlementAreaVals==1] = smthLyr[settlementAreaVals==1]

tmpVal = rasterToPolygons(settlementAreaVals, dissolve=FALSE)
writeOGR(obj = tmpVal, dsn = "settlementsLandScan",layer="settlements",driver="ESRI Shapefile",overwrite_layer=TRUE, delete_dsn=TRUE)
}}