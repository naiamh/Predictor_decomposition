# Script for cropping environmental layers to be used in modelling

#---------------------#
# Working directories #
#---------------------#
Computer = 'HP'
if(Computer == 'HP') {
  root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects/100_Postdoc/Data'
  present.orig.dir = paste(root.dir, 'Climate/ORIGINAL/HST/Summary/', sep="/")
  future.orig.dir = paste(root.dir, 'Climate/ORIGINAL/CSIRO_A1B/', sep="/")
  
    
} else if (Computer == 'Ack-lab') {
  root.dir = 'C:/Naia/Data'
  future.orig.dir = 'G:/'
}
bg.dir = paste(root.dir, 'Background_layers/PROCESSED', sep="/")
# Out folders
present.processed.dir = paste(root.dir, 'Climate/PROCESSED/HST/Summary/', sep="/")
future.processed.dir = paste(root.dir, 'Climate/PROCESSED/Future/Summary/', sep="/")

# Future climate scenarios
myScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
#climnames = c('aet', 'aprpck', 'cwd', 'pet', 'ppt', 'rch', 'run', 'tmn', 'tmx')
climnames = c('djf','jja')
years = '2070_2099'

myScenarioFiles = sapply(myScenarios, function(x) {paste(climnames, years, '_ave_',x, '.asc', sep='')})

#-----------#
# Libraries #
#-----------#
require(raster)
require(rgdal)

#-------------#
# Projections #
#-------------#
orig.project = '+proj=longlat +ellps=WGS84'
ta.project = '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'

#--------------#
# Process data #
#--------------#

## Load California extent
ca = readOGR(bg.dir, "GADM_California")
bay = readOGR(bg.dir, "GADM_BayArea")
p.ca = spTransform(ca, CRS(ta.project))
p.bay = spTransform(bay, CRS(ta.project))

# Read in predictors

if(Computer == 'HP') {
  present.env.files = list.files(path=present.orig.dir, pattern='asc', full.names=TRUE)
  future.env.files = list.files(path=future.orig.dir, pattern='asc', full.names=TRUE)
}

crop.fn = function(files, mask.layer) {
  myStack = stack(files)
  myCrop = crop(myStack, extent(mask.layer))
  res = mask(myCrop, mask=mask.layer)
  return(res)
} 

saveStack = function(d1, outdir, suffix, overwrite=F) {
  d2 <- unstack(d1)
  outputnames <- paste(names(d1), suffix,sep="")
  for(i in seq_along(d2)){writeRaster(d2[[i]], file=paste(outdir, outputnames[i], sep="/"),overwrite=overwrite)}
}


pres.clim.ca = crop.fn(present.env.files, p.ca)
saveStack(pres.clim.ca, outdir=present.processed.dir, suffix="_CA.asc")
rm(pres.clim.ca)

pres.clim.bay = crop.fn(present.env.files, p.bay)
saveStack(pres.clim.bay, outdir=present.processed.dir, suffix="_Bay.asc")
rm(pres.clim.bay)

future.clim.ca = crop.fn(future.env.files[7:9], p.ca)
saveStack(future.clim.ca, outdir=future.processed.dir, suffix="_CA.asc")
rm(future.clim.ca)

future.clim.bay = crop.fn(future.env.files, p.bay)
saveStack(future.clim.bay, outdir=future.processed.dir, suffix="_Bay.asc")
rm(future.clim.bay)

### For Ackerly lab computer
for(i in 2:ncol(myScenarioFiles)) {
  future.path = paste(future.orig.dir, myScenarios[i], '/',myScenarioFiles[,i], sep='')
  fca = crop.fn(future.path, p.bay)
  saveStack(fca, outdir=future.processed.dir, suffix="_Bay.asc")
  print(i)
}

# Special processing of jja and djf layers download from Cal Climate
future.orig.dir = 'C:/Naia/Data/Climate/ORIGINAL/Future/Summary'
# Future climate scenarios
myScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
#climnames = c('aet', 'aprpck', 'cwd', 'pet', 'ppt', 'rch', 'run', 'tmn', 'tmx')
climnames = c('djf','jja')
years = '2070_2099'

myScenarioFiles = dir(future.orig.dir)[-1]
new.names = sort(c(paste('djf',years,'_ave_',myScenarios,sep=''), paste('jja',years,'_ave_',myScenarios,sep='')))

for(i in 1:length(myScenarioFiles)) {
  future.clim = paste(future.orig.dir, myScenarioFiles[i], sep='/')
  r = raster(future.clim)
  fca = crop.fn(r, p.bay)
  names(fca) = new.names[i]
  saveStack(fca, outdir=future.processed.dir, suffix="_Bay.asc")
  print(i)  
}

# Special scripts for projecting jja and djf files to right resolution
djf.files = list.files(path=future.processed.dir, pattern='djf', full.names=TRUE)
jja.files = list.files(path=future.processed.dir, pattern='jja', full.names=TRUE)
my.files = c(djf.files, jja.files)
myStack = stack(my.files)
r = raster(list.files(path=future.processed.dir, pattern='asc', full.names=TRUE)[1])

s = resample(myStack, r, method="bilinear")

# Write to disk - wrong files erased - redo all with matching layers
saveStack(s, outdir=future.processed.dir,suffix=".asc")

## Special scripts for changing wrong 0's to NA in cwd

# Bay
files = list.files(path=present.processed.dir, pattern='Bay.img', full.names=TRUE)
b.clim = stack(files[-grep("xml",files)])
cwd = all.clim[[which(names(all.clim)=='cwd1981_2010_ave_HST_Bay')]]
# Set all to NA where there are -9999 values or where cwd=0
all.clim[which(values(cwd)==-9999)] = NA
all.clim[which(values(cwd)==0)] = NA
# Write to disk - wrong files erased - redo all with matching layers
##saveStack(all.clim, outdir=present.processed.dir,suffix=".img",overwrite=T)

# CA
files = list.files(path=present.processed.dir, pattern='CA.img', full.names=TRUE)
all.clim = stack(files[-grep("xml",files)])
cwd = all.clim[[which(names(all.clim)=='cwd1981_2010_ave_HST_CA')]]
# Set all to NA where there are -9999 values or where cwd=0
all.clim[which(values(cwd)==-9999)] = NA
all.clim[which(values(cwd)==0)] = NA
# Write to disk - wrong files erased - redo all with matching layers
##saveStack(all.clim, outdir=present.processed.dir,suffix=".img",overwrite=T)

#---------------------------
# Special scripts for calculating MAT for each time period
climnames = c("tmn","tmx")
# Load climate layers
# Present
files <- list.files(path=present.processed.dir, pattern='Bay.img', full.names=FALSE)
files = files[which(substr(files,1,3)%in%climnames)]
files = files[-grep("xml",files)]
present = stack(paste(present.processed.dir, files,sep="/"))
present[getValues(present)==-9999] = NA
present = mean(present,na.rm=T)
names(present) = "mat1981_2010_ave_HST_Bay"
writeRaster(present,file=paste(present.processed.dir, paste(names(present),".img",sep=""), sep="/"),overwrite=F)






# Future
files = list.files(path=future.processed.dir, pattern='Bay.asc', full.names=FALSE)
files = files[which(substr(files,1,3)%in%climnames)]

allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
futures = lapply(allScenarios, function(x) {
  maxmin = stack(paste(future.processed.dir, files[grep(x,files)], sep="/"))
  maxmin[getValues(maxmin[[1]])==-9999] = NA
  MAT = mean(maxmin,na.rm=T)
  return(MAT)  
})
names(futures) = paste("mat2070_2099",allScenarios, "Bay",sep="_")
for (i in 1:length(futures)) {
  writeRaster(futures[[i]],file=paste(future.processed.dir, paste(names(futures)[[i]],".img",sep=""), sep="/"),overwrite=F)
}

MAT.fut = sapply(futures,function(x) {cellStats(x,na.rm=T,"mean")})


#-------------------
# Read in CLN vegetation map and project to Flint projection
veg.orig = raster(paste(env.data.dir, "Vegetation/ORIGINAL/cln_veg10",sep="/"))
veg.p = projectRaster(veg.orig,crs=ta.project)
writeRaster(veg.p, file=paste(env.data.dir, "Vegetation/PROCESSED/cln_veg10_p.img",sep=""))

# Read in landscape units and project to Flint projection
lu = readOGR(paste(env.data.dir, "Background_layers/ORIGINAL/tbc3_landscape_units",sep="/"), "tbc3_landscape_units")
lu.p = spTransform(lu, CRS(ta.project))
saveRDS(lu.p,paste(env.data.dir,"Background_layers/PROCESSED/tbc3_landscape_units_p.rdata",sep="/"))

