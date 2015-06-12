# Script for running maxent models on current climate data, project to
# future scenarios and map leading and trailing edges

# Working directories
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
env.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/', sep="/")
future.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/Future/Summary/', sep="/")

sp.data.dir = paste(root.dir, '101_TBC3_Modelling/Data/Ecoengine', sep="/")
out.dir = paste(root.dir, '101_TBC3_Modelling/ModelResults/Maxent/V1', sep="/")


# Projections
orig.project = '+proj=longlat +ellps=WGS84'
ta.project = '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'

# What species?
allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", "Quercus douglasii")

env.files <- list.files(path=env.data.dir, pattern='CA.img', full.names=FALSE)
env.files = env.files[-grep("xml",env.files)]


runMxModel = function(mySpecies, mxModelType, env.files, env.data.dir, sp.data.dir, out.dir) {
  #give java more ram - for maxent modeling later
  options(java.parameters = "-Xmx1g" )
  
  # Libraries
  require(dismo)
  require(rJava)
  require(raster) 
  
  # Which predictor variables?
  if (mxModelType == "4vars") {
    climnames = c("djf", "jja", "ppt", "cwd")
  } else if (mxModelType == "cwd") {
    climnames = "cwd"
  } else if (mxModelType == "cwd-tmn") {
    climnames = c("cwd","tmn")
  }
    
  # Read in predictors 
  files = env.files[which(substr(env.files,1,3)%in%climnames)]
  
  predictors <- stack(paste(env.data.dir,files,sep="/"))
  
  # Read in species occurrence data
  occur = readRDS(paste(sp.data.dir,"/", mySpecies, ".rdata", sep=""))
  
  # Directory to write files to
  mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")
  if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)} 
  
  # Arguments for maxent models
  mxArgs = c("-a", "-z", "outputformat=raw", "maximumbackground=30000", "nothreshold")
  
  # Run the model!
  mx <- maxent(predictors, occur, progress='text', path=mx.dir, args=mxArgs)
  
  saveRDS(mx, file = paste(mx.dir, 'ModelObject', sep='/'))
  
}

  
  
  

for (i in 1:length(allSpecies)) {
  i = 3
  mySpecies = allSpecies[i]

  ####################
  # Maxent modelling #
  ####################
  #give java more ram - for maxent modeling later
  options(java.parameters = "-Xmx1g" )
  
  # Libraries
  require(dismo)
  require(rJava)
  require(raster)
  
  # Which predictor variables?
  mxModelType = "cwd-tmn"
  
  if (mxModelType == "4vars") {
    climnames = c("djf", "jja", "ppt", "cwd")
  } else if (mxModelType == "cwd") {
    climnames = "cwd"
  } else if (mxModelType == "cwd-tmn") {
    climnames = c("cwd","tmn")
  }
  
  
  env.files <- list.files(path=env.data.dir, pattern='CA.img', full.names=FALSE)
  env.files = env.files[-grep("xml",env.files)]
  
  # Read in predictors 
  files = env.files[which(substr(env.files,1,3)%in%climnames)]
  
  predictors <- stack(paste(env.data.dir,files,sep="/"))
  
  #---------#
  # Read in species occurrence data
  occur = readRDS(paste(sp.data.dir,"/", mySpecies, ".rdata", sep=""))
  #---------#
  
  # Directory to write files to
  mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")
  if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)} 
  
  # Arguments for maxent models
  mxArgs = c("-a", "-z", "outputformat=raw", "maximumbackground=30000", "nothreshold")
  
  # Run the model!
  mx <- maxent(predictors, occur, progress='text', path=mx.dir, args=mxArgs)
  
  saveRDS(mx, file = paste(mx.dir, 'ModelObject', sep='/'))
  
  ######################################
  # Project to present and future data #
  ######################################
  # Load mx model
#  mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))

  # Outdir
  proj.dir = paste(out.dir, mySpecies, mxModelType, 'Projections', sep="/")
  if(file.exists(proj.dir) == F) {dir.create(proj.dir, recursive=T)} 
  
  ## Present
  
  # Load data CA
  pres.clim.ca = predictors
  
  # Load data Bay
  Bay.env.files = list.files(path=env.data.dir, pattern='Bay.img', full.names=FALSE)
  Bay.env.files = Bay.env.files[-grep("xml",Bay.env.files)]
  files = Bay.env.files[which(substr(Bay.env.files,1,3)%in%climnames)]
  pres.clim.bay = stack(paste(env.data.dir,files,sep="/"))
  names(pres.clim.bay) = names(pres.clim.ca) # Rename to match variables
  
  # Project model to present in CA and Bay
  present.ca = predict(mx, pres.clim.ca)
  present.bay = predict(mx, pres.clim.bay)
  
  # Save to disk
  writeRaster(present.ca, file=paste(proj.dir, "Present_CA.img", sep="/"), overwrite=T)
  writeRaster(present.bay, file=paste(proj.dir, "Present_Bay.img", sep="/"), overwrite=T)
  
  #------------#
  ## Future
  future.env.files <- list.files(path=future.data.dir, pattern='asc', full.names=FALSE)
  files = future.env.files[which(substr(future.env.files,1,3)%in%climnames)]
  
  allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
  
  # Which scenario?
  #i=1
  for(i in 1:length(allScenarios)) {
    subfiles = files[grep(allScenarios[i],files)]
    myScenario = allScenarios[i]
    
    # Load future Bay
    future.clim.bay <- stack(paste(future.data.dir,subfiles,sep="/"))
    names(future.clim.bay) = names(predictors) # Rename to match variables
    
    # Project model to future in Bay
    future.bay = predict(mx, future.clim.bay)
    
    # Save to disk
    writeRaster(future.bay, file=paste(proj.dir, "/", myScenario, "_Suitability_Bay.img", sep=""), overwrite=T)
  }

}

