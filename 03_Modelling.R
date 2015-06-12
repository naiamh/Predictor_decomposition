# Script for running maxent models on current climate data, project to
# future scenarios and map leading and trailing edges

#---------------------#
# Working directories #
#---------------------#
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
env.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/', sep="/")
future.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/Future/Summary/', sep="/")

sp.data.dir = paste(root.dir, '101_TBC3_Modelling/Lead-trail_R-project/Data/Ecoengine', sep="/")
out.dir = paste(root.dir, '101_TBC3_Modelling/Lead-trail_R-project/ModelResults/Maxent/V2', sep="/")

#------------------#
# Source functions #
#------------------#
source(paste(root.dir, '101_TBC3_Modelling/Lead-trail_R-project/Scripts/00_Functions_trailing-edge.r', sep="/"))

#-------------------------#
# Parameters for modeling #
#-------------------------#
# Set the environmental file names (right now set to California layers)
env.files <- list.files(path=env.data.dir, pattern='CA.img', full.names=FALSE)
env.files = env.files[-grep("xml",env.files)]

# What species?
#allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", "Quercus douglasii")
#allSpecies = "Pinus ponderosa"
allSpecies = Species

#####################
# RUN MAXENT MODELS #
#####################

# Run the 4 variable models
for (i in 1:length(allSpecies)) {
  runMxModel(allSpecies[i], "djf-jja-ppt-cwd", env.files, env.data.dir, sp.data.dir, out.dir)
}

# Run a 3 variable model
for (i in 1:length(allSpecies)) {
  runMxModel(allSpecies[i], "djf-jja-cwd", env.files, env.data.dir, sp.data.dir, out.dir)
}


###################
# Project models  #
###################

# Which modeltype?
mxModelType = "djf-jja-cwd"
# Which predictor variables?
climnames = strsplit(mxModelType,"-")[[1]]


## Present
# California data
files = env.files[which(substr(env.files,1,3)%in%climnames)]
pres.clim.ca <- stack(paste(env.data.dir,files,sep="/"))

# Bay data
Bay.env.files = list.files(path=env.data.dir, pattern='Bay.img', full.names=FALSE)
Bay.env.files = Bay.env.files[-grep("xml",Bay.env.files)]
files = Bay.env.files[which(substr(Bay.env.files,1,3)%in%climnames)]
pres.clim.bay = stack(paste(env.data.dir,files,sep="/"))
names(pres.clim.bay) = names(pres.clim.ca) # Rename to match variables

# Which species?
#allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", "Quercus douglasii")

# Project to present layers
for (i in 1:length(allSpecies)) {
  mySpecies = allSpecies[i]
  # Load mx model
  mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")
  mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))
  
  # Outdir
  proj.dir = paste(out.dir, mySpecies, mxModelType, 'Projections', sep="/")
  if(file.exists(proj.dir) == F) {dir.create(proj.dir, recursive=T)} 

  #Project model to present in CA and Bay
  present.ca = predict(mx, pres.clim.ca, filename=paste(proj.dir, "Present_CA.img", sep="/"), overwrite=F)
  present.bay = predict(mx, pres.clim.bay, filename=paste(proj.dir, "Present_Bay.img", sep="/"), overwrite=F)
  
}

# --------------------------- #
## Future
# Which species?
#allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", "Quercus douglasii")

# Bay data
future.env.files <- list.files(path=future.data.dir, pattern='asc', full.names=FALSE)
files = future.env.files[which(substr(future.env.files,1,3)%in%climnames)]


# Which future climate scenarios?
#allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')

for (i in 1:length(allSpecies)) {
  mySpecies = allSpecies[i]
  # Load mx model
  mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")
  mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))
  
  proj.dir = paste(out.dir, mySpecies, mxModelType, 'Projections', sep="/")
  
  for(x in 1:length(allScenarios)) {
    subfiles = files[grep(allScenarios[x],files)]
    myScenario = allScenarios[x]
    
    # Load future Bay
    future.clim <- stack(paste(future.data.dir,subfiles,sep="/"))
    names(future.clim) = names(pres.clim.ca) # Rename to match variables
    
    # Project model to future in Bay
    future.bay = predict(mx, future.clim, file=paste(proj.dir, "/", myScenario, "_Suitability_Bay.img", sep=""), overwrite=T)
  }
}
   
    

