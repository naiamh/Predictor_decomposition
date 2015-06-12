# Working directories
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
sp.data.dir = paste(root.dir, '101_TBC3_Modelling/Data/Ecoengine', sep="/")

# Project to the projection of environmental layers
orig.project = '+proj=longlat +ellps=WGS84'
ta.project = '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'

getOccur = function(mySpecies, db = "Ecoengine", georef = T, out.dir, in.project, out.project, save=F) {
  require(rgdal)
  require(raster)
  if(db == "Ecoengine") {
    require(ecoengine)
    require(plyr)
    
    species_pages <- ee_pages(ee_observations(scientific_name = mySpecies, progress = FALSE,quiet = TRUE, georeferenced=TRUE))
    page_breaks <- split(1:species_pages, ceiling(seq_along(1:species_pages)/1000))
    results <- ldply(page_breaks, function(x) {
      results <- ee_observations(scientific_name = mySpecies, page = x, quiet = TRUE, georeferenced=TRUE)
      results$data 
      }, .progress = "text")
    
    #results <- ee_observations(scientific_name = mySpecies, progress = FALSE, georeferenced = TRUE)
    occur = results[,c("longitude","latitude")]
  } else if (db == "BIEN2") {
    # extract from BIEN2
  }
  coordinates(occur) = ~ longitude + latitude
  projection(occur) = CRS(orig.project)
  p.occur = spTransform(occur, CRS(ta.project))
  
  if(save == T) {
    # Write occurrence file to disk
    saveRDS(p.occur, file = paste(sp.data.dir,"/", mySpecies, ".rdata", sep=""))
  } else {
    return(p.occur)
  }
}


spList = c("Quercus agrifolia", "Quercus garryana", "Quercus douglasii", "Pseudotsuga menziesii")
# Get species data
for(i in 1:length(spList)) {
  getOccur(mySpecies=spList[i], db="Ecoengine", out.dir=sp.data.dir, in.project=orig.project, out.project=ta.project, save=T)
  print(spList[i])
}

#getOccur(mySpecies="Pinus ponderosa",db="Ecoengine", out.dir=sp.data.dir, in.project=orig.project, out.project=ta.project, save=T)
#ee_map(results)






