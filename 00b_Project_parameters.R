## Working directories
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
env.data.dir = paste(root.dir, '100_Postdoc/Data/', sep="/")
pres.clim.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/', sep="/")
fut.clim.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/Future/Summary/', sep="/") 
sp.data.dir = 'Data/Ecoengine'
out.dir = 'ModelResults/Maxent/V1'
bg.dir = paste(root.dir, '100_Postdoc/Data/Background_layers/PROCESSED', sep="/")


## Working projections
# projection of Ecoengine species points
orig.project = '+proj=longlat +ellps=WGS84'
# equal area projection from Flints in which all analyses are done
ta.project = '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'

## Future climate change scenarios
# CC scenarios
allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')


## All species names for project
allSpecies = c("Pseudotsuga menziesii", "Adenostoma fasciculatum", "Sequoia sempervirens", 
               "Quercus agrifolia", "Quercus garryana", "Quercus douglasii", "Umbellularia californica")
vegIDs = list(17, c(9), 43, 10, 39, c(4,5), 6) 
names(vegIDs) = allSpecies
