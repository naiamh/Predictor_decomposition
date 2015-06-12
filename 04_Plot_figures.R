# Working directories
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
env.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/', sep="/")
sp.data.dir = paste(root.dir, '101_TBC3_Modelling/Data/Ecoengine', sep="/")
out.dir = paste(root.dir, '101_TBC3_Modelling/ModelResults/Maxent/V1', sep="/")

# CC scenarios
allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')


# What species?
allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", "Quercus douglasii")
#mySpecies = allSpecies[3]
mySpecies = "Pseudotsuga menziesii"

# What model run?
allModels = c("cwd", "cwd-tmn", "djf-jja-ppt-cwd")
mxModelType = allModels[3]

# Folder with projections for the species
proj.dir = paste(out.dir, mySpecies, mxModelType, "Projections", sep="/")
# Folder for figures of the species
fig.dir =  paste(out.dir, mySpecies, "Figures", sep="/")
if(file.exists(fig.dir) == F) {dir.create(fig.dir, recursive=F)} 

# Folder with original maxent model
mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")

# Libraries
require(raster)



####################################
# Suitability against climate vars #
####################################
# Read in the original maxent model
mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))


# Suitability values from the present projections
present.ca = raster(paste(proj.dir, "Present_CA.img", sep="/"))
CA.suitability = values(present.ca)

present.bay = raster(paste(proj.dir, "Present_Bay.img", sep="/")) 
Bay.suitability = values(present.bay)

occur = readRDS(paste(sp.data.dir,"/", mySpecies, ".rdata", sep=""))

# #------- GET CLIMATE LAYERS -----------#
# 
# # Get values from original climate layers
# getClimate = function(vars=c("djf", "jja", "ppt", "cwd"), env.data.dir, domain) {
#   require(raster)
#   if(domain=='CA') {
#     mypattern='CA.img'
#   } else if (domain=='Bay') {
#     mypattern='Bay.img'
#   } 
#   env.files = list.files(path=env.data.dir, pattern=mypattern, full.names=FALSE)
#   env.files = env.files[-grep("xml",env.files)]
#   files = env.files[which(substr(env.files,1,3)%in%vars)]
#   clim <- stack(paste(env.data.dir,files,sep="/"))
#   
#   climate = values(clim)
#   climate[which(climate==-9999)] = NA
#   colnames(climate) = sort(vars)
#   return(climate)
# }
# 
# CA.climate = getClimate(env.data.dir=env.data.dir, domain='CA')
# Bay.climate = getClimate(env.data.dir=env.data.dir, domain='Bay')


# Convert occurrence points to p/a vector for CA
CA.pa = rep(0,nrow(CA.climate))
CA.pa[cellFromXY(present.ca,occur)] = 1

Bay.pa = rep(0,nrow(Bay.climate))
Bay.pa[cellFromXY(present.bay,occur)] = 1

Threshold = 'Equal.training.sensitivity.and.specificity.logistic.threshold'
mx.th = as.numeric(mx@results[Threshold,])

# # With David's function
#plotThresh(C=CA.climate,obs=CA.pa,typ='obs',fac=c('cwd','djf'),rs=1e4)
# plotThresh(C=Bay.climate,obs=Bay.pa,typ='obs',vt=1,fac=c('cwd','djf'),rs=1e4)
# 
#plotThresh(C=CA.climate,pred=CA.suitability,typ='mod',fac=c('cwd','djf'),th=mx.th,rs=1e4)
# plotThresh(C=Bay.climate,pred=Bay.suitability,typ='mod',fac=c('cwd','djf'),th=0.25,rs=1e4)
# 



jpeg(paste(fig.dir, "/Clim_space_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1400,height=2100,quality=100,res=500,pointsize=8) 
plotPanels(yvars = c("ppt","djf","jja"), xvar = "cwd",climate=CA.climate, pa=CA.pa, suitability=CA.suitability, threshold=mx.th, 
           mySpecies=mySpecies, addsub=T, subclimate=Bay.climate, subpa=Bay.pa, 
           subsuitability=Bay.suitability)
dev.off()


# # #Make legend as separate file
# # color scale
# require(fields)
# plot(rep(1,100),seq(0,20,length=100),pch=15,col=tim.colors(100),cex=1.8, axes=F,xlab=NA,ylab=NA)
# text(rep(1.05,11),seq(0,20,length=11),as.character(seq(0,1,0.1)),adj = c(0,0.5))
# title("Suitability")





