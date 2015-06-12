# Working directories
root.dir = 'C:/Users/morueta/Documents/Documents_share/Projects'
env.data.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/', sep="/")
sp.data.dir = paste(root.dir, '101_TBC3_Modelling/Data/Ecoengine', sep="/")
out.dir = paste(root.dir, '101_TBC3_Modelling/ModelResults/Maxent/V1', sep="/")
bg.dir = paste(root.dir, '100_Postdoc/Data/Background_layers/PROCESSED', sep="/")

# Background maps
orig.project = '+proj=longlat +ellps=WGS84'
ta.project = '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'

# Libraries
require(raster)
require(rgdal)
require(fields)

# CC scenarios
allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')

# What species?
allSpecies = c("Sequoia sempervirens", "Quercus agrifolia", "Quercus garryana", 
               "Quercus douglasii", "Pseudotsuga menziesii")
for(i in 1:5) {
  #mySpecies = allSpecies[1] 
  mySpecies = allSpecies[i]



# What model run?
allModels = c("cwd", "cwd-tmn", "djf-jja-ppt-cwd")
mxModelType = allModels[3]

# Folder with projections for the species
proj.dir = paste(out.dir, mySpecies, mxModelType, "Projections", sep="/")
# Folder for figures of the species
fig.dir =  paste(out.dir, mySpecies, "Figures", sep="/")
if(file.exists(fig.dir) == F) {dir.create(fig.dir, recursive=F)} 

# # Folder with original maxent model
# mx.dir = paste(out.dir, mySpecies, mxModelType, sep="/")


#-------------------------------------
# Read in the original maxent model
mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))

#-------------------------------------
# Load California extent

bay = readOGR(bg.dir, "GADM_BayArea")
p.bay = spTransform(bay, CRS(ta.project))
pepperwood = readOGR(paste(root.dir, "100_Postdoc/Data/Background_layers/ORIGINAL/Pepperwood", sep="/"), "Pepperwood")

nbay= p.bay[p.bay$NAME_2 %in% c("Marin","Napa", "Sonoma"),]
bay=p.bay
rm(p.bay)

#-------------------
# Present and Future suitability Bay Area
present.bay = raster(paste(proj.dir, "Present_Bay.img", sep="/"))
futures.bay = stack(paste(proj.dir, "/", allScenarios, "_Suitability_Bay.img",sep=""))
occur = readRDS(paste(sp.data.dir,"/", mySpecies, ".rdata", sep=""))


#############
# Plot maps #
#############
# Set threshold for presence/absence limit - alternatively, filter with vegetation map

##### NOTE: CHANGE THRESHOLD FOR RAW OUTPUT!!
Threshold = 'Equal.training.sensitivity.and.specificity.logistic.threshold'
mx.th = as.numeric(mx@results[Threshold,])

# i = 1
# FS = futures.bay[[i]]

plotSimpleEdges = function(present = present.bay, future = FS, th, crop=F, cropto,YLAB="") {
  m = c(0,th,0,th,1,1)
  rclmat = matrix(m, ncol=3, byrow=T)
  pa.now = reclassify(present, rclmat)
  pa.fut = reclassify(future, rclmat)
  
  XLIM = extent(present)[1:2]
  YLIM = extent(present)[3:4]
  if(crop) {
    XLIM = extent(cropto)[1:2]
    YLIM = extent(cropto)[3:4]
  }
  pa.change = (2*pa.fut)-pa.now # multiply future by two to differentiate 0-0 and 1-1
  
#   plot(present, main="Present", col=topo.colors(100), zlim=c(0,1), axes=F,xlim=XLIM,ylim=YLIM,legend=F)
#   if(crop) {plot(cropto,add=T)}
#   plot(future, main="Future",  col=topo.colors(100), zlim=c(0,1), axes=F,xlim=XLIM,ylim=YLIM)
#   if(crop) {plot(cropto,add=T)}
  plot(pa.change, col=c("red","grey","black","blue"), 
       zlim=c(-1,2), axes=F, main="",xlim=XLIM,ylim=YLIM,legend=F)
  mtext(YLAB,2,cex=0.8)
  if(crop) {plot(cropto,add=T,border="orange",lwd=1.5)}
}

jpeg(paste(fig.dir, "/Suite_change_map_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1400,height=2100,quality=100,res=500,pointsize=8)
par(mfrow=c(6,2), mar=c(0,2.5,0,0))

for(i in 1:6) {
  plotSimpleEdges(mySpecies,present.bay,futures.bay[[i]],th=mx.th,crop=T,cropto=nbay,YLAB=allScenarios[i])
  plotSimpleEdges(mySpecies,present.bay,futures.bay[[i]],th=mx.th,crop=T,cropto=pepperwood)
  print(i)
}
dev.off()

# #Make legend as separate file
# # color scale
# plot(rep(1,4),seq(0,4,length=4),pch=15,col=c("red","grey","black","blue"),cex=6, axes=F,xlab=NA,ylab=NA,ylim=c(-1,6))
# text(rep(1.08,4),seq(0,4,length=4),c("Expansion","Same presence","Same absence","Contraction"),adj = c(0,0.5))
# title("Suitability")

###########################
# Change mean suitability #
###########################
# # Present and Future suitability Bay Area
# present.bay = raster(paste(proj.dir, "Present_Bay.img", sep="/"))
# futures.bay = stack(paste(proj.dir, "/", allScenarios, "_Suitability_Bay.img",sep=""))
# 

# present.processed.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/HST/Summary', sep="/")
# future.processed.dir = paste(root.dir, '100_Postdoc/Data/Climate/PROCESSED/Future/Summary', sep="/")
# 
# MAT.pres = raster(paste(present.processed.dir, "mat1981_2010_ave_HST_Bay.img",sep="/"))
# 
# allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
# fut.files = paste("mat2070_2099",allScenarios, "Bay",sep="_")
# MAT.fut = stack(paste(future.processed.dir, paste(fut.files,".img",sep=""), sep="/"))
# NB.MAT.pres = mask(MAT.pres,nbay)
# PW.MAT.pres = mask(MAT.pres,pepperwood)
# NB.MAT.fut = mask(MAT.fut, nbay)
# PW.MAT.fut = mask(MAT.fut, pepperwood)


# 
# # North Bay
# # nb.mat = cellStats(stack(NB.MAT.pres, NB.MAT.fut), "mean")
# nb.suit = cellStats(mask(stack(present.bay, futures.bay),nbay), "mean")
# 
# # Pepperwood
# # pw.mat = cellStats(stack(PW.MAT.pres, PW.MAT.fut), "mean")
# pw.suit = cellStats(mask(stack(present.bay, futures.bay),pepperwood), "mean")
# 
# l1 = c("Present",allScenarios)
# 
# jpeg(paste(fig.dir, "/Mean_suitability-Scenario_",mxModelType, ".jpg", sep=""),quality=100,width=600)
# par(mfrow=c(1,2))
# XLIM = range(c(nb.mat,pw.mat))+c(0,3)
# YLIM = range(c(nb.suit,pw.suit))
# plot(x=nb.mat,y=nb.suit, xlab="MAT", ylab="North Bay mean suitability",main=mySpecies,pch=19,col="red",xlim=XLIM,ylim=YLIM)
# text(labels=l1,x=nb.mat,y=nb.suit,cex=0.8, pos=4)
# 
# plot(x=pw.mat,y=pw.suit, xlab="MAT", ylab="Pepperwood mean suitability",main=mySpecies,pch=19,col="red",xlim=XLIM,ylim=YLIM)
# text(labels=l1,x=pw.mat,y=pw.suit,cex=0.8, pos=4)
# 
# dev.off()
}
