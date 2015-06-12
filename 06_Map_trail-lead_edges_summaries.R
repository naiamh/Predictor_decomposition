
# Scripts for plotting suitability maps, vegetation, future projections
# and changes (leading/trailing edged) as well as summary plots for landscape units

# Clear workspace
rm(list=ls())

# Source project parameters and functions
source("Scripts/00b_Project_parameters.r")
source("Scripts/00_Functions_trailing-edge.r")

# Libraries
require(raster)
require(rgdal)
require(fields)

### Data
# Vegetation map
veg = raster(paste(env.data.dir, "Vegetation/PROCESSED/cln_veg10_p.img",sep="/"))
veg.names = read.csv(paste(env.data.dir, "Vegetation/ORIGINAL/veg10_type.csv",sep="/"))

# Landscape units
lu = readRDS(paste(env.data.dir, "Background_layers/PROCESSED/tbc3_landscape_units_p.rdata",sep="/"))

#---------------#
# Set variables #
#---------------#
  
  
  # Species
  #mySpecies = "Pseudotsuga menziesii"
  for(sp in 1:length(allSpecies)) {
   mySpecies = allSpecies[sp]


  # Which model type
  mxModelType = "djf-jja-ppt-cwd"
  # Vegtype
  ids = vegIDs[[which(names(vegIDs)==mySpecies)]]
  # Landscape unit name
  luName = "North East Bay Hills"
  luAcr = "NEBH"
  

#-----------------------#
# Read in specific data #
#-----------------------#
# Species specific data and outputs:
source("Scripts/06_0_Load_or_compute_species_files.r", echo=T)

# Landscape unit specific variables
  slu = subset(lu, Name==luName)
  subveg = mask(crop(veg, slu),slu)
  vpoints = xyFromCell(subveg, which(round(getValues(subveg))%in%ids),spatial=T) 
      #<----- OBSOBS! Error in veg projections causing decimals needs to be fixed!!
  
  # Present suitability within landscape unit
  ps = mask(crop(present.bay, slu),slu)
  projection(ps) = ta.project
  
  # Future suitability within landscape unit for all CC scenarios
  fs = mask(crop(futures.bay, slu),slu)
  projection(fs) = ta.project
  
  # Presence/absence threshold from maxent models
  Threshold = 'Equal.training.sensitivity.and.specificity.logistic.threshold'
  mx.th = as.numeric(mx@results[Threshold,])

#---------#
# ACTION! #
#---------#
gmask = getMask(bg=ps,region=slu, vegmap=subveg, ids=ids) 

# # Plot current obs. veg, present suitability and 6 future scenarios
# jpeg(paste(fig.dir, "/",luAcr,"_Projections_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1100,height=2200,quality=100,res=500,pointsize=8)
# plotPresFut(ps,fs,mx.th,subveg,gmask)
# dev.off()
# 
# ## Plot change in modelled suitability
# jpeg(paste(fig.dir, "/Bay_Changes_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1100,height=2200,quality=100,res=500,pointsize=8)
# plotTrailLead(ps, fs, mx.th, subveg, slu, gmask)
# dev.off()
# 
# 
## Plot summary across landscape unit
jpeg(paste(fig.dir, "/",luAcr, "_Summary_change_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1000,height=500,quality=100,res=500,pointsize=8)
par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(0,0,0,0))

oldCols = plotSummarySq(ps,fs,main="Landscape unit",giveCols=T)

#----plotSummarySq with masked ps and fs
newCols=getSubCols(ps=ps,fs=fs,gmask=gmask,mx.th=mx.th)
plot(-1,xlim=c(0.5,3.5),ylim=c(0.5,3.5),axes=F,xlab="",ylab="")
rect(xleft=c(1,2,1,2),ybottom=c(2,2,1,1),xright=c(2,3,2,3),ytop=c(3,3,2,2), col=newCols)

dev.off()

# Get the % area in each category of change
# PA map for all LU
PA = getOrPlotChange(present=ps, future=fs, th=mx.th, crop=F,Plot=F) 

vals = c(1, 0, 2, -1)
names(vals) = c("Same pres", "Same abs", "Lead", "Trail")

res = matrix(NA, nrow=nlayers(PA), ncol=4)
colnames(res) = names(vals)
for (i in 1:nlayers(PA)) {
  res[i,] = sapply(vals, function(x) {length(which(getValues(PA[[i]])==x))})  
}
res=as.data.frame(res)
res$Scenario = names(fs)
res$Class = c("hot-dry","hot-wet","hot-dry","hot-wet","warm-dry","warm-wet")


# % area change for each of 6 scenarios across all LU
change.m = t(apply(res[,1:4], 1, function(x) {x/sum(x)*100}))
mycolors = c("black","white","blue","red")

# par(mfrow=c(2,3),mar=c(0.2,4,2,0))
# for (i in 1:nlayers(PA)) {
#   barplot(change.m[i,], main=allScenarios[[i]], col=mycolors, ylim=c(0,ceiling(max(change.m))),xaxt='n',ylab="% area")
# }


# % area change for each of 4 scenario classes across all LU 
res2 = aggregate(res[,1:4], by=list(res$Class), FUN=sum)
change.m2 = t(apply(res2[,2:5],1, function(x) {x/sum(x)*100}))
scen.class = c("warm-wet", "hot-wet", "warm-dry","hot-dry")
# order by scenario class (warm-wet, hot-wet, warm-dry, hot-dry) for plotting
change.m2 = change.m2[match(scen.class,res2[,1]),]


jpeg(paste(fig.dir, "/",luAcr, "_BarPlotChange_LU_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1000,height=1000,quality=100,res=500,pointsize=8)
par(mfrow=c(2,2),mar=c(0.2,3,2,0), oma=c(0.1,0.1,0.1,0.1))
for (i in 1:4) {
  barplot(change.m2[i,], main=scen.class[i], col=mycolors, ylim=c(0,ceiling(max(change.m2))),xaxt='n',ylab='')
  mtext("% area",2,line=2,cex=0.8)
}
dev.off()

}





# 
# 
# #----------------# OLD STUFF - CLEAN UP
# # Summaries suitabilities in veg type of landscape unit
# psuit = mean(getValues(ps)[unique(cellFromXY(ps,vpoints))],na.rm=T)
# fsuits = c()
# for(i in 1:nlayers(fs)) {
#   res = mean(getValues(fs[[i]])[unique(cellFromXY(fs[[i]],vpoints))],na.rm=T)
#   fsuits=c(fsuits,res)
# }
# 
# # Summary suitabilities across whole landscape unit (maybe at least exclude human sites)
# psuit_all = mean(getValues(ps),na.rm=T)
# fsuits_all = c()
# for(i in 1:nlayers(fs)) {
#   res = mean(getValues(fs[[i]]),na.rm=T)
#   fsuits_all=c(fsuits_all,res)
# }
# 
# #Define scenarios to lump:
# wawet = "PCM_B1"
# howet = c("CCSM4_rcp85","CNRM_rcp85")
# wadry = "GFDL_B1"
# hodry = c("GFDL_A2","MIROC_rcp85")
# Gr = list(wawet,howet,wadry,hodry)
# 
# AllScenarios = gsub("_Suitability_Bay","",names(fs))
# 
# Means = sapply(Gr,function(x) {mean(fsuits_all[which(AllScenarios%in%x)])})
# Means_v = sapply(Gr,function(x) {mean(fsuits[which(AllScenarios%in%x)])})
# 
# 
# 
# 
# # Plot colored summaries of change
# jpeg(paste(fig.dir, "/",luAcr, "_Summary_change_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1000,height=550,quality=100,res=500,pointsize=8)
# par(mfrow=c(1,2),mar=c(0,0.5,1,0),oma=c(0,0,0,0))
# 
# myCols = sapply(Means,function(x) {pickCol(psuit=psuit,fsuit=x,mx.th=mx.th)})
# plot(x=c(1,2,1,2),y=c(2,2,1,1), pch=22,cex=6,xlim=c(0,3),ylim=c(0,3),axes=F,bg=myCols,main="Landscape unit",cex.main=0.8)
# text(x=c(1,2),y=c(2.5,2.5),c("Warm","Hot"),cex=0.6,pos=3)
# text(x=c(0.3,0.3),y=c(1,2),c("Dry","Wet"),cex=0.6,srt=90)
# 
# myCols = sapply(Means_v,function(x) {pickCol(psuit=psuit,fsuit=x,mx.th=mx.th)})
# plot(x=c(1,2,1,2),y=c(2,2,1,1), pch=22,cex=6,xlim=c(0,3),ylim=c(0,3),axes=F,bg=myCols,main="Within veg. type",cex.main=0.8)
# text(x=c(1,2),y=c(2.5,2.5),c("Warm","Hot"),cex=0.6,pos=3)
# text(x=c(0.3,0.3),y=c(1,2),c("Dry","Wet"),cex=0.6,srt=90)
# dev.off()
# 
# #######
# ## Leading - trailing edges for all Bay Area
# subveg = veg
# ps = present.bay
# fs = futures.bay
# ids=c(4,5)
# vpoints = xyFromCell(subveg, which(getValues(subveg)%in%ids),spatial=T)
# slu = lu
# 
# jpeg(paste(fig.dir, "/Bay_Projections_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1100,height=2200,quality=100,res=500,pointsize=8)
# par(mfrow=c(4,2),mar=c(0.1,0.1,1.5,0))
# plot(subveg,axes=F,legend=F)
# title("Obs. vegetation")
# points(vpoints,pch=0,cex=0.1)
# greyOutTh(ps, mx.th, legend=T)
# title("Present suitability",cex=0.8)
# for(i in 1:nlayers(fs)){
#   greyOutTh(fs[[i]],mx.th)
#   scen=strsplit(names(fs[[i]]),"_")[[1]][1:2]
#   title(paste(scen[1],scen[2],sep="_"),cex=0.8)  
# }
# dev.off()
# 
# # Cells with ok vegetation visible
# vc = unique(cellFromXY(ps,vpoints))
# gmask = ps
# gmask[] = 1
# gmask[vc] = NA
# 
# 
# 
# 
# jpeg(paste(fig.dir, "/Bay_Changes_",mxModelType, "_t-",mx.th, ".jpg", sep=""),width=1100,height=2200,quality=100,res=500,pointsize=8)
# par(mfrow=c(4,2),mar=c(0.1,0.1,1.5,0))
# plot(subveg,axes=F,legend=F)
# plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.9))
# #plot(slu,add=T)
# title("Obs. vegetation")
# greyOutTh(ps, mx.th, legend=T)
# #plot(slu,add=T)
# plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.5))
# title("Present suitability",cex=0.8)
# for(i in 1:nlayers(fs)){
#   plotSimpleEdges(ps,fs[[i]],mx.th)
#   plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.5))
#   #plot(slu,add=T)
#   if(i == 1) {
#     x0 = xmin(fs[[i]])+(xmax(fs[[i]])-xmin(fs[[i]]))/100
#     xt = x0 + (xmax(fs[[i]])-xmin(fs[[i]]))/50
#     y0 = ymin(fs[[i]])+(ymax(fs[[i]])-ymin(fs[[i]]))/10
#     y1 = ymin(fs[[i]])+(ymax(fs[[i]])-ymin(fs[[i]]))/2
#     trail.legend(x0=x0,xt=xt,y0=y0,y1=y1)
#   }
#   scen=strsplit(names(fs[[i]]),"_")[[1]][1:2]
#   title(paste("Change",scen[1],scen[2],sep="_"),cex=0.8)  
# }
# dev.off()
# 
