## Scripts for swapping predictor layers and assess importance of each variable
## for changes in suitability

#---------------------#
# Working directories #
#---------------------#
setwd("C:/Users/morueta/Documents/Documents_share/Projects/")
env.data.dir = '100_Postdoc/Data/Climate/PROCESSED/HST/Summary/'
future.data.dir = '100_Postdoc/Data/Climate/PROCESSED/Future/Summary/'
model.dir = paste('101_TBC3_modelling/Lead-trail_R-project/ModelResults/Maxent/V2',sep="")

#-----------#
# Libraries #
#-----------#
require(dismo)
require(raster)

#----------------#
# Set parameters #
#----------------#
mySpecies = "Pinus ponderosa" # TRY WITH OTHER SPP. BLUE OAK, CA BAY, ETC.
mxModelType = "djf-jja-cwd"

allScenarios = c('MIROC_rcp85', 'CNRM_rcp85', 'GFDL_A2', 'CCSM4_rcp85', 'GFDL_B1', 'PCM_B1')
myScenario = allScenarios[1]

mx.dir = paste(model.dir, mySpecies, mxModelType, sep="/")

#---------#
# Action! #
#---------#

# Load mx model
mx = readRDS(paste(mx.dir, 'ModelObject', sep='/'))

# Layer combinations to project to
#
# # (MOVE TO 01_Prepare_env.layers file)
climnames = strsplit(mxModelType,"-")[[1]]
env.files = list.files(path=env.data.dir, pattern='Bay.img', full.names=FALSE)
env.files = env.files[-grep("xml",env.files)]
files = env.files[which(substr(env.files,1,3)%in%climnames)]

pres.clim = stack(paste(env.data.dir,files,sep="/"))
names(pres.clim) = sub("Bay","CA",names(pres.clim))

env.files = list.files(path=future.data.dir, pattern='Bay.asc', full.names=FALSE)
files = env.files[which(substr(env.files,1,3)%in%climnames)]
files = files[grep(myScenario,files)]
fut.clim = stack(paste(future.data.dir,files,sep="/"))
names(fut.clim) = names(pres.clim)
fut.clim[getValues(fut.clim[[1]])==-9999]=NA # --- QUICK FIX - ERASE FOR NEW LAYERS

rm(env.files,files)

# Make data.frames with present and future layers
mp = getValues(pres.clim)
mp = data.frame(mp[which(!is.na(mp[,1])),])
mf = getValues(fut.clim)
mf = data.frame(mf[which(!is.na(mf[,1])),])

# Suitability vectors for present and future layers
pres = predict(mx, mp)
fut = predict(mx, mf)

subcells = sample(nrow(mp), 1000)

# Check correlations among predictors
vars = c("cwd", "djf", "jja") #--------- WILL CHANGE FOR OTHER VARS
test=mp
names(test)=vars

pairs(~cwd+djf+jja,data=test[subcells,],col=rgb(0,0,0,0.1),pch=21)


# Stein-Alpert factor separation for 3 vars
mpf = list(mp[,1:3],mf[,1:3]) #only the first 3 vars
FAC = list()
FAC[[1]] = cbind(mpf[[1]][,1],mpf[[1]][,2],mpf[[1]][,3]) #v0 (all present values)
FAC[[2]] = cbind(mpf[[2]][,1],mpf[[1]][,2],mpf[[1]][,3]) #v1
FAC[[3]] = cbind(mpf[[1]][,1],mpf[[2]][,2],mpf[[1]][,3]) #v2
FAC[[4]] = cbind(mpf[[1]][,1],mpf[[1]][,2],mpf[[2]][,3]) #v3
FAC[[5]] = cbind(mpf[[2]][,1],mpf[[2]][,2],mpf[[1]][,3]) #v12
FAC[[6]] = cbind(mpf[[2]][,1],mpf[[1]][,2],mpf[[2]][,3]) #v13
FAC[[7]] = cbind(mpf[[1]][,1],mpf[[2]][,2],mpf[[2]][,3]) #v23
FAC[[8]] = cbind(mpf[[2]][,1],mpf[[2]][,2],mpf[[2]][,3]) #v123 all future values

for(i in 1:length(FAC)) {
  colnames(FAC[[i]]) = colnames(mpf[[1]])
}

preds = lapply(FAC, function(x) {predict(mx, x)})

# Create a stack of maps with the predictions
layer0 = pres.clim[[1]] #raster to fill in values
FS = stack(lapply(preds, function(x) {
  res=layer0
  res[!is.na(getValues(pres.clim[[1]]))]=x
  return(res)
}))

stackNames = c("v0", "v1", "v2", "v3", "v12","v13", "v23"  ,"v123")
names(FS) = stackNames

# Stein-Albert decomposition
DS=list()
DS[['d1']] = FS[['v1']] - FS[['v0']]
DS[['d2']] = FS[['v2']] - FS[['v0']]
DS[['d3']] = FS[['v3']] - FS[['v0']]
DS[['d12']] = FS[['v12']] - (FS[['v1']] + FS[['v2']]) + FS[['v0']]
DS[['d13']] = FS[['v13']] - (FS[['v1']] + FS[['v3']]) + FS[['v0']]
DS[['d23']] = FS[['v23']] - (FS[['v2']] + FS[['v3']]) + FS[['v0']]
DS[['d123']] = FS[['v123']] - (FS[['v12']] + FS[['v13']] + FS[['v23']]) + 
  (FS[['v1']] + FS[['v2']] + FS[['v3']]) - FS[['v0']]

DS = stack(DS)
dsf = getValues(DS)[which(!is.na(getValues(DS[[1]]))),]
subdsf = dsf[subcells,]
submp = mp[subcells,]
colnames(submp)=paste(vars,"present",sep="_")
submf = mf[subcells,]
colnames(submf)=paste(vars,"future",sep="_")
rownames(subdsf)=rownames(submp)=rownames(submf)=NULL

myData = cbind(submp,submf,subdsf)

#----------#
# Plotting #
#----------#

# Map of suitability changes for each variable and interaction term
maxrange = max(abs(c(min(getValues(DS),na.rm=T), max(getValues(DS),na.rm=T))))
ZLIM = c(-ceiling(maxrange*100)/100, ceiling(maxrange*100)/100)
plot(DS,zlim=ZLIM, col=colorRampPalette(c("red","grey","blue"))(10))

# Correlation between variable contributions?
pairs(~d1+d2+d3+d12+d13+d23+d123,data=subdsf,col=rgb(0,0,0,0.1))


# Plot maps of the difference in suitability for each combination relative to present

pdf("Plots/Test_maps_all_effects.pdf")
par(mfrow=c(3,3), mar=c(0.1,0.2,2,1),oma=c(0.1,0.1,0.1,1))
plot(FS[[1]],legend=T,main="Present",zlim=c(0,1),axes=F)
plot(FS[[8]],main="Future",zlim=c(0,1),axes=F,legend=F)
for(i in 1:nlayers(DS)) {
  if(i==1) {
    plot(DS[[i]], zlim=ZLIM, col=colorRampPalette(c("red","grey","blue"))(10), 
         legend=T, axes=F,main=names(DS)[i])
  } else {
  plot(DS[[i]], zlim=ZLIM, col=colorRampPalette(c("red","grey","blue"))(10),
       legend=F, axes=F,main=names(DS)[i])
  }
}
dev.off()

# Single variable effects
# Plot delta suitability ~ delta predictor over time (each point a pixel)
pdf("Plots/Test_single_var_effect.pdf")
par(mfrow=c(2,3),mar=c(4,4,0.5,0.5),oma=c(0,0,0,0))
for(i in 1:3) {
  var = vars[i]
  delta.var = submf[,i] - submp[,i]
  delta.suit = subdsf[,paste("d",i,sep="")]
  plot(delta.var, delta.suit,col=rgb(0,0,0,0.1),pch=21,
       ylab="delta.suitability",xlab=paste("Change in",var), ylim=ZLIM)
  abline(0,0,col="green")
}

# Plot corresponding maps of predicted change
for(i in 1:3) {
  var = vars[i]
  plot(DS[[paste("d",i,sep="")]], zlim=ZLIM, col=colorRampPalette(c("red","grey","blue"))(10),
       legend=T, axes=F)
}
dev.off()


deltavars = submf - submp
colnames(deltavars) = paste(vars,"delta",sep="_")

# Plot change in variable 1 and variable 2 with each point colored 
# by change in suitability based on their interaction

pdf("Plots/Test_plots_bivariate_v1.pdf")

  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  myx = "cwd"
  myy = "dfj"

  xv = myData[,paste(myx,"future",sep="_")]-myData[,paste(myx,"present",sep="_")]
  yv = myData[,paste(myy,"future",sep="_")]-myData[,paste(myy,"present",sep="_")]
  zv = myData[,"d12"]
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=range(zv))
  p1 = ggplot(myData, aes(x=xv, y=yv, colour=zv)) + 
            geom_point(size=3) + sc +
            xlab(paste("Delta",vars[1])) +
            ylab(paste("Delta",vars[2]))
  print(p1)
  
  xv = myData[,"cwd_future"]-myData[,"cwd_present"]
  yv = myData[,"jja_future"]-myData[,"jja_present"]
  zv = myData[,"d13"]
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=range(zv))
  p2 = ggplot(myData, aes(x=xv, y=yv, colour=zv)) + 
    geom_point(size=3) + sc +
    xlab(paste("Delta",vars[1])) +
    ylab(paste("Delta",vars[3]))
  print(p2)
  
  xv = myData[,"djf_future"]-myData[,"djf_present"]
  yv = myData[,"jja_future"]-myData[,"jja_present"]
  zv = myData[,"d23"]
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=range(zv))
  p3 = ggplot(myData, aes(x=xv, y=yv, colour=zv)) + 
    geom_point(size=3) + sc +
    xlab(paste("Delta",vars[2])) +
    ylab(paste("Delta",vars[3]))
  print(p3)


dev.off()


# Plot value of each predictor and change in predictor with each point
# colored by change in suitability
pdf("Plots/Test_plots1.pdf")
for(i in 1:3) {
  xv = myData[,paste(vars[i],"future",sep="_")]-myData[,paste(vars[i],"present",sep="_")]
  yv = myData[,paste(vars[i],"present",sep="_")]
  zv = myData[,paste("d",i,sep="")]
  
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-0.4, 0.4),
                               name = paste("d",i,sep=""))
  
  print(ggplot(myData, aes(x=xv, y=yv, colour=zv)) + 
    geom_point(size=3) + sc +
    xlab(paste("Delta",vars[i])) +
    ylab(paste("Present",vars[i])))
  
}
dev.off()



#----------------#
# Plotting - OLD #
#----------------#
subcells = sample(nrow(mp),10000)

# Plot full future layers set against set excluding change in each layer in turn
par(mfrow=c(2,2))
for(i in 1:4) {
  var = colnames(res)[i]
  plot(res[subcells,i],fut[subcells],col=rgb(0,0,0,0.1),pch=20,
       ylab="All predictors", xlab=paste("Without",var))
  abline(0,1,col="green")
}

# Plot delta suitability ~ delta predictor over time (each point a pixel)
par(mfrow=c(2,2),mar=c(5,5,0.5,0.5))
for(i in 1:4) {
  var = colnames(res)[i]
  delta.var = mf[,i] - mp[,i]
  delta.suit = res[,i] - pres
  plot(delta.var[subcells], delta.suit[subcells],col=rgb(0,0,0,0.1),pch=21,
       ylab="delta.suitability",xlab=paste("Change in",var))
  abline(0,0,col="green")
}

# Plot suitability in the present as a function of each variable
par(mfrow=c(2,2))
for(i in 1:4) {
  var = colnames(res)[i]
  pres.suit = pres
  myvar = mp[,i]
  plot(myvar[subcells], pres.suit[subcells],col=rgb(0,0,0,0.1),pch=21,
       ylab="present suitability",xlab=var)
}

# Plot vectors of change of suitability as a function of variable change
par(mfrow=c(2,2))
for(i in 1:4) {
  var = colnames(res)[i]
  pres.suit = pres[subcells]
  fut.suit = fut[subcells]
  pvar = mp[subcells,i]
  fvar = mf[subcells,i]
  plot(0,0,ylab="suitability",xlab=var,xlim=range(c(pvar,fvar),na.rm=T),ylim=c(0,0.9))
  segments(x0=pvar,x1=fvar,y0=pres.suit,y1=fut.suit,col=rgb(0,0,0,0.1))
  #points(pvar,pres.suit,col=rgb(0,0,1,0.1),pch=20,cex=0.1)
  #points(fvar,fut.suit,col=rgb(1,0,0,0.1),pch=20,cex=0.1)
}

