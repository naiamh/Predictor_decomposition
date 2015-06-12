# Functions for the trailing-leading edge project

#----------------------#
# Extract species data #
#----------------------#

# Function to extract occurrence data (currently only from Ecoengine)
# The data can be optionally be written to disk. Input and output projection
# desired must be given
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
    saveRDS(p.occur, file = paste(out.dir,"/", mySpecies, ".rdata", sep=""))
  } else {
    return(p.occur)
  }
}

#----------------------#
# Extract climate data #
#----------------------#
## Function to load in either historic or future climate data for "CA" or "Bay" region
# and future scenario "scen" for specified variables (by mxModelType). 
# Returns a stack of climate layers
getClim = function(env.dir, region, period, mxModelType, scen) {
  require(raster)
  climnames = strsplit(mxModelType,"-")[[1]]
  if(period == "HST") {
    env.files = list.files(path=env.dir, pattern=paste(region,'.img',sep=""), full.names=FALSE)
    env.files = env.files[-grep("xml",env.files)]
    files = env.files[which(substr(env.files,1,3)%in%climnames)]
    clim = stack(paste(env.dir,files,sep="/"))
    names(clim) = sort(climnames) # fix to match maxent var names
    
    rm(env.files,files)
  } else if(period == "Future") {
    env.files = list.files(path=env.dir, pattern='asc', full.names=FALSE)
    files = env.files[which(substr(env.files,1,3)%in%climnames)]
    subfiles = files[grep(scen,files)]
    clim = stack(paste(env.dir,subfiles,sep="/"))
    names(clim) = sort(climnames) # fix to match maxent var names
  }
  return(clim)
}


# Function for running a maxent model for a single species. Output format now is set to raw,
# with 30,000 background points chosen and no threshold features
runMxModel = function(mySpecies, mxModelType, env.files, env.data.dir, sp.data.dir, out.dir) {
  #give java more ram - for maxent modeling later
  options(java.parameters = "-Xmx1g" )
  
  # Libraries
  require(dismo)
  require(rJava)
  require(raster) 
  
  # Which predictor variables?
  climnames = strsplit(mxModelType,"-")[[1]]
  
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

#----------#
# Plotting #
#----------#

# David's function for plotting thresholds adjusted
plotThresh <- function(C,pred=fit,obs=obsveg,typ='obs',fac=c('cwd','jja'),th=0.5,rs=1e4)
{
  require(geometry)
  require(classInt)
  require(fields)
  allpts <- C[,fac]
  allpts <- allpts[complete.cases(allpts),]
  yr <- range(allpts[,2])
  xr <- range(allpts[,1])
  rSamp <- sample(nrow(allpts),rs)
  plot(allpts[rSamp,],pch='.',xlim=xr,ylim=yr,col="darkgrey") 
  ch <- convhulln(allpts)
  for (i in 1:nrow(ch)) lines(allpts[ch[i,],],col='black')
  if (typ=='obs') {
    vs <- which(obs==1)
    pts <- C[vs,fac]
    pts <- pts[complete.cases(pts),]
    points(pts,col='red',pch=19,cex=0.1)
    ch <- convhulln(pts)
    for (i in 1:nrow(ch)) lines(pts[ch[i,],],col='red')
  } else if (typ=='mod') {
    myData <- cbind(C[,fac], pred)
    myData <- as.data.frame(myData[complete.cases(myData),])
    
    #plotclr = colorRampPalette(c("blue","red"),space="Lab")(100)
    plotclr = tim.colors(100)
    cuts = classIntervals(myData[,'pred'], style="fixed", fixedBreaks=seq(0,1,0.01))
    myData$colcode <- findColours(cuts, plotclr)
    
    vs <- which(myData[,'pred']>=th)
    pts <- myData[vs,]
    rownames(pts) <- NULL
    rSamp <- sample(nrow(pts),rs)
    points(pts[rSamp,fac],col=pts[rSamp,'colcode'],pch=19,cex=0.1) 
    ch <- convhulln(pts[,fac])
    for (i in 1:nrow(ch)) lines(pts[ch[i,],fac],col='red')
  }
}

addPlotTh = function(C,fac=c('cwd','jja'))
{
  require(geometry)
  allpts <- C[,fac]
  allpts <- allpts[complete.cases(allpts),]
  yr <- range(allpts[,2])
  xr <- range(allpts[,1])
  ch <- convhulln(allpts)
  for (i in 1:nrow(ch)) lines(allpts[ch[i,],],col='black',lty="dashed")
}

plotPanels = function(yvars = c("ppt","djf","jja"), xvar = "cwd", climate, pa, suitability, threshold, mySpecies=mySpecies, subclimate, subpa, subsuitability, addsub=F) {
  par(mfrow=c(length(yvars),2),mar=c(2,2,0,0), oma=c(2,2,3,1))
  for(y in 1:length(yvars)) {
    plotThresh(C=climate,obs=pa,typ='obs',fac=c(xvar,yvars[y]),rs=1e4)
    axis(2,ylab=yvars[y])
    mtext(yvars[y], side = 2, line= 2.5,las = 3)
    if(y==1) {
      mtext(paste(mySpecies,"- obs"), side = 3, line=1.5, cex=0.8)
    }
    if(y==length(yvars)) {
      mtext(xvar, side = 1, line= 2.5, cex=0.8)
    }
    if(addsub==T) {
      addPlotTh(C=subclimate,fac=c(xvar,yvars[y]))
    }
    
    plotThresh(C=climate,pred=suitability,typ='mod',fac=c(xvar,yvars[y]),th=threshold,rs=1e4)
    if(y==1) {
      mtext(paste("Pred (threshold = ",threshold, ")",sep=''), side = 3, line=1.5, cex=0.8)
    }
    if(y==length(yvars)) {
      mtext(xvar, side = 1, line= 2.5, cex=0.8)
    }
    if(addsub==T) {
      addPlotTh(C=subclimate,fac=c(xvar,yvars[y]))
    }    
  }
}

# Function to plot change in mean suitability of a species in a region as a function of
# change in MAT (as proxy of climate change scenario)
plotChangeMeanSuit = function(MAT.pres, MAT.fut, suit.pres, suit.fut) {
  require(raster)
  mat0 = cellStats(MAT.pres,"mean")
  mat1 = cellStats(MAT.fut, "mean")
  suit0 = cellStats(suit.pres, "mean")
  suit1 = cellStats(suit.fut, "mean")
  points(x=c(mat0,mat1), y=c(suit0,suit1), xlab="Change MAT", ylab="Change suitability")
}

###### Trailing edge functions
# Function to calculate temporal change in suitability. Input is the present
# rasterof suitability (present), the future raster of suitability, maxent 
# threshold (th), whether it should return the difference or plot it (Plot)
# whether the map should be cropped (crop). If crop=T, supply cropping
# raster or shapefile (cropto). The ylab can be set (YLAB).
# NOTE: leading/trailing edges and no-change areas are based on whether suitability
# is above or below the presence/absence threhold
getOrPlotChange = function(present, future, th, crop=F, cropto,YLAB="",Plot=T) {
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
  
  if(Plot){
    plot(pa.change, col=c("red","white","black","blue"), 
         zlim=c(-1,2), axes=F, main="",xlim=XLIM,ylim=YLIM,legend=F)
    mtext(YLAB,2,cex=0.8)
    if(crop) {plot(cropto,add=T,border="orange",lwd=1.5)}
  } else {
    return(pa.change)
  }
  
}



# Cells with ok vegetation visible
getMask = function(bg=ps, region=slu, vegmap=subveg, ids=c(4,5)) {
  ### OBSOBS! Added round() due to decimal erros in vegmap!! Needs to be fixed
  vpoints = xyFromCell(vegmap, which(round(getValues(vegmap))%in%ids),spatial=T)
  vc = unique(cellFromXY(bg,vpoints))
  gmask = bg
  gmask[] = 1
  gmask[vc] = NA
  gmask = mask(gmask,slu)
  return(gmask)  
}

# Function to plot suitability map with regions under threshold in grey scale
greyOutTh = function(map, th, legend=F) {
  m1 = c(0,th,NA) # above threshold
  r1 = reclassify(map, m1)
  m2 = c(th,1,NA) # bellow threshold
  r2 = reclassify(map,m2)
  Greys = two.colors(100,"grey90","black","grey45")
  Colors = tim.colors(100)
  plot(r2, col=Greys, zlim=c(0,1),legend=F,axes=F)
  plot(r1, col=Colors,zlim=c(0,1),add=T,legend=F,axes=F)
  if(legend==T) {
    xpos = extent(map)[1]+2000
    ypos1 = extent(map)[3]+5000
    ypos2 = ypos1+2e4
    points(rep(xpos,100),seq(ypos1,ypos2,length=100),pch=15, col="black",cex=2.1)
    points(rep(xpos,100),seq(ypos1,ypos2,length=100),pch=15,
           col=c(Greys[0:round(th*100)],Colors[(round(th*100)+1):100]),cex=2)
    text(rep(xpos,3),seq(ypos1,ypos2,length=3), as.character(seq(0,1,0.5)),pos=4,cex=0.8)
  }
}

# Add legend to existing plot (ok for landscape units, but not pretty!)
trail.legend = function(x0=-225e3,xt=-22e4,y0=2e4,y1=4e4) {
  points(rep(x0,4),seq(y0,y1,length=4),pch=22,col="black",bg=c("red","black","white","blue"))
  text(rep(xt,4),seq(y0,y1,length=4),c("Trail","Same pres","Same abs","Lead"),adj = c(0,0.5))  
}

# Plot current veg, present suitability and 6 future scenarios leading/trailing edges
plotTrailLead = function(ps=ps, fs=fs, mx.th=mx.th, vegmap=subveg, slu=slu, gmask=gmask) {
  par(mfrow=c(4,2),mar=c(0.1,0.1,1.5,0))
  plot(vegmap,axes=F,legend=F)
  plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.9))
  #plot(slu,add=T)
  title("Obs. vegetation")
  greyOutTh(ps, mx.th, legend=T)
  #plot(slu,add=T)
  plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.5))
  title("Present suitability",cex=0.8)
  for(i in 1:nlayers(fs)){
    getOrPlotChange(ps,fs[[i]],mx.th,Plot=T)
    plot(gmask,legend=F,add=T,col=rgb(1,1,1,0.5))
    #plot(slu,add=T)
    if(i == 1) {
      x0 = xmin(fs[[i]])+(xmax(fs[[i]])-xmin(fs[[i]]))/100
      xt = x0 + (xmax(fs[[i]])-xmin(fs[[i]]))/50
      y0 = ymin(fs[[i]])+(ymax(fs[[i]])-ymin(fs[[i]]))/10
      y1 = ymin(fs[[i]])+(ymax(fs[[i]])-ymin(fs[[i]]))/2
      trail.legend(x0=x0,xt=xt,y0=y0,y1=y1)
    }
    scen=strsplit(names(fs[[i]]),"_")[[1]][1:2]
    title(paste("Change",scen[1],scen[2],sep="_"),cex=0.8)  
  }
}

# Function to plot current obs. veg, present suitability and 6 future scenarios
plotPresFut = function(ps=ps, fs=fs, mx.th=mx.th, vegmap=subveg, gmask=gmask) {
  par(mfrow=c(4,2),mar=c(0.1,0.1,1.5,0))
  plot(vegmap,axes=F,legend=F)
  title("Obs. vegetation")
  plot(gmask,add=T,col=rgb(1,1,1,0.5),legend=F)
  greyOutTh(ps, mx.th, legend=T)
  title("Present suitability",cex=0.8)
  for(i in 1:nlayers(fs)){
    greyOutTh(fs[[i]],mx.th)
    scen=strsplit(names(fs[[i]]),"_")[[1]][1:2]
    title(paste(scen[1],scen[2],sep="_"),cex=0.8)  
  }
}

# Function to choose color code based on the mean suitability across a region
# in the future as a fraction of present suitability 
# (0-25%: red, 25-75%: orange, 75-125%: grey, >125%: green)
pickCol = function(psuit,fsuit) {
  if(fsuit < 0.25*psuit) {
    return("red") # if mean future suitability less than 25% of present suitability
  } else if (fsuit < 0.75*psuit) {
    return("orange") # if mean future suitability less than 75% than present suitability
  } else if (fsuit <= 1.25*psuit) {
    return("grey") # if mean future suit. between 75-125% of present suit
  } else {
    return("green") # if mean future suit. larger than present
  }   
}


# Function to plot out summary suitabilities across a region
# The region can be a lanscape unit, a subset of it (e.g. excluding regions where
# the vegetation type is absent), or any other. 
plotSummarySq_old = function(p=ps, f=fs,main) {
  psuit = mean(getValues(p),na.rm=T)
  fsuits = c()
  for(i in 1:nlayers(f)) {
    res = mean(getValues(f[[i]]),na.rm=T)
    fsuits=c(fsuits,res)
  }
  
  #Define scenarios to lump:
  wawet = "PCM_B1"
  howet = c("CCSM4_rcp85","CNRM_rcp85")
  wadry = "GFDL_B1"
  hodry = c("GFDL_A2","MIROC_rcp85")
  Gr = list(wawet,howet,wadry,hodry)
  
  AllScenarios = gsub("_Suitability_Bay","",names(fs))
  
  Means = sapply(Gr,function(x) {mean(fsuits[which(AllScenarios%in%x)])})
  
  # Plot colored summaries of change
  
  myCols = sapply(Means,function(x) {pickCol(psuit=psuit,fsuit=x)})
  plot(x=c(1,2,1,2),y=c(2,2,1,1), pch=22,cex=6,xlim=c(0,3),ylim=c(0,3),axes=F,
       bg=myCols,main=main,cex.main=0.8)
  text(x=c(1,2),y=c(2.5,2.5),c("Warm","Hot"),cex=0.6,pos=3)
  text(x=c(0.3,0.3),y=c(1,2),c("Dry","Wet"),cex=0.6,srt=90)  
}

plotSummarySq = function(p=ps, f=fs,main, giveCols=F) {
  psuit = mean(getValues(p),na.rm=T)
  fsuits = c()
  for(i in 1:nlayers(f)) {
    res = mean(getValues(f[[i]]),na.rm=T)
    fsuits=c(fsuits,res)
  }
  
  #Define scenarios to lump:
  wawet = "PCM_B1"
  howet = c("CCSM4_rcp85","CNRM_rcp85")
  wadry = "GFDL_B1"
  hodry = c("GFDL_A2","MIROC_rcp85")
  Gr = list(wawet,howet,wadry,hodry)
  
  AllScenarios = gsub("_Suitability_Bay","",names(fs))
  
  Means = sapply(Gr,function(x) {mean(fsuits[which(AllScenarios%in%x)])})
  
  # Plot colored summaries of change
  
  myCols = sapply(Means,function(x) {pickCol(psuit=psuit,fsuit=x)})
  plot(-1,xlim=c(0.5,3.5),ylim=c(0.5,3.5),axes=F,xlab="",ylab="")
  rect(xleft=c(1,2,1,2),ybottom=c(2,2,1,1),xright=c(2,3,2,3),ytop=c(3,3,2,2), col=myCols)
  #text(x=c(1,2),y=c(2.5,2.5),c("Warm","Hot"),cex=0.6,pos=3)
  #text(x=c(0.3,0.3),y=c(1,2),c("Dry","Wet"),cex=0.6,srt=90)  
  if(giveCols) {
    return(myCols)
  }
}


# Special function for getting the change colors for areas covered by veg type
# in addition to leading edge areas under each scenario
getSubCols = function(ps,fs,gmask,mx.th)
{
  ChangeMaps = list()
  for (i in 1:nlayers(fs)) {
    ChangeMaps[[i]] = getOrPlotChange(ps,fs[[i]],th=mx.th,Plot=F)
  }
  
  vegcells = which(is.na(getValues(gmask)))
  # create list of masking maps that only include current veg type and leading
  # edges under each future scenario
  Masks = list()
  for (i in 1:length(ChangeMaps)) {
    Change = ChangeMaps[[i]]
    leadcells = which(getValues(Change)==2)
    ocells = which(!(1:ncell(ps)%in%unique(c(vegcells,leadcells))))
    Masks[[i]] = Change
    Masks[[i]][] = 1
    Masks[[i]][ocells] = NA
  }
  
  # make a table with the mean suitabilities
  psuit = sapply(Masks, function(x) {
    subp = mask(ps, x)
    res = mean(getValues(subp),na.rm=T)
    return(res)
  })
  fsuit = c()
  for(i in 1:length(Masks)) {
    subf = mask(fs[[i]], Masks[[i]])
    res = mean(getValues(subf),na.rm=T)
    fsuit = c(fsuit,res)
  }
  suitabilities = data.frame(pres=psuit, fut=fsuit, scen=names(fs), 
                             agscen = c("hodry","howet","hodry","howet","wadry","wawet"))
  
  
  agsuit = aggregate(.~agscen,data=suitabilities, FUN=mean)[,c("pres","fut","agscen")]
  agsuit$color = apply(agsuit[,c("pres","fut")], 1,function(x) {pickCol(x[1],x[2])})
  
  # order as "warmwet", "hotwet", "warmdry", "hotdry" for plotting
  myCols = agsuit$color[match(c("wawet","howet","wadry","hodry"),agsuit$agscen)]
  
}

