# Scripts to load in species specific data for plotting trailing and leading 
# edges.
# If files are missing, relevant models and projections are run
# The script requires that directories have been defined already, as well as species name


# Load occurrence data. Create if it doesn't exist.
occurData = paste(sp.data.dir,"/", mySpecies, ".rdata", sep="")
  if(!file.exists(occurData)) {
    print(paste("No occurrence data. Extracting from Ecoengine species...",mySpecies))
    getOccur(mySpecies=mySpecies, db="Ecoengine", out.dir=sp.data.dir, 
             in.project=orig.project, out.project=ta.project, save=T)
  }
occur = readRDS(occurData)
rm(occurData)

# Load maxent model. Run model if it doesn't exist.
mx.obj = paste(out.dir, mySpecies, mxModelType,'ModelObject', sep="/")
  if(!file.exists(mx.obj)) {
    print(paste("No maxent model found. Computing for...", mxModelType))
    env.files <- list.files(path=pres.clim.data.dir, pattern='CA.img', full.names=FALSE)
    env.files = env.files[-grep("xml",env.files)]
    runMxModel(mySpecies, mxModelType, env.files, pres.clim.data.dir, sp.data.dir, out.dir)
    rm(env.files)
  }
mx = readRDS(mx.obj)
rm(mx.obj)

# Load present projection for Bay Area. Project if it doesn't exist.
proj.dir = paste(out.dir, mySpecies, mxModelType, "Projections", sep="/")
pres.pred = paste(proj.dir, "Present_Bay.img", sep="/")
  if(!file.exists(pres.pred)) {
    print(paste("No present projection. Getting climate data..."))
    pres.clim = getClim(env.dir=pres.clim.data.dir, region="Bay", period="HST",mxModelType=mxModelType)
    print("Projecting...")
    names(pres.clim) = paste(names(pres.clim),"1981_2010_ave_HST_CA",sep="") #match names for mx object
    if(!file.exists(proj.dir)) {dir.create(proj.dir, recursive=F)}
    predict(mx, pres.clim, filename=paste(proj.dir, "Present_Bay.img", sep="/"), overwrite=F)
    rm(pres.clim)
  }
present.bay = raster(pres.pred)
rm(pres.pred)

# Load future projections (all scenarios) for Bay Area. Project if they don't exist.
fut.preds = paste(proj.dir, "/", allScenarios, "_Suitability_Bay.img",sep="")
  if(!file.exists(fut.preds[1])){
    for (i in 1:length(allScenarios)) {
      myScenario = allScenarios[i]
      future.clim = getClim(env.dir=fut.clim.data.dir, period="Future",
                                   mxModelType=mxModelType, scen=myScenario)
      names(future.clim) = paste(names(future.clim),"1981_2010_ave_HST_CA",sep="") #match names for mx object
      print(paste("Projecting...", myScenario))
      if(!file.exists(proj.dir)) {dir.create(proj.dir, recursive=F)}
      predict(mx, future.clim, file=paste(proj.dir, "/", myScenario, "_Suitability_Bay.img", sep=""), overwrite=F)
      rm(future.clim,myScenario)
    }
  }
futures.bay = stack(fut.preds)
rm(fut.preds)

# Set directory to write figures to.
fig.dir =  paste(out.dir, mySpecies, "Figures", sep="/")
if(file.exists(fig.dir) == F) {dir.create(fig.dir, recursive=F)} 
