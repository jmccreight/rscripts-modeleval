###################################################################################################
#' # Setup
#' 
#' Load the rwrfhydro package. 
## ------------------------------------------------------------------------
library("rwrfhydro")

# Create geogrid from predefined masks?
man2geo <- TRUE
# Where geogrid basin mask files live
# Assumes these are already at the geogrid resolution
mskgeoPath <- '../DOMAIN/MASKS'
# ID value for basin mask files
basid <- 1
# South->north orientation?
# Geogrids are south->north, hydrogrids are north->south
south_north <- FALSE

# Where the high-res routing grid lives
hydPath <- '../DOMAIN/Fulldom_hires_netcdf_file.nc'
# Create geogrid masks from high-res masks?
hyd2geo <- TRUE
# Aggregation factor between hydrogrid and geogrid
# Only necessary if hyd2geo is TRUE
aggfact <- 10

# Where to save the R workspace
rimgPath <- 'urg_masks_NEW.Rdata'

#' Create a lookup list to match low-res basin mask ID to gage ID.
## ------------------------------------------------------------------------
file2gageList <- list("alamosa_1k"="ALATERCO", 
               "conejos_1k"="CONMOGCO", 
               "s_frk_1k"="RIOSFKCO", 
               "rio_wagw_1k"="RIOWAGCO", 
               "saguache_1k"="SAGSAGCO", 
               "trinchera_1k"="TRITURCO",
               "pinos_1k"="08248000",
               "plat_inflo_1k"="PLATERO_INFLOW",
               "rio_deln_1k"="RIODELCO",
               "sn_anton_1k"="SANORTCO")

#' Create a lookup list to match high-res basin mask ID to gage ID.
## ------------------------------------------------------------------------
gage2basinList <- list("ALATERCO"=c(9,10,12), 
               "CONMOGCO"=c(2,35,36,37,38,39,40,3,33,34), 
               "RIOSFKCO"=c(17), 
               "RIOWAGCO"=c(15,18), 
               "SAGSAGCO"=c(21), 
               "TRITURCO"=c(27),
               "08248000"=c(4),
               "CONPLACO"=c(3,33,34),
               "RIODELCO"=c(14,15,17,18),
               "SANORTCO"=c(8))
# Reverse this list
basin2gageList <- list()
for (i in unique(unlist(gage2basinList))) {
  idlist <- c()
  for (j in names(gage2basinList)) {
    if (i %in% gage2basinList[[j]]) idlist <- c(idlist, j)
    }
  basin2gageList[[paste0(i)]] <- idlist
  }

#' Create a lookup list to match high-res frxst pt ID to gage ID.
## ------------------------------------------------------------------------
stid2gageList <- list("21"="ALATERCO",
                     "4"="CONMOGCO",
                     "35"="RIODELCO",
                     "34"="RIOSFKCO",
                     "39"="RIOWAGCO",
                     "46"="SAGSAGCO",
                     "2"="SANORTCO",
                     "22"="TRITURCO",
                     "0"="08248000",
                     "20"="CONPLACO")

###################################################################################################
#' # Run
#' 

#' ### Create manual GEOGRID basin mask objects
if (man2geo) {
  
  mskgeo.fileList <- list.files(mskgeoPath, pattern=glob2rx("*.nc"))

  # Initialize
  mskgeoman.List <- list()
  mskgeoman.nameList <- list()
  mskgeoman.areaList <- list()
  mskgeoman.minInds <- data.frame(x=integer(0), y=integer(0))
  mskgeoman.maxInds <- data.frame(x=integer(0), y=integer(0))
  mskgeoman.countInds <- data.frame(x=integer(0), y=integer(0))

  # Loop through basin masks
  for (i in 1:length(mskgeo.fileList)) {
    print(paste0("Basin: ", mskgeo.fileList[i]))
    # Setup mask
    ncid <- ncdf4::nc_open(paste0(mskgeoPath, "/", mskgeo.fileList[i]))
    mskvar <- ncdf4::ncvar_get(ncid, "basn_msk")
    ncdf4::nc_close(ncid)
    # Subset to basinID
    mskvar[which(mskvar != basid)] <- 0.0
    mskvar[which(mskvar == basid)] <- 1.0
    # Reverse y-direction for N->S hydro grids to S->N
    if (!south_north) {
      mskvar <- mskvar[,order(ncol(mskvar):1)]
    }
    # Calculate cell indices
    mskgeoman.minInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=min)
    mskgeoman.maxInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=max)
    mskgeoman.countInds[i,1:2] <- mskgeoman.maxInds[i,1:2] - mskgeoman.minInds[i,1:2] + 1
    # Calculate basin area as a cell count
    basarea <- sum(mskvar)
    mskgeoman.List[[i]] <- mskvar[mskgeoman.minInds[i,1]:mskgeoman.maxInds[i,1], mskgeoman.minInds[i,2]:mskgeoman.maxInds[i,2]]
    mskgeoman.nameList[[i]] <- file2gageList[[unlist(strsplit(mskgeo.fileList[i], "[.]"))[1]]]
    mskgeoman.areaList[[i]] <- basarea
  }
  
  # Add names/IDs
  names(mskgeoman.List) <- unlist(mskgeoman.nameList)
  names(mskgeoman.areaList) <- unlist(mskgeoman.nameList)
  mskgeoman.minInds$id <- unlist(mskgeoman.nameList)
  mskgeoman.maxInds$id <- unlist(mskgeoman.nameList)
  mskgeoman.countInds$id <- unlist(mskgeoman.nameList)
  
}

#' ### Create HYDROGRID basin mask objects

# Initialize
mskhyd.List <- list()
mskhyd.nameList <- list()
mskhyd.areaList <- list()
mskhyd.minInds <- data.frame(x=integer(0), y=integer(0))
mskhyd.maxInds <- data.frame(x=integer(0), y=integer(0))
mskhyd.countInds <- data.frame(x=integer(0), y=integer(0))

if (hyd2geo) {
  # Initialize
  mskgeo.List <- list()
  mskgeo.nameList <- list()
  mskgeo.areaList <- list()
  mskgeo.minInds <- data.frame(x=integer(0), y=integer(0))
  mskgeo.maxInds <- data.frame(x=integer(0), y=integer(0))
  mskgeo.countInds <- data.frame(x=integer(0), y=integer(0))
}

ncid <- ncdf4::nc_open(hydPath)
mskvarAll <- ncdf4::ncvar_get(ncid, "basn_msk")
ncdf4::nc_close(ncid)

# Loop through basin masks
for (i in 1:length(gage2basinList)) {
  print(paste0("Basin: ", names(gage2basinList)[i]))
  # Subset to basinID
  mskvar <- mskvarAll
  mskvar[which(!(mskvar %in% gage2basinList[[i]]))] <- 0.0
  mskvar[which(mskvar %in% gage2basinList[[i]])] <- 1.0
  # Reverse y-direction for N->S hydro grids to S->N
  mskvar <- mskvar[,order(ncol(mskvar):1)]
  # Calculate cell indices
  mskhyd.minInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=min)
  mskhyd.maxInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=max)
  mskhyd.countInds[i,1:2] <- mskhyd.maxInds[i,1:2] - mskhyd.minInds[i,1:2] + 1
  # Calculate basin area as a cell count
  basarea <- sum(mskvar)
  mskhyd.List[[i]] <- mskvar[mskhyd.minInds[i,1]:mskhyd.maxInds[i,1], mskhyd.minInds[i,2]:mskhyd.maxInds[i,2]]
  mskhyd.nameList[[i]] <- names(gage2basinList)[i]
  mskhyd.areaList[[i]] <- basarea
  if (hyd2geo) {
    # Resample the high-res grid to the low-res LSM
    if (aggfact > 1) {
      mskvar <- raster::as.matrix(raster::aggregate(raster::raster(mskvar), fact=aggfact, fun=mean))
    }
    # Calculate cell indices
    mskgeo.minInds[i,1:2] <- apply(which(mskvar>0, arr.ind=TRUE), MARGIN=2, FUN=min)
    mskgeo.maxInds[i,1:2] <- apply(which(mskvar>0, arr.ind=TRUE), MARGIN=2, FUN=max)
    mskgeo.countInds[i,1:2] <- mskgeo.maxInds[i,1:2] - mskgeo.minInds[i,1:2] + 1
    # Calculate basin area as a cell count
    basarea <- sum(mskvar)
    mskgeo.List[[i]] <- mskvar[mskgeo.minInds[i,1]:mskgeo.maxInds[i,1], mskgeo.minInds[i,2]:mskgeo.maxInds[i,2]]
    mskgeo.nameList[[i]] <- names(gage2basinList)[i]
    mskgeo.areaList[[i]] <- basarea
  }
} 

# Add names/IDs
# HYDRO
names(mskhyd.List) <- unlist(mskhyd.nameList)
names(mskhyd.areaList) <- unlist(mskhyd.nameList)
mskhyd.minInds$id <- unlist(mskhyd.nameList)
mskhyd.maxInds$id <- unlist(mskhyd.nameList)
mskhyd.countInds$id <- unlist(mskhyd.nameList)
# GEO
if (hyd2geo) {
  names(mskgeo.List) <- unlist(mskgeo.nameList)
  names(mskgeo.areaList) <- unlist(mskgeo.nameList)
  mskgeo.minInds$id <- unlist(mskgeo.nameList)
  mskgeo.maxInds$id <- unlist(mskgeo.nameList)
  mskgeo.countInds$id <- unlist(mskgeo.nameList)
}

# Save
if (exists("mskgeoman.List")) {
  save(list=c("mskgeo.nameList", "mskgeo.areaList", "mskgeo.List", 
              "mskgeo.minInds", "mskgeo.maxInds", "mskgeo.countInds",
              "mskhyd.nameList", "mskhyd.areaList", "mskhyd.List", 
              "mskhyd.minInds", "mskhyd.maxInds", "mskhyd.countInds",
              "mskgeoman.nameList", "mskgeoman.areaList", "mskgeoman.List", 
              "mskgeoman.minInds", "mskgeoman.maxInds", "mskgeoman.countInds",
              "gage2basinList", "basin2gageList", "stid2gageList"), 
       file=rimgPath)
} else {
  save(list=c("mskgeo.nameList", "mskgeo.areaList", "mskgeo.List", 
            "mskgeo.minInds", "mskgeo.maxInds", "mskgeo.countInds",
            "mskhyd.nameList", "mskhyd.areaList", "mskhyd.List", 
            "mskhyd.minInds", "mskhyd.maxInds", "mskhyd.countInds",
            "gage2basinList", "basin2gageList", "stid2gageList"), 
    file=rimgPath)
}



