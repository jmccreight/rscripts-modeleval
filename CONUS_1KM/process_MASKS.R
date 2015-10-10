###################################################
##               Generate masks                  ##
###################################################

##################### Setup #######################

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/nral0008/zhangyx/CONUS1km_LSMOnly_daily_snowmods/DOMAIN/geo_em.d01.nc.conus_1km'

## Select which masks/points to create:

# Basin masks:
createBasMask <- TRUE

    # Specify the mask file
    maskFile <- '/glade/p/ral/RHAP/adugger/CONUS_1km/DOMAIN/gages2_1km.nc'
    # If relevant, specify variable name
    maskVar <- "Band1"
    # Specify the aggregation factor between the input mask grid and geogrid
    aggfact <- 1
    # Reverse the y direction from N->S to S->N?
    ns2sn <- FALSE

# Point masks:
    # Ameriflux points
    createAmfMask <- TRUE
    # SNOTEL points
    createSnoMask <- TRUE
    # MET station points
    createMetMask <- FALSE
    # MET station sites (must include columns: id, lat, lon)
    metSites <- NULL

## If available, specify the list that matches basin ID to gage ID
#  EXAMPLE: gage2basinList <- list("ALATERCO"=c(9,10,12), "CONMOGCO"=c(2,35,36,37,38,39,40,3,33,34), "SANORTCO"=c(8))
gage2basinList <- NULL

## If available, specify a lookup list to match high-res frxst pt ID to gage ID.
#  EXAMPLE: stid2gageList <- list("21"="ALATERCO", "4"="CONMOGCO", "2"="SANORTCO")
stid2gageList <- NULL

# Specify the .Rdata file to create
maskFileOut <- "conus1km_masks.Rdata"

###################################################################################################
## Run 

library(rwrfhydro)
saveList <- c("gage2basinList", "basin2gageList")

if (is.null(maskVar)) maskVar <- "basn_msk"
if (is.null(aggfact)) aggfact <- 1
if (!is.null(stid2gageList)) saveList <- c(saveList, "stid2gageList")

## Create masks for point obs
if (createAmfMask) {
    ptgeo.amf <- GetGeogridIndex(amfMeta, geoFile, id="id_txt")
    ptgeo.amf <- subset(ptgeo.amf, !is.na(ptgeo.amf$sn))
    saveList <- c(saveList, "ptgeo.amf")
}
if (createSnoMask) {
    ptgeo.sno <- GetGeogridIndex(snotelMeta, geoFile, id="site_id")
    ptgeo.sno <- subset(ptgeo.sno, !is.na(ptgeo.sno$sn))
        saveList <- c(saveList, "ptgeo.sno")
}
if (createMetMask) {
    ptgeo.met <- GetGeogridIndex(metSites, geoFile, id="id")
    ptgeo.met <- subset(ptgeo.met, !is.na(ptgeo.met$sn))
        saveList <- c(saveList, "ptgeo.met")
}

## Create HYDROGRID basin mask objects

if (createBasMask) {
    # Initialize - GEO
    mskgeo.List <- list()
    mskgeo.nameList <- list()
    mskgeo.areaList <- list()
    mskgeo.minInds <- data.frame(x=integer(0), y=integer(0))
    mskgeo.maxInds <- data.frame(x=integer(0), y=integer(0))
    mskgeo.countInds <- data.frame(x=integer(0), y=integer(0))
    if (aggfact > 1) {
        # Initialize - HYDRO
        mskhyd.List <- list()
        mskhyd.nameList <- list()
        mskhyd.areaList <- list()
        mskhyd.minInds <- data.frame(x=integer(0), y=integer(0))
        mskhyd.maxInds <- data.frame(x=integer(0), y=integer(0))
        mskhyd.countInds <- data.frame(x=integer(0), y=integer(0))
    }
    # Open basn_msk variable

    ncid <- tryCatch(suppressWarnings(ncdf4::nc_open(maskFile)),  
        error=function(cond) {message(cond); return(c())})
    if (length(ncid)>0) {
        mskvarAll <- ncdf4::ncvar_get(ncid, maskVar)
        ncdf4::nc_close(ncid)
    } else {
        stop("Input mask must be NetCDF format.")
    }
    mskvarAll[which(is.na(mskvarAll))] <- 0.0
    if (ns2sn) {
        # Reverse y-direction for N->S hydro grids to S->N
        mskvarAll <- mskvarAll[,order(ncol(mskvarAll):1)]
    }
    # Get unique ID list
    if (is.null(gage2basinList)) {
        mskidList <- na.exclude(unique(c(mskvarAll)))
        mskidList <- subset(mskidList, mskidList>0)
        mskidList <- mskidList[order(mskidList)]
        gage2basinList <- as.list(mskidList)
        names(gage2basinList) <- mskidList
    }
    # Reverse the list
    basin2gageList <- list()
    for (i in unique(unlist(gage2basinList))) {
        idlist <- c()
        for (j in names(gage2basinList)) {
            if (i %in% gage2basinList[[j]]) idlist <- c(idlist, j)
        }
        basin2gageList[[paste0(i)]] <- idlist
    }
    # Loop through basin masks
    if (length(gage2basinList) > 0) {
        for (i in 1:length(gage2basinList)) {
            print(paste0("Basin: ", names(gage2basinList)[i]))
            # Subset to basinID
            mskvar <- mskvarAll
            mskvar[which(!(mskvar %in% gage2basinList[[i]]))] <- 0.0
            mskvar[which(mskvar %in% gage2basinList[[i]])] <- 1.0
            ##Reverse y-direction for N->S hydro grids to S->N
            ##mskvar <- mskvar[,order(ncol(mskvar):1)]
            # Resample the high-res grid to the low-res LSM
            if (aggfact > 1) {
                # Calculate cell indices
                mskhyd.minInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=min)
                mskhyd.maxInds[i,1:2] <- apply(which(mskvar==1, arr.ind=TRUE), MARGIN=2, FUN=max)
                mskhyd.countInds[i,1:2] <- mskhyd.maxInds[i,1:2] - mskhyd.minInds[i,1:2] + 1
                # Calculate basin area as a cell count
                basarea <- sum(mskvar)
                mskhyd.List[[i]] <- mskvar[mskhyd.minInds[i,1]:mskhyd.maxInds[i,1], mskhyd.minInds[i,2]:mskhyd.maxInds[i,2]]
                mskhyd.nameList[[i]] <- names(gage2basinList)[i]
                mskhyd.areaList[[i]] <- basarea
                # Resample the high-res grid to the low-res LSM
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
        # Add names/IDs
        # HYDRO
        if (aggfact > 1) {
            names(mskhyd.List) <- unlist(mskhyd.nameList)
            names(mskhyd.areaList) <- unlist(mskhyd.nameList)
            mskhyd.minInds$id <- unlist(mskhyd.nameList)
            mskhyd.maxInds$id <- unlist(mskhyd.nameList)
            mskhyd.countInds$id <- unlist(mskhyd.nameList)
            saveList <- c(saveList, "mskhyd.nameList", "mskhyd.areaList", "mskhyd.List",
                            "mskhyd.minInds", "mskhyd.maxInds", "mskhyd.countInds")
        }
        # GEO
        names(mskgeo.List) <- unlist(mskgeo.nameList)
        names(mskgeo.areaList) <- unlist(mskgeo.nameList)
        mskgeo.minInds$id <- unlist(mskgeo.nameList)
        mskgeo.maxInds$id <- unlist(mskgeo.nameList)
        mskgeo.countInds$id <- unlist(mskgeo.nameList)
        # Add to save object list
        saveList <- c(saveList, "mskgeo.nameList", "mskgeo.areaList", "mskgeo.List",
                            "mskgeo.minInds", "mskgeo.maxInds", "mskgeo.countInds")
        } else {
            warning("No basins specified in the high-res domain file.")
        }
}

# Save all relevant objects
save(list=saveList, file=maskFileOut)

quit(save="no")
