###################################################
##      Main Control Script to process obs       ##
###################################################

## Post-process Ameriflux data?
postprocAmf <- TRUE

        # If TRUE, specify Ameriflux data file:
        dataAmfIn <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/AMF/'
	# Specify output data file (if NULL, will overwrite in file):
	dataAmfOut <- NULL


## Post-process USGS streamflow data?
postprocUsgs <- TRUE

        # If TRUE, specify Ameriflux data file:
        dataUsgs <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/'
        # Specify output data file (if NULL, will overwrite in file):
        dataUsgsOut <- NULL

## Use parallel processing? Specify number of cores:
ncores <- 15


###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)

# Multi-core
parallelFlag <- FALSE
if (ncores>1) {
        library(doParallel)
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
        parallelFlag <- TRUE
}


###----------- AMERIFLUX ------------###
if (postprocAmf) {
	load(dataAmfIn)
	if (is.null(dataAmfOut)) dataAmfOut <- dataAmfIn 
	# Derived vars
	obsAmfData$RgNet <- with(obsAmfData, Rg - RgOut)
	obsAmfData$RglNet <- with(obsAmfData, Rgl - RglOut)
	obsAmfData$EnResid <- with(obsAmfData, 
     		(Rg - RgOut) + (Rgl - RglOut) - LE - H - FG)
	obsAmfData$H2O_mmps <- with(obsAmfData, ifelse(TA>0, LE/2510/1000, LE/2844/1000))

	# Date
	obsAmfData$UTC_date <- as.Date(trunc(as.POSIXct(format(obsAmfData$POSIXct, tz="UTC"), tz="UTC"), "days"))

	# Daily means
	obsAmfData.d <- plyr::ddply(obsAmfData, plyr::.(site_id, UTC_date), 
                         plyr::summarise, 
                            TA_mean=mean(TA, na.rm=TRUE),
                            WS_mean=mean(WS, na.rm=TRUE),
                            NEE_mean=mean(NEE, na.rm=TRUE),
                            FC_mean=mean(FC, na.rm=TRUE),
                            H_mean=mean(H, na.rm=TRUE),
                            SH_mean=mean(SH, na.rm=TRUE),
                            LE_mean=mean(LE, na.rm=TRUE),
                            RH_mean=mean(RH, na.rm=TRUE),
                            PRESS_mean=mean(PRESS, na.rm=TRUE),
                            CO2_mean=mean(CO2, na.rm=TRUE),
                            VPD_mean=mean(VPD, na.rm=TRUE),
                            SWC1_mean=mean(SWC1, na.rm=TRUE),
                            SWC2_mean=mean(SWC2, na.rm=TRUE),
                            Rn_mean=mean(Rn, na.rm=TRUE),
                            Rg_mean=mean(Rg, na.rm=TRUE),
                            RgOut_mean=mean(RgOut, na.rm=TRUE),
                            Rgl_mean=mean(Rgl, na.rm=TRUE),
                            RglOut_mean=mean(RglOut, na.rm=TRUE),
                            H2O_mean=mean(H2O, na.rm=TRUE),
                            RE_mean=mean(RE, na.rm=TRUE),
                            GPP_mean=mean(GPP, na.rm=TRUE),
                            ZL_mean=mean(ZL, na.rm=TRUE),
                            H2Ommps_mean=mean(H2O_mmps, na.rm=TRUE),
                         .parallel=parallelFlag)

	# Add a POSIXct date
	obsAmfData.d$POSIXct <- as.POSIXct(paste0(obsAmfData.d$UTC_date," 00:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")

	# Unit conversions
	obsAmfData.d$wy <- CalcWaterYear(obsAmfData.d$POSIXct)
	obsAmfData.d$H2O_mmpd <- with(obsAmfData.d, H2Ommps_mean*86400)

	# Save outputs
	save(obsAmfData, obsAmfMeta, obsAmfData.d, file=dataAmfOut)
}


###----------- USGS ------------###

if (postprocUsgs) {
	load(dataUsgsIn)
	if (is.null(dataUsgsOut)) dataUsgsOut <- dataUsgsIn

	# Aggregate the observed data to a daily timestep.
	obsStrData$UTC_date <- as.Date(trunc(as.POSIXct(format(obsStrData$POSIXct, tz="UTC"), tz="UTC"), "days"))
	obsStrData.d <- plyr::ddply(obsStrData, plyr::.(site_no, UTC_date),
                         plyr::summarise, mean_qcms=mean(q_cms, na.rm=TRUE),
                         .parallel=parallelFlag)
	# Unit conversion: m^3/s -> m^3/dy -> ft^3/dy -> ac-ft/dy
	obsStrData.d$qvol_acft <- obsStrData.d$mean_qcms * 86400 / (0.3048^3) / 43560

	# Add a POSIXct column for ease of calculations and plotting.
	obsStrData.d$POSIXct <- as.POSIXct(paste0(obsStrData.d$UTC_date," 00:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")

	# Calculate a cumulative volume
	obsStrData.d$wy <- CalcWaterYear(obsStrData.d$POSIXct)
	obsStrData.d <- obsStrData.d[order(obsStrData.d$site_no, obsStrData.d$UTC_date),]
	wyList <- unique(obsStrData.d$wy)
	gageList <- unique(obsStrData.d$site_no)
	obsStrData.d$cumqvol_acft <- 0
	obsStrData.d$cumqvol_mm <- 0
	for (i in 1:length(gageList)) {
  		print(paste0("index ",i, " gage ", gageList[[i]]))
  		tmpgage <- subset(obsStrData.d, obsStrData.d$site_no==gageList[i])
    		obsStrData.d$cumqvol_acft[obsStrData.d$site_no==gageList[i]] <- CumsumNa(tmpgage$qvol_acft)
    		obsStrData.d$cumqvol_mm[obsStrData.d$site_no==gageList[i]] <- CumsumNa(tmpgage$qvol_acft) /
      						obsStrMeta$area_sqmi[obsStrMeta$site_no==gageList[i]] /
      						(5280^2) * 43560 * 25.4 * 12
    		rm(tmpgage)
  	}

	# Save outputs
	save(obsStrData, obsStrMeta, obsStrData.d, file=dataUsgsOut)
}

quit()
