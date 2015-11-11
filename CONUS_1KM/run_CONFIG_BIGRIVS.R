###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doMC package installed
ncores <- 15

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/Fulldom_hires_netcdf_file_1km.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/geo_em.d01.nc.conus_1km'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 1

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/CONUS_1km/DOMAIN/bigrivs_MASKS.Rdata'


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- NULL 

## Path to SNOTEl data .Rdata file
#SNOfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG_NEW.Rdata" 
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_1km/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
#STRfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG_NEW.Rdata"
STRfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/obsStrData.Rdata"

## Do these need updated? TRUE/FALSE ...to be worked on later


################ Model Output Reads ###############

## Read model output?
readMod <- TRUE

## If TRUE, specify the following to read in model output:

	# Specify the model run output directory or directories
	modPathList <- c('/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/NHDPLUS_Run/output_new/cold_start_no_terr_rtg')

	# Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('SPINUP2013_1')

	# Specify the output .Rdata file to create
	modReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_spinup2013_eval_1.Rdata'
	# Append to existing file? FALSE = create new file (or overwrite existing!)
	modAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- TRUE  
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- FALSE

		# Snotel sites
		readSnoLdasout <- FALSE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- FALSE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- TRUE
	varsLdasoutNFIE <- FALSE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- NULL
	readModEnd <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	#forcPathList <- c('/glade/scratch/zhangyx/WRF-Hydro/RioGrande/NLDAS2.data')
	forcPathList <- c('/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/forcing', 
			'/glade/scratch/zhangyx/WRF-Hydro/RioGrande/NLDAS2.data')

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
        #forcTagList <- c('nldas2dwnsc')
	forcTagList <- c('nldas2', 'nldas2dwnsc')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_forcing_raw2.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        forcAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means
		readBasinLdasin <- TRUE

		# SNOTEL sites
		readSnoLdasin <- TRUE

		# Ameriflux sites
		readAmfLdasin <- FALSE

		# MET sites
		readMetLdasin <- TRUE

        # Specify start and end dates if you do NOT want to read all files
        readForcStart <- as.POSIXct("2013-01-01", format="%Y-%m-%d", tz="UTC")
        readForcEnd <- NULL


############# Model Performance Stats #############

## Calculate stats?
calcStats <- FALSE

	## Calculate streamflow performance stats?
	strProc <- TRUE

	## Calculate SNOTEL performance stats?
	snoProc <- TRUE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- TRUE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- NULL

        # Specify the stats output .Rdata file to create
        statsFileOut <- NULL

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- NULL



################### Plotting ######################

## Create plots and/or maps?
createPlots <- FALSE

## If TRUE, specify output directory
writePlotDir <- NULL

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- TRUE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- as.POSIXct("2015-04-01", format="%Y-%m-%d", tz="UTC")

		# Specify end date
		accflowEndDate <- NULL

	## Generate hydrographs?
	hydroPlot <- TRUE

        	# Specify which run tags to plot
        	hydroTags <- NULL
        
        	# Specify start date
        	hydroStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	hydroEndDate <- NULL

	## Generate accumulated precip plots?
	accprecipPlot <- FALSE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- FALSE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	flowsweEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- TRUE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

        	# Specify end date
        	sweEndDate <- NULL

	########### MAPS #############

	## Generate STRFLOW bias maps?
	strBiasMap <- TRUE

        	# Specify which run tags to plot
        	strBiasTags <- NULL

        	# Specify which run seasons to plot
        	strBiasSeas <- NULL

	## Generate STRFLOW correlation maps?
	strCorrMap <- TRUE

        	# Specify which run tags to plot
        	strCorrTags <- NULL

        	# Specify which run seasons to plot
        	strCorrSeas <- NULL

	## Generate SNOTEL SWE error maps?
	snosweErrMap <- TRUE

        	# Specify which run tags to plot
        	snosweErrTags <- NULL

        	# Specify which run seasons to plot
        	snosweErrSeas <- NULL

	## Generate SNOTEL Precip error maps?
	snoprecipErrMap <- TRUE

        	# Specify which run tags to plot
        	snoprecipErrTags <- NULL

        	# Specify which run seasons to plot
        	snoprecipErrSeas <- NULL



###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)

load(maskFile)
source("util_FUNC.R")

# Multi-core
parallelFlag <- FALSE
if (ncores>1) {
	library(doParallel)
	cl <- makeForkCluster(ncores)
	registerDoParallel(cl)
	parallelFlag <- TRUE
}

# Obs
if (!is.null(AMFfile) & amfProc & exists("ptgeo.amf")) {
	if (file.exists(AMFfile)) {
		load(AMFfile)
		obsAmfData <- subset(obsAmfData, obsAmfData$site_id %in% ptgeo.amf$id)
	} else {
		stop(paste("Ameriflux obs file specified but does not exist:", AMFfile))
	}
}
if (!is.null(SNOfile) & (snoProc | swePlot) & exists("ptgeo.sno")) {
        if (file.exists(SNOfile)) {
                load(SNOfile)
		obsSnoData <- subset(obsSnoData, obsSnoData$site_id %in% ptgeo.sno$id)
	} else {
                stop(paste("SNOTEL obs file specified but does not exist:", SNOfile))
        }
}
if (!is.null(METfile) & metProc & exists("ptgeo.met")) {
        if (file.exists(METfile)) {
                load(METfile)
		obsMetData <- subset(obsMetData, obsMetData$site_id %in% ptgeo.met$id)
	} else {
                stop(paste("MET obs file specified but does not exist:", METfile))
        }
}
if (!is.null(STRfile) & (strProc | accflowPlot | hydroPlot | flowswePlot) & exists("stid2gageList")) {
	obsStrData_FINAL <- data.frame()
	obsStrMeta_FINAL <- data.frame()
	for (i in STRfile) {
		if (file.exists(i)) {
			load(i)
			obsStrData <- remapData(obsStrData, obsStrData.map)
			obsStrMeta <- remapData(obsStrMeta, obsStrMeta.map)
			#if (!("site_no" %in% names(obsStrData))) names(obsStrData)[which(names(obsStrData)=="Station")] <- "site_no"
			#if (!("site_no" %in% names(obsStrMeta))) names(obsStrMeta)[which(names(obsStrMeta)=="Station")] <- "site_no"
			obsStrData_TMP <- subset(obsStrData, obsStrData$site_no %in% stid2gageList)
			obsStrData_FINAL <- plyr::rbind.fill(obsStrData_FINAL, obsStrData_TMP)
			obsStrMeta_TMP <- subset(obsStrMeta, obsStrMeta$site_no %in% stid2gageList)
                        obsStrMeta_FINAL <- plyr::rbind.fill(obsStrMeta_FINAL, obsStrMeta_TMP)
		} else {
			stop(paste("Streamflow obs file specified but does not exist:", STRfile))
		}
	}
	obsStrData <- obsStrData_FINAL
	obsStrMeta <- obsStrMeta_FINAL
	rm(obsStrData_FINAL, obsStrMeta_FINAL, obsStrData_TMP, obsStrMeta_TMP)
}

stop()
# Model Reads 
if (readMod | readForc) {
	source("read_MODELOUT.R")
}

# Stats Calculations
if (calcStats & (strProc | snoProc | amfProc | metProc)) {
	message("Calculating stats")
	if (is.null(modReadFileOut)) {
		if (file.exists(modReadFileIn)) {
			load(modReadFileIn)
		}
	} else {
		if (is.null(modReadFileIn)) {
			if (file.exists(modReadFileOut)) {
				load(modReadFileOut)
			}
		} else {
			if (file.exists(modReadFileIn)) {
				load(modReadFileIn)
			}
		}
	}
	source("calc_PERFSTATS.R")
}

# Plots
if (createPlots & (accflowPlot | hydroPlot | accprecipPlot | swePlot | 
	strBiasMap | strCorrMap | 
	snosweErrMap | snoprecipErrMap)) {
        message("Generating plots")
	load(statsFileOut)
        if (is.null(modReadFileOut)) {
                if (file.exists(modReadFileIn)) {
                        load(modReadFileIn)
                }
        } else {
                if (is.null(modReadFileIn)) {
                        if (file.exists(modReadFileOut)) {
                                load(modReadFileOut)
                        }
                } else {
                        if (file.exists(modReadFileIn)) {
                                load(modReadFileIn)
                        }
                }
        }
        source("calc_PLOTS.R")
}

# EXIT
stopCluster(cl)
quit("no")

