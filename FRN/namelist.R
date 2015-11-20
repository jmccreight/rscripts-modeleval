###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doParallel package installed
ncores <- 15

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/DOMAIN/AD_KS_151117/Fulldom_hires_netcdf_file_151117.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/DOMAIN/AD_KS_151117/geo_em.d02_151117.nc'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 10

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/DOMAIN/AD_KS_151117/frn_4basns_MASKS.Rdata'

## Specify manual link2gage, if it exists
link2gage.man <- data.frame(site_no=c("06721000","SVCPLACO","BIGLASCO","09058000"), link=c(225621,2890908,14962,1233053))

################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- NULL

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
STRfile <- c("/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/OBS/obsStrData_USGS.Rdata",
		"/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/OBS/obsStrData_CODWR.Rdata")


################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
        modPathList <- c('/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/test_basn4_gochis_reach')

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
        modTagList <- c('BASN4 Reach')

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/ANALYSIS/151118_frn_modelreads.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        modAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- TRUE
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- FALSE

		# Channel routing
		readChrtout <- TRUE

		# Snotel sites
		readSnoLdasout <- FALSE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- FALSE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- FALSE
	varsLdasoutNFIE <- TRUE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- NULL
	readModEnd <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/scratch/zhangyx/WRF-Hydro/CONUS1km.ForcingData/2013') 

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
	forcTagList <- c('NLDAS2-Downscaled')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_spinup2013_forcings.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        forcAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means
		readBasinLdasin <- FALSE

		# SNOTEL sites
		readSnoLdasin <- TRUE

		# Ameriflux sites
		readAmfLdasin <- TRUE

		# MET sites
		readMetLdasin <- FALSE

        # Specify start and end dates if you do NOT want to read all files
        readForcStart <- as.POSIXct("2013-01-01", format="%Y-%m-%d", tz="UTC")
        readForcEnd <- NULL


############# Model Performance Stats #############

## Calculate stats?
calcStats <- FALSE

	## Calculate streamflow performance stats?
	strProc <- TRUE

	## Calculate SNOTEL performance stats?
	snoProc <- FALSE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- FALSE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/ANALYSIS/151118_frn_modelreads.Rdata'
	forcReadFileIn <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_forcingreads.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/ANALYSIS/151118_frn_stats.Rdata'

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- NULL

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2013-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2013-10-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/ANALYSIS/spinup2013_PLOTS'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- TRUE

## Create HTML files?
writeHtml <- TRUE

## If TRUE, specify output directory
writePlotDir <- '/glade/p/ral/RHAP/adugger/FRNG_MASTER_15Min_LONG_16Cores/ANALYSIS/spinup2013_PLOTS'

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- FALSE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- NULL

		# Specify end date
		accflowEndDate <- NULL

	## Generate hydrographs?
	hydroPlot <- TRUE

        	# Specify which run tags to plot
        	hydroTags <- NULL
        
        	# Specify start date
        	hydroStartDate <- NULL
        
        	# Specify end date
        	hydroEndDate <- NULL

	## Generate accumulated precip plots?
	accprecipPlot <- FALSE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- NULL
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- FALSE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- NULL
        
        	# Specify end date
        	flowsweEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- FALSE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- NULL

        	# Specify end date
        	sweEndDate <- NULL

        ## Generate MET station plots?
        metPlot <- FALSE

                # Specify which run tags to plot
                metTags <- NULL

                # Specify start date
                metStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

                # Specify end date
                metEndDate <- NULL


	########### MAPS #############

	## Generate STRFLOW bias maps?
	strBiasMap <- FALSE

        	# Specify which run tags to plot
        	strBiasTags <- NULL

        	# Specify which run seasons to plot
        	strBiasSeas <- NULL

	## Generate STRFLOW correlation maps?
	strCorrMap <- FALSE

        	# Specify which run tags to plot
        	strCorrTags <- NULL

        	# Specify which run seasons to plot
        	strCorrSeas <- NULL

	## Generate SNOTEL SWE error maps?
	snosweErrMap <- FALSE

        	# Specify which run tags to plot
        	snosweErrTags <- NULL

        	# Specify which run seasons to plot
        	snosweErrSeas <- NULL

	## Generate SNOTEL Precip error maps?
	snoprecipErrMap <- FALSE

        	# Specify which run tags to plot
        	snoprecipErrTags <- NULL

        	# Specify which run seasons to plot
        	snoprecipErrSeas <- NULL



###########################################################################################
## RUN (do not change anything below this line)

source("run_CONFIG.R")

