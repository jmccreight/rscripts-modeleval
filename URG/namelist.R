###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doParallel package installed
ncores <- 15

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/updated_Nov_5_2014/Fulldom_hires_netcdf_file.nc' 

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/geo_em.d02.nc' 

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 10

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/DOMAIN/urg_MASKS_NEWDOMAIN.Rdata' 

## Specify whether the model run used NHD reach-based routing (otherwise gridded routing assumed)
# If TRUE, mask file should also contain rtLinks dataframe.
reachRting <- FALSE

## Temp directory to write intermediate files
tmpDir <- '/glade/scratch/adugger'


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- NULL 

## Path to SNOTEl data .Rdata file
SNOfile <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata'

## Path to meteorological station data .Rdata file
METfile <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/MET/met_URG_NEW.Rdata'

## Path to streamflow data .Rdata file
STRfile <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/obsStrData.Rdata'
 

################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
	modPathList <- c('/glade/p/ral/RHAP/gochis/Col_SWCol_URG_UppGunn_Animas/test_debug',
                        '/glade/p/ral/RHAP/karsten/Col_SWCol_URG_UppGunn_Animas_NSSL') 

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('NLDAS2-Downscaled', 'NSSL') 

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_modelreads_RETEST.Rdata' 
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        modAppend <- TRUE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- FALSE
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- TRUE

		# Channel routing
		readChrtout <- FALSE
			# Read only links with gages?
			readChrtout_GAGES <- FALSE
			# Read specified subset? Provide object with link and site_no columns
			readLink2gage <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/link2gage_bigbigrivs.txt", sep="\t", header=TRUE, colClasses=c("integer","character"))

		# Snotel sites
		readSnoLdasout <- TRUE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- TRUE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- TRUE
	varsLdasoutNFIE <- FALSE
	varsLdasoutIOC0 <- FALSE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- NULL 
	readModEnd <- NULL 


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/scratch/zhangyx/WRF-Hydro/RioGrande/NLDAS2.data') 

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
	forcTagList <- c('NLDAS2-Downscaled')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_forcingreads_RETEST.Rdata'
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
        readForcEnd <- as.POSIXct("2013-12-31", format="%Y-%m-%d", tz="UTC")


############# Model Performance Stats #############

## Calculate stats?
calcStats <- FALSE

	## Calculate streamflow performance stats?
	strProc <- TRUE
		# Read specified subset? Provide object with link and site_no columns
		statsLink2gage <- NULL 
		# Calculate daily stats?
		strProcDaily <- FALSE

	## Calculate SNOTEL performance stats?
	snoProc <- TRUE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- TRUE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
        modReadFileIn <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_modelreads_RETEST.Rdata'
        forcReadFileIn <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_forcingreads.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_urg_stats_RETEST.Rdata'

	# Range dates for main stats
        stdate_stats <- NULL
        enddate_stats <- as.POSIXct("2015-09-30 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
        stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
        enddate_stats_sub <- as.POSIXct("2015-09-30 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_PLOTS_RETEST'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- TRUE

## Create HTML files?
writeHtml <- TRUE

## If TRUE, specify output directory
writePlotDir <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/151111_PLOTS_RETEST'

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- TRUE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- as.POSIXct("2015-04-01", format="%Y-%m-%d", tz="UTC")

		# Specify end date
		accflowEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

	## Generate hydrographs?
	hydroPlot <- TRUE

        	# Specify which run tags to plot
        	hydroTags <- NULL 
 
        	# Specify start date
        	hydroStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	hydroEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

	## Generate accumulated precip plots?
	accprecipPlot <- TRUE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC") 
        
        	# Specify end date
        	accprecipEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- TRUE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	flowsweEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

        ## Generate Streamflow and Basin-mean LSM Runoff plots?
        flowlsmPlot <- FALSE

                # Specify which run tags to plot
                flowlsmTags <- NULL

                # Specify start date
                flowlsmStartDate <- NULL

                # Specify end date
                flowlsmEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- TRUE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

        	# Specify end date
        	sweEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

        ## Generate MET station plots?
        metPlot <- TRUE

                # Specify which run tags to plot
                metTags <- NULL

                # Specify start date
                metStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

                # Specify end date
                metEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")


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

        ## Generate Ameriflux ET error maps?
        amfetErrMap <- FALSE

                # Specify which run tags to plot
                amfetErrTags <- NULL

                # Specify which run seasons to plot
                amfetErrSeas <- NULL

        ## Generate Ameriflux ET correlation maps?
        amfetCorrMap <- FALSE

                # Specify which run tags to plot
                amfetCorrTags <- NULL

                # Specify which run seasons to plot
                amfetCorrSeas <- NULL

	## Include summary stats tables?
	statsMapTables <- TRUE


###########################################################################################
## RUN (do not change anything below this line)

source("run_CONFIG.R")

