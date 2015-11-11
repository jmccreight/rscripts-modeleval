###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doMC package installed
ncores <- 15

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/updated_Nov_5_2014/Fulldom_hires_netcdf_file.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/geo_em.d02.nc'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 10

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/DOMAIN/urg_MASKS_NEWDOMAIN.Rdata'


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- NULL 

## Path to SNOTEl data .Rdata file
#SNOfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG_NEW.Rdata" 
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/MET/met_URG_NEW.Rdata"

## Path to streamflow data .Rdata file
#STRfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG_NEW.Rdata"
STRfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/obsStrData.Rdata"


################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

	# Specify the model run output directory or directories
	modPathList <- c('/glade/p/ral/RHAP/gochis/Col_SWCol_URG_UppGunn_Animas/results/spinup_WY_2014_thru_2015_cold_start_full_rtng',
			'/glade/scratch/karsten/Col_SWCol_URG_UppGunn_Animas_NSSL')

	# Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('NLDAS2-Downscaled', 'NSSL')

	# Specify the output .Rdata file to create
	modReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_urg_modelreads.Rdata'
	# Append to existing file? FALSE = create new file (or overwrite existing!)
	modAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- TRUE  
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- TRUE

		# Snotel sites
		readSnoLdasout <- TRUE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- TRUE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- TRUE
	varsLdasoutNFIE <- FALSE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- NULL
	readModEnd <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- TRUE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/scratch/zhangyx/WRF-Hydro/RioGrande/NLDAS2.data')

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
	forcTagList <- c('NLDAS2-Downscaled')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_urg_forcingreads.Rdata'
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
	modReadFileIn <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_urg_modelreads.Rdata'
	forcReadFileIn <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_urg_forcingreads.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_urg_stats.Rdata'

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_PLOTS'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- FALSE

## Create HTML files?
writeHtml <- TRUE

## If TRUE, specify output directory
writePlotDir <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/aTest_PLOTS'

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
	accprecipPlot <- TRUE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- TRUE

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

        ## Generate MET station plots?
        metPlot <- TRUE

                # Specify which run tags to plot
                metTags <- NULL

                # Specify start date
                metStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

                # Specify end date
                metEndDate <- NULL


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

source("run_CONFIG.R")

