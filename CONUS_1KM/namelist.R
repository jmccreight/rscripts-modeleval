###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doParallel package installed
ncores <- 15

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/Fulldom_hires_netcdf_file_1km.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/geo_em.d01.nc.conus_1km'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 1

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/gagesII_MASKS.Rdata'

## Specify whether the model run used NHD reach-based routing (otherwise gridded routing assumed)
# If TRUE, mask file should also contain rtLinks dataframe.
reachRting <- TRUE

## Temp directory to write intermediate files
tmpDir <- '/glade/scratch/adugger'


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/AMF/obs_AMF_1998_current.Rdata" 

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
STRfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/obsStrData_GAGESII_2010_2014_DV.Rdata"

################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
	modPathList <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/IOC_calib_runs/terr_rtg/v1.1_w_calib_no_oCONUS_no_res_new_route_link_Dec_26_2015'

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('SPINUP5YR_Full_Routing_v1.1')

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_GAGESII/151231_conus_bigrivs_su2010v11_modelout_STR_UPDATE.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        modAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- FALSE
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- FALSE

		# Channel routing
		readChrtout <- TRUE
			# Read only links with gages?
			readChrtout_GAGES <- FALSE
			# Read specified subset? Provide object with link and site_no columns
			readLink2gage <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/link2gage_gagesII.txt", sep="\t", header=TRUE, colClasses=c("integer", "character"))

		# Snotel sites
		readSnoLdasout <- FALSE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- FALSE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- FALSE
	varsLdasoutNFIE <- FALSE
	varsLdasoutIOC0 <- TRUE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- as.POSIXct("2010-01-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	readModEnd <- as.POSIXct("2015-12-31 23:59", format="%Y-%m-%d %H:%M", tz="UTC")


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/forcing/') 

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
	forcTagList <- c('NLDAS2-Downscaled')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_FORCINGS/conus_nldas_forcings_2010.Rdata'
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
        readForcStart <- as.POSIXct("2010-01-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
        readForcEnd <- as.POSIXct("2010-12-31 23:59", format="%Y-%m-%d %H:%M", tz="UTC")


############# Model Performance Stats #############

## Calculate stats?
calcStats <- FALSE

	## Calculate streamflow performance stats?
	strProc <- TRUE
		# Read specified subset? Provide object with link and site_no columns
                statsLink2gage <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/link2gage_gagesII_CORRECTED.txt",
                                        sep="\t", header=TRUE, colClasses=c("integer","character"))
                # Calculate daily stats?
                strProcDaily <- TRUE

	## Calculate SNOTEL performance stats?
	snoProc <- FALSE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- FALSE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_GAGESII/160106_conus_gagesII_su2010v11_modelout_STR.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_GAGESII/160106_conus_gagesII_su2010v11_stats_STR.Rdata'

	# Range dates for main stats
	stdate_stats <- as.POSIXct("2011-01-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats <- as.POSIXct("2013-12-31 23:59", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2011-10-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2013-09-30 23:59", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_GAGESII/160106_gagesII_su2010v11_PLOTS_ALL'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- TRUE

	## Create HTML files?
	writeHtml <- TRUE

	## If TRUE, specify output directory
	writePlotDir <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/OUT_GAGESII/160106_gagesII_su2010v11_PLOTS_ALL'

	## Plot specified subset? Provide object with link and site_no columns. HYDRO PLOTS ONLY!
        plotLink2gage <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/link2gage_gagesII_CORRECTED.txt",
                                        sep="\t", header=TRUE, colClasses=c("integer","character"))

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- FALSE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- as.POSIXct("2014-04-01", format="%Y-%m-%d", tz="UTC")

		# Specify end date
		accflowEndDate <- NULL

	## Generate hydrographs?
	hydroPlot <- FALSE

        	# Specify which run tags to plot
        	hydroTags <- c("SPINUP5YR_No_Routing_v1.1", "SPINUP5YR_Full_Routing_v1.1")

        	# Specify start date
        	hydroStartDate <- as.POSIXct("2011-01-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	hydroEndDate <- as.POSIXct("2012-12-31", format="%Y-%m-%d", tz="UTC")

                # Plot daily values?
                hydroPlotDaily <- TRUE

	## Generate accumulated precip plots?
	accprecipPlot <- FALSE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- FALSE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	flowsweEndDate <- NULL

        ## Generate Streamflow and Basin-mean LSM Runoff plots?
        flowlsmPlot <- FALSE

                # Specify which run tags to plot
                flowlsmTags <- NULL

                # Specify start date
                flowlsmStartDate <- NULL

                # Specify end date
                flowlsmEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- FALSE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")

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

        ## Threshold for "completeness" of sites to include. Multiplier on the max n in the set.
        #  Set to 0 to plot all.
        nThresh <- 0.75

        ## Trusted gages table (if used, otherwise NULL). 
        #  The table should have a "site_no" column and a "fractPerfect" column with the "trust" metric (0-1).
        trustGages <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/trustGages.txt",
                                        sep="\t", header=TRUE, colClasses=c("character", "numeric"))
        # Threshold for "trusted" gages. Set to 0 to plot all.
        trustThresh <- 0.75

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

        ## Generate SNOTEL SWE error maps?
        amfetErrMap <- FALSE

                # Specify which run tags to plot
                amfetErrTags <- NULL

                # Specify which run seasons to plot
                amfetErrSeas <- NULL

        ## Generate SNOTEL SWE error maps?
        amfetCorrMap <- FALSE
        
                # Specify which run tags to plot
                amfetCorrTags <- NULL

                # Specify which run seasons to plot
                amfetCorrSeas <- NULL

	## Include summary stats tables?
	statsMapTables <- FALSE


###########################################################################################
## RUN (do not change anything below this line)

source("run_CONFIG.R")

