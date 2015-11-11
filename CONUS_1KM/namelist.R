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
AMFfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/AMF/obs_AMF_1998_current.Rdata" 

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
STRfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/gageData_BIGRIVSAMPLE.Rdata"


################ Model Output Reads ###############

## Read model output?
readMod <- TRUE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
        modPathList <- c('/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/NHDPLUS_Run/output_new/cold_start_no_terr_rtg')

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
        modTagList <- c('SPINUP2013_GWPT')

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_spinup2013_modelout.Rdata'
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
	snoProc <- TRUE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- TRUE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_spinup2013_modelout.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_spinup2013_stats.Rdata'

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- NULL

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2014-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2014-10-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/spinup2013_PLOTS'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- FALSE

## Create HTML files?
writeHtml <- TRUE

## If TRUE, specify output directory
writePlotDir <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/spinup2013_PLOTS'

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- TRUE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- as.POSIXct("2014-04-01", format="%Y-%m-%d", tz="UTC")

		# Specify end date
		accflowEndDate <- NULL

	## Generate hydrographs?
	hydroPlot <- TRUE

        	# Specify which run tags to plot
        	hydroTags <- NULL
        
        	# Specify start date
        	hydroStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	hydroEndDate <- NULL

	## Generate accumulated precip plots?
	accprecipPlot <- TRUE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- TRUE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	flowsweEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- TRUE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")

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

source("run_CONFIG.R")

