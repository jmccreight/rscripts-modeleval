###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doMC package installed
ncores <- 16

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/updated_Nov_5_2014/Fulldom_hires_netcdf_file.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/geo_em.d02.nc'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 10

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/DOMAIN/urg_MASKS.Rdata'


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- NULL 

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG.Rdata" 

## Path to meteorological station data .Rdata file
METfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/MET/met_URG.Rdata"

## Path to streamflow data .Rdata file
STRfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG.Rdata"

## Do these need updated? TRUE/FALSE ...to be worked on later


################ Model Output Reads ###############

## Read model output?
readMod <- TRUE

## If TRUE, specify the following to read in model output:

	# Specify the model run output directory or directories
	modPathList <- c('/glade/p/ral/RHAP/adugger/Upper_RioGrande/RUN.NEWMP/OUTPUT_NLDAS_SPINUP',
			'/glade/p/ral/RHAP/adugger/Upper_RioGrande/RUN.NEWMP/OUTPUT_NLDASDWNSC_SPINUP')

	# Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('su2013_15_NLDAS_newmodel', 
			'su2013_15_NLDASdwnsc_newmodel')

	# Specify the output .Rdata file to create
	modReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_spinup_eval_raw.Rdata'
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
	readStart <- NULL
	readEnd <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/scratch/zhangyx/WRF-Hydro/RioGrande/NLDAS2.data')

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
        forcTagList <- c('nldas2dwnsc')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_forcing_raw.Rdata'
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
        startForc <- NULL
        endForc <- NULL


############# Model Performance Stats #############

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
	modReadFileIn <- NULL

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_stats_all.Rdata'

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/EVAL_151013'



################### Plotting ######################
#(to be added after everything else works)



###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)

load(maskFile)
source("util_FUNC.R")

# Multi-core
if (ncores>1) {
	library(doMC)
	registerDoMC(ncores)
	parallelFlag <- TRUE
}

# Obs
if (!is.null(AMFfile)) {
	if (file.exists(AMFfile)) {
		load(AMFfile)
		#stop(paste("Ameriflux obs file specified but does not exist:", AMFfile))
	}
}
if (!is.null(SNOfile)) {
        if (file.exists(SNOfile)) {
                load(SNOfile)
                #stop(paste("SNOTEL obs file specified but does not exist:", SNOfile))
        }
}
if (!is.null(METfile)) {
        if (file.exists(METfile)) {
                load(METfile)
                #stop(paste("MET obs file specified but does not exist:", METfile))
        }
}
if (!is.null(STRfile)) {
        if (file.exists(STRfile)) {
                load(STRfile)
                #stop(paste("Streamflow obs file specified but does not exist:", STRfile))
        }
}

# Model Reads 
if (readMod | readForc) {
	source("read_MODELOUT.R")
}

# Model Output processing
if (strProc | snoProc | amfProc | metProc) {
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

