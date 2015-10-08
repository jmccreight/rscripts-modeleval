###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doMC package installed
ncores <- 16

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/nral0008/zhangyx/CONUS1km_LSMOnly_daily_snowmods/Fulldom_hires_netcdf_file.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/nral0008/zhangyx/CONUS1km_LSMOnly_daily_snowmods/geo_em.d01.nc.conus_1km'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 4


##################### Masks #######################
## Create mask objects from hydro domain basins?
createMask <- FALSE

## If FALSE, specify location of .Rdata file containing pre-processed mask objects
maskFileIn <- '/glade/p/ral/RHAP/adugger/CONUS_3km/ANALYSIS/conus_masks_ALL_NEW.Rdata'

## If TRUE, specify the following:

	# Select which masks/points to create:
		# Basin masks
		createBasMask <- TRUE
		# Ameriflux points
		createAmfMask <- TRUE
		# SNOTEL points
		createSnoMask <- TRUE
		# MET station points
		createMetMask <- FALSE
		# MET station sites (must include columns: id, lat, lon)
		metSites <- NULL

	# Specify the .Rdata file to create
	maskFileOut <- NULL


################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- "/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/AMF/amf_URG.Rdata" 

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
	modPathList <- c('/glade/p/nral0008/zhangyx/CONUS1km_LSMOnly_daily_snowmods')

	# Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('su2011_13_LSMonly_NLDASdwnsc_oldmodel')

	# Specify the output .Rdata file to create
	modReadFileOut <- 'conus1km_gagesII_raw.Rdata'

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- TRUE  
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- FALSE

		# Snotel sites
		readSnoLdasout <- TRUE

		# Ameriflux sites
		readAmfLdasout <- TRUE

		# MET sites
		readMetLdasout <- TRUE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- TRUE

	# Specify start and end dates if you do NOT want to read all files
	startRead <- NULL
	endRead <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- TRUE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPath <- '/glade/scratch/zhangyx/WRF-Hydro/CONUS1km.ForcingData'

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- 'conus1k_nldas2dwnsc.Rdata'

	# Select what aggregations/imports to run:

		# Basin means
		readBasinLdasin <- TRUE

		# SNOTEL sites
		readSnoLdasin <- TRUE

		# Ameriflux sites
		readAmfLdasin <- TRUE

		# MET sites
		readMetLdasin <- TRUE

        # Specify start and end dates if you do NOT want to read all files
        startForc <- NULL
        endForc <- NULL


############# Model Output Processing #############

## Process basin data?
basinProc <- TRUE

## Process point data?
pointProc <- TRUE

## If either is TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- NULL

	# Range dates to restrict analysis
	stdate <- NULL
	enddate <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for main stats
	stdate_stats <- NULL
	enddate_stats <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

	# Calculate stats?
	runStats <- TRUE
	writeStatsFile <- TRUE



################### Plotting ######################
#(to be added after everything else works)



###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)

# Multi-core
if (ncores>1) {
	library(doMC)
	registerDoMC(ncores)
}

# Masks
if (createMask) {
	source("create_MASKS.R")
} else {
	if (file.exists(maskFileIn)) {
		load(maskFileIn)
	} else if (readMod) {
		stop(paste("Model read requested but mask file does not exist:", maskFileIn, "Check path or change createMask flag to TRUE to create."))
	}
}

# Obs
if (!isnull(AMFfile)) {
	if (file.exists(AMFfile)) {
		load(AMFfile)
		#stop(paste("Ameriflux obs file specified but does not exist:", AMFfile))
	}
}
if (!isnull(SNOfile)) {
        if (file.exists(SNOfile)) {
                load(SNOfile)
                #stop(paste("SNOTEL obs file specified but does not exist:", SNOfile))
        }
}
if (!isnull(METfile)) {
        if (file.exists(METfile)) {
                load(METfile)
                #stop(paste("MET obs file specified but does not exist:", METfile))
        }
}
if (!isnull(STRfile)) {
        if (file.exists(STRfile)) {
                load(STRfile)
                #stop(paste("Streamflow obs file specified but does not exist:", STRfile))
        }
}

# Model Reads 
if (readMod) {
	source("read_MODELOUT.R")
} else {
	if (file.exists(modReadFileIn)) {
		load(modReadFileIn)
		#stop(paste("Model read file specified but does not exist:", modReadFileIn))
	}
}

# Model Output processing
if (basinProc | pointProc) {
	if (isnull(modReadFileOut)) {
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
}
if (basinProc) {source("process_BASIN.R")}
if (pointProc) {source("process_PTS.R")}
