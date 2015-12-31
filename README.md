# rscripts-modeleval

WRF_Hydro model evaluation R scripts.

# Main control script

## namelist.R

This file is the main argument definition file, and it links to the other "common" scripts (the general idea is that the namelist would change for each project, but the processing scripts would mostly be shared). To run the tool, go into namelist.R and update the paths to your various files (run output dirs, obs pathnames, where you want to write files, etc.). Then there are various flags you can turn on and off for what you want to run. Most important are the run output directories, and each output directory should have an associated "tag" that the script will use to ID those outputs (and it will use that tag in plot labels and plot filenames, so be cognizant of tag length, spaces, etc.).

General workflow is to do the model reads first, then if that works successfully run the stats and plots options. The model reads are the most tedious, so depending on the number of files and model output directories you might want to chunk things up. You can turn on the "append" option to keep writing to the same Rdata workspace if you do need to chunk.

Note that there are some upper level boolean flags that turn off entire sections (generally denoted by indents) so it is easy to turn on/off reads, stats, and plots without turning off all of the sub-options.


# Observation data

## Streamflow

The streamflow R dataset should contain at least two dataframes called obsStrData and obsStrMeta. It will look for the following fields in particular:

obsStrData:
site_no - gage ID
POSIXct - observation timestep in POSIXct
q_cms - streamflow in m3/sec

obsStrMeta:
site_no - gage ID
site_name - gage name
area_sqmi - gage drainage area in square miles

These are the current minimum requirements, but other fields can be included for potential future use.

If you don't want to rename your headers to match those listed above, you can create a header "map" file for each (obsStrData.map, obsStrMeta.map), which is just a list that matches the header name that the tool wants with your data header name.

See: /glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/obsStrData_GAGESII_2010_2014_DV.Rdata for an example with USGS data.

## SNOTEL

## Ameriflux


# Outputs

modChrtout (data table) -  If you do LDASOUT reads those are lists of data tables. Stats are output as data frames.
modChrtout.d (data table) - Daily aggregation of modChrtout (by UTC day)
modLdasout - List of data tables, including the native model output timestep (modLdasout[["native"]]), UTC day (modLdasout[["utcday"]]), SNOTEL (PST) day (modLdasout[["snoday"]]), and UTC month (modLdasout[["utcmonth"]])
stats_ (data frame) - Statistics for each spatial unit (e.g., basin/gage, cell/station)


# Utilities

## Process Masks

Create properly formatted mask objects:
* Update the parameters in process_MASKS.R:
	* Optional info to process existing mask files
	* Pathname to the high-res grid and aggregation factor if creating from scratch
	* Output R dataset filename
* The script will create a number of mask objects:
	* msk.List - the basin mask matrices
	* msk.areaList - the total area of each basin
	* msk.nameList
	* msk.minInds
	* msk.maxInds
	* msk.countInds
	where mskgeo are masks at the geogrid resolution and maskhyd are masks at the high-res routing resolution
* The script also stores any relevant lookup tables (that you must manually specify)

## Process Observations

Download and prep observations. See:
* prep_AMF.R
* prep_SNOTEL.R
* prep_MET.R
* prep_STR.R

# EXAMPLE

```sh
###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doParallel package installed
ncores <- 15
```
Specify the number of cores you want to use. Assumes you have doParallel installed. The tool will initiate a doParallel session and setup a cluster when needed, then kill hte cluster after processing is complete.
```sh
## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/Fulldom_hires_netcdf_file_1km.nc'

## Specify the low-resolution geogrid file
geoFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/geo_em.d01.nc.conus_1km'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 1
```
This is where you specify the domain files used in your run. These are only used for reference, dimensions, etc. so exact content is not critical.
```sh
## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/gagesII_MASKS.Rdata'
```
This spcifies the location of the R dataset containing the mask information. This dataset is created by the process_MASK.R script. Mask options include: basin masks to match gages, cell indices to match observation points.
```sh
## Specify whether the model run used NHD reach-based routing (otherwise gridded routing assumed)
# If TRUE, mask file should also contain rtLinks dataframe.
reachRting <- TRUE
```
Specify whether or not your model uses reach-based routing. If it does, the model reads will rely on CHRTOUT files and the tool will look for a route link dataframe in your mask file. If not, the model reads will rely on Frxstout for streamflow output and the route link file is not needed.
```sh
## Temp directory to write intermediate files
tmpDir <- '/glade/scratch/adugger'
```
Provide a directory where the tool can write temp files. The tool will periodically write out an interim R dataset, particularly during model reads, so that a crash in one section does not necessarily mean you lose all data. These temp datasets are NOT deleted when the tool finishes, so you should periodically clean these up manually.

```sh
################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/AMF/obs_AMF_1998_current.Rdata"

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
STRfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/USGS/obsStrData_GAGESII_2010_2014_DV.Rdata"
```
Specify the paths to the R datasets containing various observation datasets that you want to do statistical analysis or plots using. These should follow the format specified above in the "OBSERVATIONS" section.

```sh
################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
        modPathList <- c('/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/NHDPLUS_Run_5yr_no_terr_rtg',
                        '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/NHDPLUS_Run/')

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
        modTagList <- c('SPINUP5YR_No_Routing',
                        'SPINUP5YR_All_Routing')

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/151217_conus_bigrivs_su2010allrt_modelout_STR.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        modAppend <- FALSE
```




# OLDER INSTRUCTIONS

## Process Model Output

1. Read in model output	
* Update the parameters in the urg_modelreads_ALL.R:
	* Path to the model output directory
	* Path to the forcing data (if needed)
	* Pathname for the output R dataset
	* Suffix to add to the objects to identify them
	* Options to run
* Run the urg_modelreads_ALL.R script: ./jobIS_reads.sh
* The script will create a new R dataset with "raw" data objects. Possible objects (depending on which options you choose):
	* modFrxstout - streamflow time series
	* modGwout - groundwater outflow time series
	* modLdasin - climate forcings
	* modLdasout - LSM output
	* Objects will be appended with "_BAS" for basin means and "_SNO" for point values.
* Import these new "raw" dataframes into the master "raw" dataframe containing output from all model runs.

2. Process the model output by basin:
* Update the parameters in the urg_process_BASIN.R:
	* Pathname to the raw dataset from step 1
	* Pathname for the output R dataset
	* Suffixes to process (all or some subset of those from step 1)
	* Dates that each of the model runs stopped putting out meaningful data (in case the model ran out with static forcings)
	* Start date to restrict the output data to (applies to all models)
	* Start and end dates to restrict the full stats calcs
	* Start and end dates to restrict the subset stats calcs
* Run the urg_process_BASIN.R script: ./jobIS_procBAS.sh
* The script will create a new R dataset with processed data objects at the basin scale. The objects will be same as in step 1 with new fields added. New objects include:
	* stats_str - summary statistics for each model run
	* stats_str_all - combined object with summary statistics for all runs

3. Process the model output by station:
* Update the parameters in the urg_process_SNOTEL.R:
	* Pathname to the raw dataset from step 1
	* Pathname for the output R dataset
	* Suffixes to process (all or some subset of those from step 1)
	* Dates that each of the model runs stopped putting out meaningful data (in case the model ran out with static forcings)
	* Start date to restrict the output data to (applies to all models)
	* Start and end dates to restrict the full stats calcs
	* Start and end dates to restrict the subset stats calcs
	* Options to run
* Run the urg_process_SNOTEL.R script: ./jobIS_procSNO.sh
* The script will create a new R dataset with processed data objects at the point station scale. The objects will be same as in step 1 with new fields added. New objects include:
	* .metd - objects as above but aggregated to the MET station day
	* .snod - objects as above but aggregated to the SNOTEL station day
	* stats_met_all - combined summary statistics for each model run compared to the MET stations
	* stats_sno_all - combined summary statistics for each model run compared to the SNOTEL stations

4. Create plots:
* Update the parameters in urg_plots_BASIN.R or urg_plots_SNOTEL.R:
	* Pathname to the processed R dataset
	* Pathname to the directory to output plots
	* Start and end date to restrict the plots
* Run the plot script to generate plot images. NOTE: plot function are stored in a separate R script file that is sourced by this script. This file specifies the loops over basins/sites, the plot function calls and the file exports.

