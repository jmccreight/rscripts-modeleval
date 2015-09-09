# rscripts-modeleval

WRF_Hydro model evaluation R scripts.

# Main control script

TO BE DETERMINED

# Process Masks

Create properly formatted mask objects:
* Update the parameters in urg_process_MASKS.R:
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

# Process Observations

Download and prep observations. See:
* urg_prep_AMF.R
* urg_prep_SNOTEL.R
* urg_prep_MET.R
* urg_prep_STR.R

# Process Model Output

1. Read in model output	
* Update the parameters in the urg_modelreads_ALL.R:
	* Path to the model output directory
	* Path to the forcing data (if needed)
	* Pathname for the output R dataset
	* Suffix to add to the objects to identify them
	* Options to run
* Run the urg_modelreads_ALL.R script.
* The script will create a new R dataset with "raw" data objects. Possible objects (depending on which options you choose):
	* modFrxstout - streamflow time series
	* modGwout - groundwater outflow time series
	* modLdasin - climate forcings
	* modLdasout - LSM output
	* Objects will be appended with "_BAS" for basin means and "_SNO" for point values.

2. Process the model output by basin:
* Update the parameters in the urg_process_BASIN.R:
	* Pathname to the raw dataset from step 1
	* Pathname for the output R dataset
	* Suffixes to process (all or some subset of those from step 1)
	* Dates that each of the model runs stopped putting out meaningful data (in case the model ran out with static forcings)
	* Start date to restrict the output data to (applies to all models)
	* Start and end dates to restrict the full stats calcs
	* Start and end dates to restrict the subset stats calcs
* Feed the R dataset created in step 1 into this processing function.
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
* Feed the R dataset created in step 1 into this processing function.
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

