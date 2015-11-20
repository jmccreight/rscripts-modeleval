###################################################
##  Main Control Script to process model output  ##
###################################################

###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)
library(data.table)

load(maskFile)
source("util_FUNC.R")


# Model Reads 
if (readMod | readForc) {
        source("read_MODELOUT.R")
}

# Obs
if (calcStats | createPlots) {
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
		if (exists("link2gage.man")) {
			obsStrData <- plyr::join(obsStrData, link2gage.man, by="site_no")
		}
		rm(obsStrData_FINAL, obsStrMeta_FINAL, obsStrData_TMP, obsStrMeta_TMP)
	}
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
	if (metProc) {
        	if (is.null(forcReadFileOut)) {
                	if (file.exists(forcReadFileIn)) {
                        	load(forcReadFileIn)
                	}
        	} else {
                	if (is.null(forcReadFileIn)) {
                        	if (file.exists(forcReadFileOut)) {
                                	load(forcReadFileOut)
                        	}
                	} else {
                        	if (file.exists(forcReadFileIn)) {
                                	load(forcReadFileIn)
                        	}
                	}
        	}
	}
	source("calc_PERFSTATS.R")
}

# Plots
if (createPlots) {
	if (accflowPlot | hydroPlot | accprecipPlot | swePlot | 
			strBiasMap | strCorrMap | 
			snosweErrMap | snoprecipErrMap) {
        	message("Generating plots")
		#load(statsFileOut)
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
	}
	if (metPlot) {
                if (is.null(forcReadFileOut)) {
                        if (file.exists(forcReadFileIn)) {
                                load(forcReadFileIn)
                        }
                } else {
                        if (is.null(forcReadFileIn)) {
                                if (file.exists(forcReadFileOut)) {
                                        load(forcReadFileOut)
                                }
                        } else {
                                if (file.exists(forcReadFileIn)) {
                                        load(forcReadFileIn)
                                }
                        }
                }
        }
        source("calc_PLOTS.R")
}

# EXIT
quit("no")

