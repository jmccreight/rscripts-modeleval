###################################################
##         CALCULATE PERFORMANCE STATS           ##
###################################################


################## General Setup ##################

saveList <- c()
if (writeStatsFile) {
	dir.create(writeDir, showWarnings = FALSE)
}

## ------------------------------------------------------------------------
# Setup processing functions

CalcStrStats <- function(modDf, obsDf, stid2gageList.=stid2gageList, 
                      stdate=NULL, enddate=NULL) {
  sites <- names(stid2gageList.)
  results <- data.frame()
  results <- foreach(n=1:length(sites), .combine=rbind, .inorder=FALSE, .errorhandling='remove') %dopar% {
  out <- tryCatch(suppressWarnings( CalcModPerfMulti( subset(modDf, modDf$st_id==sites[n]), 
                                                      subset(obsDf, obsDf$Station==stid2gageList[[sites[n]]]), 
                                                      flxCol.obs="q_cms_adj", flxCol.mod="q_cms",
                                                      stdate=stdate,
                                                      enddate=enddate) ), 
                  error=function(cond) {message(cond); return(NA)})
  if ( !is.na(out) ) {
      out$st_id <- sites[n]
      out$STAID <- stid2gageList.[[sites[n]]]
      out
      }
  }
  #results<-plyr::join(results, stid2gageList, by="site_id")
  results[results=="Inf"]<-NA
  results[results=="-Inf"]<-NA
  results
}

CalcVarStats <- function(modDf, siteDf, obsDf,
                        stdate=NULL, enddate=NULL,
                        flxCol.obs, flxCol.mod,
                        idCol.obs="site_id", idCol.mod="statArg",
			overwriteDate=FALSE) {
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    modDf <- subset(modDf, modDf$POSIXct >= stdate & modDf$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    modDf <- subset(modDf, modDf$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    modDf <- subset(modDf, modDf$POSIXct <= enddate)
  }
  # Set date if flagged
  if (overwriteDate) {
  	modDf$POSIXct <- as.POSIXct(CalcDateTrunc(modDf$POSIXct, timeZone = "UTC"), tz="UTC")
  	obsDf$POSIXct <- as.POSIXct(CalcDateTrunc(obsDf$POSIXct, timeZone = "UTC"), tz="UTC")  
  }
  # Run stats
  sites <- unique(siteDf[,idCol.obs])
  results <- data.frame()
  results <- foreach(n=1:length(sites), .combine=plyr::rbind.fill, .inorder=FALSE, .errorhandling='remove') %dopar% {
  	out <- tryCatch(suppressWarnings( CalcModPerfMulti( subset(modDf, modDf[,idCol.mod]==sites[n]),
                                                      subset(obsDf, obsDf[,idCol.obs]==sites[n]),
                                                      flxCol.obs=flxCol.obs, flxCol.mod=flxCol.mod) ),
                                    error=function(cond) {message(cond); return(NA)})
       if ( !is.na(out) ) {
         out$site_id <- sites[n]
         out
         }
       }
#  if (idCol.mod != idCol.obs) {
#       siteDf[,idCol.mod] <- siteDf[,idCol.obs]
#       }
  results<-plyr::join(results, siteDf, by=idCol.obs)

  results[results=="Inf"]<-NA
  results[results=="-Inf"]<-NA
  return(results)
}

## ------------------------------------------------------------------------
# Calculate Streamflow Stats
if (strProc) {
	# Initialize
	stats_str <- data.frame()
	runTagList <- unique(modFrxstout$tag)
	# Loop through runs
	for (j in 1:length(runTagList)) {
		runTag <- runTagList[j]
		# Subset
		modFrxstout_tmp <- subset(modFrxstout, modFrxstout$tag==runTag)
		# Stats
		results <- CalcStrStats(modFrxstout_tmp, obsStr, stid2gageList, stdate=stdate_stats, enddate=enddate_stats)
		results$tag <- runTag
		results$seas <- "Full"
		stats_str <- rbind(stats_str, results)
		results <- CalcStrStats(modFrxstout_tmp, obsStr, stid2gageList, stdate=stdate_stats_sub, enddate=enddate_stats_sub)
		results$tag <- runTag
		results$seas <- "Sub"
		stats_str <- rbind(stats_str, results)
	}
	# Output
	if (writeStatsFile) {
		# Change NAs to large negative for QGIS so data type is not affected
		stats_str_tmp <- stats_str
		stats_str_tmp[is.na(stats_str_tmp)]<-(-1e+30)
		write.table(stats_str_tmp, file=paste0(writeDir, "/stats_str.txt"), sep="\t", row.names=FALSE)
	}
saveList <- c(saveList, stats_str)
}

## -----------------------------------------------------------------------
# Calculate SNOTEL Stats
if (snoProc) {
	# Forcing stats
	if (exists("modLdasin_SNO")) {
	        # Initialize
        	stats_ldasin_sno <- data.frame()
        	runTagList <- unique(modLdasin_SNO[["native"]]$tag)
		modLdasin_tmp.snod <- modLdasin_SNO[["snoday"]]
        	# Loop through runs
        	for (j in 1:length(runTagList)) {
                	runTag <- runTagList[j]
  			# Subset
			modLdasin_tmp.snod <- subset(modLdasin_tmp.snod, modLdasin_tmp.snod$tag==runTag)
                        # SNOTEL (daily)
                        # Full run
                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="Tavg_K", flxCol.mod="T2D_mean")
                        results$tag <- runTag
                        results$var <- "Tmean"
                        results$seas <- "Full"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)

                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="Tmin_K", flxCol.mod="T2D_min")
                        results$tag <- runTag
                        results$var <- "Tmin"
                        results$seas <- "Full"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)

                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="Tmax_K", flxCol.mod="T2D_max")
                        results$tag <- runTag
                        results$var <- "Tmax"
                        results$seas <- "Full"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)

                        # Subset (e.g., spring)
                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tavg_K", flxCol.mod="T2D_mean")
                        results$tag <- runTag
                        results$var <- "Tmean"
                        results$seas <- "Sub"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)

                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tmin_K", flxCol.mod="T2D_min")
                        results$tag <- runTag
                        results$var <- "Tmin"
                        results$seas <- "Sub"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)

                        results <- CalcVarStats(modLdasin_tmp.snod, sno.sites, sno.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="Tmax_K", flxCol.mod="T2D_max")
                        results$tag <- runTag
                        results$var <- "Tmax"
                        results$seas <- "Sub"
                        stats_ldasin_sno <- plyr::rbind.fill(stats_ldasin_sno, results)
		} # end loop
        	# Output
        	if (writeStatsFile) {
                	# Change NAs to large negative for QGIS so data type is not affected
                	stats_ldasin_sno_tmp <- stats_ldasin_sno
                	stats_ldasin_sno_tmp[is.na(stats_ldasin_sno_tmp)]<-(-1e+30)
                	write.table(stats_ldasin_sno_tmp, file=paste0(writeDir, "/stats_ldasin_sno.txt"), sep="\t", row.names=FALSE)
        	}
		saveList <- c(saveList, stats_ldasin_sno)
	} # end if modLdasin

        # Output stats
        if (exists("modLdasout_SNO")) {
                # Initialize
                stats_ldasout_sno <- data.frame()
                runTagList <- unique(modLdasout_SNO[["native"]]$tag)
		if ("snoday" %in% names(modLdasout_SNO)) {
                	modLdasout_tmp.snod <- modLdasout_SNO[["snoday"]]
			overwriteDate_flag <- FALSE
                } else {
			modLdasout_tmp.snod <- modLdasout_SNO[["native"]]
			overwriteDate_flag <- TRUE
		}
		# Loop through runs
                for (j in 1:length(runTagList)) {
                        runTag <- runTagList[j]
                        # Subset
                        modLdasout_tmp.snod <- subset(modLdasout_tmp.snod, modLdasout_tmp.snod$tag==runTag)
                	# SNOTEL (daily)
                	# Full time period
                	results <- CalcVarStats(modLdasout_tmp.snod, sno.sites, sno.data, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="Prec_mm", flxCol.mod="DEL_ACCPRCP", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "Precip"
                	results$seas <- "Full"
                	stats_ldasout_sno <- plyr::rbind.fill(stats_ldasout_sno, results)

                	results <- CalcVarStats(modLdasout_tmp.snod, sno.sites, sno.data, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="SWE_mm", flxCol.mod="SNEQV_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SWE"
                	results$seas <- "Full"
                	stats_ldasout_sno <- plyr::rbind.fill(stats_ldasout_sno, results)

	                # Subset time period
        	        results <- CalcVarStats(modLdasout_tmp.snod, sno.sites, sno.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                	                flxCol.obs="Prec_mm", flxCol.mod="DEL_ACCPRCP", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "Precip"
                	results$seas <- "Sub"
                	stats_ldasout_sno <- plyr::rbind.fill(stats_ldasout_sno, results)

                	results <- CalcVarStats(modLdasout_tmp.snod, sno.sites, sno.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                       		        flxCol.obs="SWE_mm", flxCol.mod="SNEQV_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SWE"
                	results$seas <- "Sub"
                	stats_ldasout_sno <- plyr::rbind.fill(stats_ldasout_sno, results)
                } # end loop
                # Output
                if (writeStatsFile) {
                        # Change NAs to large negative for QGIS so data type is not affected
                        stats_ldasout_sno_tmp <- stats_ldasout_sno
                        stats_ldasout_sno_tmp[is.na(stats_ldasout_sno_tmp)]<-(-1e+30)
                        write.table(stats_ldasout_sno_tmp, file=paste0(writeDir, "/stats_ldasout_sno.txt"), sep="\t", row.names=FALSE)
                }
                saveList <- c(saveList, stats_ldasout_sno)
        } # end if modLdasout
} # end snoProc


## -----------------------------------------------------------------------
# Calculate MET Stats
if (metProc) {
        # Forcing stats
        if (exists("modLdasin_MET")) {
                # Initialize
                stats_ldasin_met <- data.frame()
                runTagList <- unique(modLdasin_MET[["native"]]$tag)
                modLdasin_tmp.meth <- modLdasin_MET[["native"]]
                # Loop through runs
                for (j in 1:length(runTagList)) {
                        runTag <- runTagList[j]
                        # Subset
                        modLdasin_tmp.meth <- subset(modLdasin_tmp.meth, modLdasin_tmp.meth$tag==runTag)
                        # MET (hourly)
                        # Full run
                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats, enddate=enddate_stats, 
					flxCol.obs="Tmean_K", flxCol.mod="T2D")
                        results$tag <- runTag
                        results$var <- "Temp"
                        results$seas <- "Full"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="RH_mean", flxCol.mod="RelHum")
                        results$tag <- runTag
                        results$var <- "RelHum"
                        results$seas <- "Full"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="SurfPressmean_Pa", flxCol.mod="PSFC")
                        results$tag <- runTag
                        results$var <- "SurfPress"
                        results$seas <- "Full"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="Wind_mean", flxCol.mod="Wind")
                        results$tag <- runTag
                        results$var <- "Wind"
                        results$seas <- "Full"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="SWRad_mean", flxCol.mod="SWDOWN")
                        results$tag <- runTag
                        results$var <- "SWdown"
                        results$seas <- "Full"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        # Subset (e.g., spring)
                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tmean_K", flxCol.mod="T2D")
                        results$tag <- runTag
                        results$var <- "Temp"
                        results$seas <- "Sub"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="RH_mean", flxCol.mod="RelHum")
                        results$tag <- runTag
                        results$var <- "RelHum"
                        results$seas <- "Sub"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="SurfPressmean_Pa", flxCol.mod="PSFC")
                        results$tag <- runTag
                        results$var <- "SurfPress"
                        results$seas <- "Sub"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="Wind_mean", flxCol.mod="Wind")
                        results$tag <- runTag
                        results$var <- "Wind"
                        results$seas <- "Sub"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        results <- CalcVarStats(modLdasin_tmp.meth, met.sites, met.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="SWRad_mean", flxCol.mod="SWDOWN")
                        results$tag <- runTag
                        results$var <- "SWdown"
                        results$seas <- "Sub"
                        stats_ldasin_met <- plyr::rbind.fill(stats_ldasin_met, results)

                        } # end modldasin loop
                # Output
                if (writeStatsFile) {
                        # Change NAs to large negative for QGIS so data type is not affected
                        stats_ldasin_met_tmp <- stats_ldasin_met
                        stats_ldasin_met_tmp[is.na(stats_ldasin_met_tmp)]<-(-1e+30)
                        write.table(stats_ldasin_met_tmp, file=paste0(writeDir, "/stats_ldasin_met.txt"), sep="\t", row.names=FALSE)
                }
                saveList <- c(saveList, stats_ldasin_met)
        } # end if modLdasin

        # Output stats
        if (exists("modLdasout_MET")) {
                # Initialize
                stats_ldasout_met <- data.frame()
                runTagList <- unique(modLdasout_MET[["native"]]$tag)
		if ("utcday" %in% names(modLdasout_MET)) {
                        modLdasout_tmp.metd <- modLdasout_MET[["utcday"]]
                        overwriteDate_flag <- FALSE
                } else {
                        modLdasout_tmp.metd <- modLdasout_MET[["native"]]
                        overwriteDate_flag <- TRUE
                }
                # Loop through runs
                for (j in 1:length(runTagList)) {
                        runTag <- runTagList[j]
                        # Subset
                        modLdasout_tmp.metd <- subset(modLdasout_tmp.metd, modLdasout_tmp.metd$tag==runTag)
                        # MET (daily)
                        # Full time period
                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="PrecTot", flxCol.mod="DEL_ACCPRCP", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "Precip"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                      		          flxCol.obs="SnoDepmean_m", flxCol.mod="SNOWH_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SnowDepth"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="SM1_mean", flxCol.mod="SOIL_M1_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilM1"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="SM2_mean", flxCol.mod="SOIL_M2_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilM2"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="SoilT1mean_K", flxCol.mod="SOIL_T1_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilT1"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="SoilT2mean_K", flxCol.mod="SOIL_T2_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilT2"
                	results$seas <- "Full"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	# Subset time period
                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                	flxCol.obs="PrecTot", flxCol.mod="DEL_ACCPRCP", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "Precip"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                        	        flxCol.obs="SnoDepmean_m", flxCol.mod="SNOWH_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SnowDepth"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                        	        flxCol.obs="SM1_mean", flxCol.mod="SOIL_M1_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilM1"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                        	        flxCol.obs="SM2_mean", flxCol.mod="SOIL_M2_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilM2"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                        	        flxCol.obs="SoilT1mean_K", flxCol.mod="SOIL_T1_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilT1"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                	results <- CalcVarStats(modLdasout_tmp.metd, met.sites, met.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                        	        flxCol.obs="SoilT2mean_K", flxCol.mod="SOIL_T2_mean", overwriteDate=overwriteDate_flag)
                	results$tag <- runTag
                	results$var <- "SoilT2"
                	results$seas <- "Sub"
                	stats_ldasout_met <- plyr::rbind.fill(stats_ldasout_met, results)

                } # end loop
                # Output
                if (writeStatsFile) {
                        # Change NAs to large negative for QGIS so data type is not affected
                        stats_ldasout_met_tmp <- stats_ldasout_met
                        stats_ldasout_met_tmp[is.na(stats_ldasout_met_tmp)]<-(-1e+30)
                        write.table(stats_ldasout_met_tmp, file=paste0(writeDir, "/stats_ldasout_met.txt"), sep="\t", row.names=FALSE)
                }
                saveList <- c(saveList, stats_ldasout_met)
	} # end if modLdasin

} # end if metProc

## -----------------------------------------------------------------------
# Calculate AMF Stats
if (amfProc) {
        # Forcing stats
        if (exists("modLdasin_AMF")) {
                # Initialize
                stats_ldasin_amf <- data.frame()
                runTagList <- unique(modLdasin_AMF[["native"]]$tag)
                modLdasin_tmp.amfh <- modLdasin_AMF[["native"]]
                # Loop through runs
                for (j in 1:length(runTagList)) {
                        runTag <- runTagList[j]
                        # Subset
                        modLdasin_tmp.amfh <- subset(modLdasin_tmp.amfh, modLdasin_tmp.amfh$tag==runTag)
                        # AMF (hourly)
                        # Full run
                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="Rg", flxCol.mod="SWFORC")
                        results$tag <- runTag
                        results$var <- "SWdown"
                        results$seas <- "Full"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="Rgl", flxCol.mod="LWFORC")
                        results$tag <- runTag
                        results$var <- "LWdown"
                        results$seas <- "Full"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="RH", flxCol.mod="RelHum")
                        results$tag <- runTag
                        results$var <- "RelHum"
                        results$seas <- "Full"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="WS", flxCol.mod="Wind")
                        results$tag <- runTag
                        results$var <- "Wind"
                        results$seas <- "Full"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        # Subset (e.g., spring)
                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Rg", flxCol.mod="SWFORC")
                        results$tag <- runTag
                        results$var <- "SWdown"
                        results$seas <- "Sub"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="Rgl", flxCol.mod="LWFORC")
                        results$tag <- runTag
                        results$var <- "LWdown"
                        results$seas <- "Sub"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="RH", flxCol.mod="RelHum")
                        results$tag <- runTag
                        results$var <- "RelHum"
                        results$seas <- "Sub"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        results <- CalcVarStats(modLdasin_tmp.amfh, amf.sites, amf.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="WS", flxCol.mod="Wind")
                        results$tag <- runTag
                        results$var <- "Wind"
                        results$seas <- "Sub"
                        stats_ldasin_amf <- plyr::rbind.fill(stats_ldasin_amf, results)

                        } # end modldasin loop
                # Output
                if (writeStatsFile) {
                        # Change NAs to large negative for QGIS so data type is not affected
                        stats_ldasin_amf_tmp <- stats_ldasin_amf
                        stats_ldasin_amf_tmp[is.na(stats_ldasin_amf_tmp)]<-(-1e+30)
                        write.table(stats_ldasin_amf_tmp, file=paste0(writeDir, "/stats_ldasin_amf.txt"), sep="\t", row.names=FALSE)
                }   
                saveList <- c(saveList, stats_ldasin_amf)
        } # end if modLdasin

        # Output stats
        if (exists("modLdasout_AMF")) {
                # Initialize
                stats_ldasout_amf <- data.frame()
                runTagList <- unique(modLdasout_AMF[["native"]]$tag)
                if ("utcday" %in% names(modLdasout_AMF)) {
                        modLdasout_tmp.amf <- modLdasout_AMF[["utcday"]]
                        overwriteDate_flag <- FALSE
                } else {
                        modLdasout_tmp.amf <- modLdasout_AMF[["native"]]
                        overwriteDate_flag <- TRUE
                }   
                # Loop through runs
                for (j in 1:length(runTagList)) {
                        runTag <- runTagList[j]
                        # Subset
                        modLdasout_tmp.amf <- subset(modLdasout_tmp.amf, modLdasout_tmp.amf$tag==runTag)
                        # AMF (daily)
                        # Full time period
                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="RgNet", flxCol.mod="FSA", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "SWnet"
                        results$seas <- "Full"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                          flxCol.obs="RglNet", flxCol.mod="FIRAsurf", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "LWnet"
                        results$seas <- "Full"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="Rn", flxCol.mod="Rnet", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "Rnet"
                        results$seas <- "Full"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="LE", flxCol.mod="LH", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "LH"
                        results$seas <- "Full"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="H", flxCol.mod="HFX", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "SH"
                        results$seas <- "Full"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        # Subset time period
                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="RgNet", flxCol.mod="FSA", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "SWnet"
                        results$seas <- "Sub"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                          flxCol.obs="RglNet", flxCol.mod="FIRAsurf", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "LWnet"
                        results$seas <- "Sub"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Rn", flxCol.mod="Rnet", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "Rnet"
                        results$seas <- "Sub"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="LE", flxCol.mod="LH", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "LH"
                        results$seas <- "Sub"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                        results <- CalcVarStats(modLdasout_tmp.amf, amf.sites, amf.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="H", flxCol.mod="HFX", overwriteDate=overwriteDate_flag)
                        results$tag <- runTag
                        results$var <- "SH"
                        results$seas <- "Sub"
                        stats_ldasout_amf <- plyr::rbind.fill(stats_ldasout_amf, results)

                } # end loop
                # Output
                if (writeStatsFile) {
                        # Change NAs to large negative for QGIS so data type is not affected
                        stats_ldasout_amf_tmp <- stats_ldasout_amf
                        stats_ldasout_amf_tmp[is.na(stats_ldasout_amf_tmp)]<-(-1e+30)
                        write.table(stats_ldasout_amf_tmp, file=paste0(writeDir, "/stats_ldasout_amf.txt"), sep="\t", row.names=FALSE)
                }
                saveList <- c(saveList, stats_ldasout_amf)
        } # end if modLdasin

} # end if amfProc


## ------------------------------------------------------------------------
# Cleanup`
save(list=saveList, file=statsFileOut)

proc.time()
quit("no")

