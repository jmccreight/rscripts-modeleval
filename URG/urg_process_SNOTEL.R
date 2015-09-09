###################################################################################################
# Setup
# 
# Load the rwrfhydro package. 
## ------------------------------------------------------------------------
library("rwrfhydro")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_masks_NEW.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/MET/met_URG.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/AMF/amf_URG.Rdata")

# If you want to use R's multi-core capability (make sure  doMC is installed) specify the number 
# of cores.
## ------------------------------------------------------------------------
ncores <- 16
library(doMC)
registerDoMC(ncores)

# Where the raw data lives
#inImg <- 'urg_wy2015_NLDAS2dwnsc_fullrtng_SNO.Rdata'
#inImg <- 'urg_wy2015_NLDAS2dwnsc_NSSL_fullrtng_SNO.Rdata'
#inImg <- 'urg_wy2015_NLDASdwnsc_fullrtng_snowmod_SNO.Rdata'
#inImg <- 'urg_wy2015_NLDASdwnsc_NSSL_fullrtng_snowmod_SNO.Rdata'
#inImg <- 'urg_wy2015_NLDAS2dwnsc_SIMGM_BATSalb_fullrtng_ALL.Rdata'

inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_ALL_RAW.Rdata'

# Where to output the processed data
#outImg <- 'urg_wy2015_NLDAS2dwnsc_fullrtng_SNO_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_NSSL_fullrtng_SNO_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDASdwnsc_fullrtng_snowmod_SNO_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDASdwnsc_NSSL_fullrtng_snowmod_SNO_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_SIMGM_BATSalb_fullrtng_SNO_PROCESSED.Rdata'

outImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_SNOCLIM_PROCESSED.Rdata'

# Suffix for the output objects
objSuffixList <- c('_wy2015_NLDAS2dwnsc_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_fullrtng',
                   '_wy2015_NLDAS2dwnsc_snowmod_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng',
                   '_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng',
                   '_wy2015_NLDAS2dwnsc_SIMGM_BATSalb_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_nlcd11_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng',
                   '_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist1_canresist05_fullrtng',
                   '_wy2015_NLDAS2dwnsc_snowmod_mikerec_snowresist50_fullrtng')

stopDates <- c(as.POSIXct("2015-06-11 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-06-13 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-06-26 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-06-26 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-07-13 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-07-15 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-06-04 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-07-15 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-07-13 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC"))

# Range dates to restrict analysis
stdate <- NULL
#enddate <- NULL
enddate <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

# Range dates for main stats
stdate_stats <- as.POSIXct("2014-10-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
enddate_stats <- as.POSIXct("2015-05-31 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

# Range dates for seasonal stats (e.g., spring)
stdate_stats_sub <- as.POSIXct("2014-12-12 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
enddate_stats_sub <- as.POSIXct("2015-04-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

# What to run
runLdasin <- TRUE
runLdasout <- TRUE
runStats <- TRUE
writeStatsFile <- TRUE

###################################################################################################
# Run

load(inImg)


## ------------------------------------------------------------------------
# Setup processing functions

ProcessLdasin <- function(modLdasin, stdate=NULL, enddate=NULL) { 
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    modLdasin <- subset(modLdasin, modLdasin$POSIXct >= stdate & modLdasin$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    modLdasin <- subset(modLdasin, modLdasin$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    modLdasin <- subset(modLdasin, modLdasin$POSIXct <= enddate)
  }
  ## ------------------------------------------------------------------------
  # Calculate daily forcings
  # Adjust to dates to match SNOTEL daily report.
  # SNOTEL daily reports are derived from prior day's hourlies, and all 
  # times are PST (no daylight savings).
  # Adjust UTC to PST
  modLdasin$PST_time <- modLdasin$POSIXct - 8*3600
  # Calculate truncated date from PST time
  modLdasin$PST_date <- CalcDateTrunc(modLdasin$PST_time)
  # Shift by 1 day so aggregations match daily report
  modLdasin$PST_dateP1 <- modLdasin$PST_date + 1
  # Unit conversions
  modLdasin$RelHum <- 0.01 * with(modLdasin,                    
                            0.263*PSFC*Q2D*(exp((17.67*(T2D-273.16))/(T2D-29.65)))^(-1))
  modLdasin$Wind <- with(modLdasin, sqrt(U2D^2 + V2D^2))
  # Run daily aggs
  modLdasin.snoday <- plyr::ddply(modLdasin, plyr::.(statArg, PST_dateP1), plyr::summarize, 
                                 T2D_mean=mean(T2D), T2D_min=min(T2D), T2D_max=max(T2D),
                                 Q2D_mean=mean(Q2D), Q2D_min=min(Q2D), Q2D_max=max(Q2D),
                                 U2D_mean=mean(U2D), U2D_min=min(U2D), U2D_max=max(U2D),
                                 V2D_mean=mean(V2D), V2D_min=min(V2D), V2D_max=max(V2D),
                                 PSFC_mean=mean(PSFC), PSFC_min=min(PSFC), PSFC_max=max(PSFC),
                                 SWDOWN_mean=mean(SWDOWN), SWDOWN_min=min(SWDOWN), SWDOWN_max=max(SWDOWN),
                                 LWDOWN_mean=mean(LWDOWN), LWDOWN_min=min(LWDOWN), LWDOWN_max=max(LWDOWN),
                                 RelHum_mean=mean(RelHum), RelHum_min=min(RelHum), RelHum_max=max(RelHum),
                                 Wind_mean=mean(Wind), Wind_min=min(Wind), Wind_max=max(Wind),
                                 .parallel=TRUE)  
  # Add dummy POSIXct for ease of plotting
  modLdasin.snoday$POSIXct <- as.POSIXct(paste0(modLdasin.snoday$PST_dateP1, " 00:00"), tz="UTC")

  # Calculate truncated date from UTC time
  modLdasin$UTC_date <- CalcDateTrunc(modLdasin$POSIXct)
  # Run daily aggs
  modLdasin.utcday <- plyr::ddply(modLdasin, plyr::.(statArg, UTC_date), plyr::summarize,
                                 T2D_mean=mean(T2D), T2D_min=min(T2D), T2D_max=max(T2D),
                                 Q2D_mean=mean(Q2D), Q2D_min=min(Q2D), Q2D_max=max(Q2D),
				 U2D_mean=mean(U2D), U2D_min=min(U2D), U2D_max=max(U2D),
                                 V2D_mean=mean(V2D), V2D_min=min(V2D), V2D_max=max(V2D),
				 PSFC_mean=mean(PSFC), PSFC_min=min(PSFC), PSFC_max=max(PSFC),
                                 SWDOWN_mean=mean(SWDOWN), SWDOWN_min=min(SWDOWN), SWDOWN_max=max(SWDOWN),
				 LWDOWN_mean=mean(LWDOWN), LWDOWN_min=min(LWDOWN), LWDOWN_max=max(LWDOWN),
                                 RelHum_mean=mean(RelHum), RelHum_min=min(RelHum), RelHum_max=max(RelHum),
				 Wind_mean=mean(Wind), Wind_min=min(Wind), Wind_max=max(Wind),
                                 .parallel=TRUE)
  # Add dummy POSIXct for ease of plotting
  modLdasin.utcday$POSIXct <- as.POSIXct(paste0(modLdasin.utcday$UTC_date, " 00:00"), tz="UTC")

  return(list(modLdasin, modLdasin.snoday, modLdasin.utcday))
}


ProcessLdasout <- function(modLdasout, stdate=NULL, enddate=NULL) {
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    modLdasout <- subset(modLdasout, modLdasout$POSIXct >= stdate & modLdasout$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    modLdasout <- subset(modLdasout, modLdasout$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    modLdasout <- subset(modLdasout, modLdasout$POSIXct <= enddate)
  }
  modLdasout
}

CalcVarStats <- function(modDf=modLdasout, siteDf=sno.URG.sites, obsDf=sno.URG.data, 
			stdate=NULL, enddate=NULL, 
			flxCol.obs, flxCol.mod,
			idCol.obs="site_id", idCol.mod="statArg") {
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
#	siteDf[,idCol.mod] <- siteDf[,idCol.obs]
#	}
  results<-plyr::join(results, siteDf, by=idCol.obs)

  results[results=="Inf"]<-NA
  results[results=="-Inf"]<-NA
  return(results)
}


## ------------------------------------------------------------------------
# Run Processing

saveList <- c()
stats_sno_all <- data.frame()
stats_met_all <- data.frame()
for (j in 1:length(objSuffixList)) {
  # Get files
  objSuffix <- objSuffixList[j]
  # Generate daily forcing values
  if (runLdasin) {
	if (exists(paste0("modLdasin", objSuffix, "_SNO"))) {
		modLdasin <- get(paste0("modLdasin", objSuffix, "_SNO"))
  		modLdasinList <- ProcessLdasin(modLdasin, stdate=stdate, enddate=enddate)
  		assign(paste0("modLdasin", objSuffix, "_SNO"), modLdasinList[[1]])
  		assign(paste0("modLdasin", objSuffix, "_SNO.snod"), modLdasinList[[2]])
		assign(paste0("modLdasin", objSuffix, "_SNO.metd"), modLdasinList[[3]])
		saveList <- c(saveList, paste0("modLdasin", objSuffix, "_SNO"), 
					paste0("modLdasin", objSuffix, "_SNO.snod"),
					paste0("modLdasin", objSuffix, "_SNO.metd"))
		if (runStats) {
			# SNOTEL (daily)
			# Full run
         		results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="Tavg_K", flxCol.mod="T2D_mean")
          		results$run <- objSuffix
          		results$var <- "Tmean"
			results$seas <- "All"
          		stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

          		results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats, enddate=enddate_stats,
                        	        flxCol.obs="Tmin_K", flxCol.mod="T2D_min")
          		results$run <- objSuffix
          		results$var <- "Tmin"
			results$seas <- "All"
          		stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

          		results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats, enddate=enddate_stats,
                       		         flxCol.obs="Tmax_K", flxCol.mod="T2D_max")
          		results$run <- objSuffix
          		results$var <- "Tmax"
			results$seas <- "All"
          		stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

                        # Subset (e.g., spring)
                        results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tavg_K", flxCol.mod="T2D_mean")
                        results$run <- objSuffix
                        results$var <- "Tmean"
                        results$seas <- "Sub"
                        stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

                        results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tmin_K", flxCol.mod="T2D_min")
                        results$run <- objSuffix
                        results$var <- "Tmin"
                        results$seas <- "Sub"
                        stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

                        results <- CalcVarStats(modLdasinList[[2]], sno.URG.sites, sno.URG.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="Tmax_K", flxCol.mod="T2D_max")
                        results$run <- objSuffix
                        results$var <- "Tmax"
                        results$seas <- "Sub"
                        stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)


			# MET (hourly)
			# Full run
                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                        flxCol.obs="Tmean_K", flxCol.mod="T2D")
                        results$run <- objSuffix
                        results$var <- "Temp"
			results$seas <- "All"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="RH_mean", flxCol.mod="RelHum")
                        results$run <- objSuffix
                        results$var <- "RelHum"
			results$seas <- "All"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="SurfPressmean_Pa", flxCol.mod="PSFC")
                        results$run <- objSuffix
                        results$var <- "SurfPress"
			results$seas <- "All"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="Wind_mean", flxCol.mod="Wind")
                        results$run <- objSuffix
                        results$var <- "Wind"
			results$seas <- "All"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats, enddate=enddate_stats,
                                         flxCol.obs="SWRad_mean", flxCol.mod="SWDOWN")
                        results$run <- objSuffix
                        results$var <- "SWRad"
			results$seas <- "All"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        # Subset (e.g., spring)
                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                        flxCol.obs="Tmean_K", flxCol.mod="T2D")
                        results$run <- objSuffix
                        results$var <- "Temp"
                        results$seas <- "Sub"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="RH_mean", flxCol.mod="RelHum")
                        results$run <- objSuffix
                        results$var <- "RelHum"
                        results$seas <- "Sub"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="SurfPressmean_Pa", flxCol.mod="PSFC")
                        results$run <- objSuffix
                        results$var <- "SurfPress"
                        results$seas <- "Sub"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="Wind_mean", flxCol.mod="Wind")
                        results$run <- objSuffix
                        results$var <- "Wind"
                        results$seas <- "Sub"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                        results <- CalcVarStats(modLdasinList[[1]], met.URG.sites, met.URG.data.hr, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                         flxCol.obs="SWRad_mean", flxCol.mod="SWDOWN")
                        results$run <- objSuffix
                        results$var <- "SWRad"
                        results$seas <- "Sub"
                        stats_met_all <- plyr::rbind.fill(stats_met_all, results)

			}
		}
	}
  # Process LDASOUT
  if (runLdasout) {
	modLdasout <- get(paste0("modLdasout", objSuffix, "_SNO"))
        modLdasout <- ProcessLdasout(modLdasout, stdate=stdate, enddate=stopDates[j])
        assign(paste0("modLdasout", objSuffix, "_SNO"), modLdasout)
	saveList <- c(saveList, paste0("modLdasout", objSuffix, "_SNO"))
  	if (runStats) {
		# SNOTEL
		# Full time period
		results <- CalcVarStats(modLdasout, sno.URG.sites, sno.URG.data, stdate=stdate_stats, enddate=enddate_stats, 
				flxCol.obs="Prec_mm", flxCol.mod="DEL_ACCPRCP")
		results$run <- objSuffix
	  	results$var <- "Precip"
		results$seas <- "All"
	  	stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

	  	results <- CalcVarStats(modLdasout, sno.URG.sites, sno.URG.data, stdate=stdate_stats, enddate=enddate_stats, 
				flxCol.obs="SWE_mm", flxCol.mod="SNEQV")
	  	results$run <- objSuffix
	  	results$var <- "SWE"
		results$seas <- "All"
	  	stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)
		# Subset time period
                results <- CalcVarStats(modLdasout, sno.URG.sites, sno.URG.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub, 
                                flxCol.obs="Prec_mm", flxCol.mod="DEL_ACCPRCP")
                results$run <- objSuffix
                results$var <- "Precip"
                results$seas <- "Sub"
                stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)

                results <- CalcVarStats(modLdasout, sno.URG.sites, sno.URG.data, stdate=stdate_stats_sub, enddate=enddate_stats_sub, 
                                flxCol.obs="SWE_mm", flxCol.mod="SNEQV")
                results$run <- objSuffix
                results$var <- "SWE"
                results$seas <- "Sub"
                stats_sno_all <- plyr::rbind.fill(stats_sno_all, results)


		# MET
		# Full time period
                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="PrecTot", flxCol.mod="DEL_ACCPRCP")
                results$run <- objSuffix
                results$var <- "Precip"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="SnoDepmean_m", flxCol.mod="SNOWH")
                results$run <- objSuffix
                results$var <- "SnowDepth"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="SM1_mean", flxCol.mod="SOIL_M1")
                results$run <- objSuffix
                results$var <- "SoilM1"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="SM2_mean", flxCol.mod="SOIL_M2")
                results$run <- objSuffix
                results$var <- "SoilM2"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="SoilT1mean_K", flxCol.mod="SOIL_T1")
                results$run <- objSuffix
                results$var <- "SoilT1"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats, enddate=enddate_stats,
                                flxCol.obs="SoilT2mean_K", flxCol.mod="SOIL_T2")
                results$run <- objSuffix
                results$var <- "SoilT2"
		results$seas <- "All"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                # Subset time period
                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="PrecTot", flxCol.mod="DEL_ACCPRCP")
                results$run <- objSuffix
                results$var <- "Precip"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="SnoDepmean_m", flxCol.mod="SNOWH")
                results$run <- objSuffix
                results$var <- "SnowDepth"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="SM1_mean", flxCol.mod="SOIL_M1")
                results$run <- objSuffix
                results$var <- "SoilM1"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="SM2_mean", flxCol.mod="SOIL_M2")
                results$run <- objSuffix
                results$var <- "SoilM2"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="SoilT1mean_K", flxCol.mod="SOIL_T1")
                results$run <- objSuffix
                results$var <- "SoilT1"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)

                results <- CalcVarStats(modLdasout, met.URG.sites, met.URG.data.dy, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                                flxCol.obs="SoilT2mean_K", flxCol.mod="SOIL_T2")
                results$run <- objSuffix
                results$var <- "SoilT2"
                results$seas <- "Sub"
                stats_met_all <- plyr::rbind.fill(stats_met_all, results)


		}
  }
}

if (runStats) {
	saveList <- c(saveList, "stats_sno_all", "stats_met_all")
	if (writeStatsFile) {
      	  	# Change NAs to large negative for QGIS so data type is not affected
      		stats_tmp <- stats_sno_all
		stats_tmp[is.na(stats_tmp)]<-(-1e+30)
    		write.table(stats_tmp, file="stats_sno_all.txt", sep="\t", row.names=FALSE)
                stats_tmp <- stats_met_all
                stats_tmp[is.na(stats_tmp)]<-(-1e+30)
                write.table(stats_tmp, file="stats_met_all.txt", sep="\t", row.names=FALSE)
    		}
	}



## ------------------------------------------------------------------------
# Cleanup
save(list=saveList, file=outImg)

proc.time()
quit("no")
