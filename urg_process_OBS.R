#' ---
#' title: "ANALYSIS (OBS): Evaluate streamflow simulations over multiple basins with rwrfhydro"
#' author: "Aubrey Dugger"
#' ---
#' 
#' # Background
#' We are using WRF-Hydro to predict streamflow for multiple basins in the Upper Rio Grande for 
#' 2004-2014. We ran WRF-Hydro in LSM-only mode (no routing) with NoahMP as the LSM for the 10-year 
#' period with daily output. We want to evaluate model performance at various gage stations in the 
#' domain.

###################################################################################################
#' # Setup
#' 
#' Load the rwrfhydro package. 
## ------------------------------------------------------------------------
#library("rwrfhydro")
library(doMC)
registerDoMC(3)

#' Set the data paths for the Upper Rio Grande (test case general, model output, streamflow  
#' observations, basin masks).
## ------------------------------------------------------------------------
# Mask dataset so we can get station IDs
maskPath <- 'urg_masks_NEW.Rdata'
# Where to save the R workspace
rimgPath <- '../OBS/strflow_URG.Rdata'

###################################################################################################
#' # Download streamflow data from CO DWR website
#' 
#' Download and process data from the CO DWR website. First, we build a 
#' list of stations.
## ------------------------------------------------------------------------
load(maskPath)
stnList <- names(gage2basinList)

#' Then, we grab stage height and discharge data for the 2015 water year 
#' (to date) using GetCoDwrData.
## ------------------------------------------------------------------------
obsStr <- GetCoDwrData(siteIDs=stnList, 
                       paramCodes=c("GAGE_HT", "DISCHRG"), 
                       timeInt="raw", 
                       startDate="01/01/04", endDate="07/31/15")

#' We also bring in the manually created reservoir storage data.
## ------------------------------------------------------------------------
obsRes.plat <- GetCoDwrData("PLARESCO", 
                            paramCodes=c("STORAGE", "SURF_AC"), 
                            timeInt="raw", 
                            startDate="01/01/04", endDate="07/31/15")
#obsRes.plat <- read.table("../OBS/RESERVOIR/Platoro_storage_change_WY_2015.csv", 
#                          header=F, skip=1, sep=",", 
#                          stringsAsFactors=FALSE)
#names(obsRes.plat) <- c("Date.Time", "delstor_cms")
#obsRes.plat$POSIXct <- as.POSIXct(format(as.POSIXct(obsRes.plat$Date.Time, 
#                                                    format="%m/%d/%Y %H:%M", 
#                                                    tz="America/Denver"),
#                                  tz="UTC"), tz="UTC")
#obsRes.plat <- subset(obsRes.plat, !is.na(obsRes.plat$POSIXct))
obsRes.plat$Station <- "CONMOGCO"
# Remove erroneous values
obsRes.plat$stor_acft_out1 <- FillOutliers(obsRes.plat$STORAGE..AF., 1000)
# Smooth the reservoir volume time series to remove some of the noise.
# We use a 24-hr (48-hr when 30-min) smoother
obsRes.plat$stor_acft <- CalcRunningMean(obsRes.plat$stor_acft_out1, 96)
# Calculate timestep
obsRes.plat$deltime_secs <- 0
obsRes.plat$deltime_secs[2:nrow(obsRes.plat)] <- 
  as.integer(difftime(obsRes.plat$POSIXct[2:nrow(obsRes.plat)], 
                      obsRes.plat$POSIXct[1:(nrow(obsRes.plat)-1)], units="secs"))
# Recalculate storage in m3
obsRes.plat$stor_m3 <- obsRes.plat$stor_acft * 43560 * (0.3048^3)
# And then as a flowrate
obsRes.plat$delstor_m3 <- 0
obsRes.plat$delstor_m3[2:nrow(obsRes.plat)] <- diff(obsRes.plat$stor_m3)
obsRes.plat$delstor_cms <- obsRes.plat$delstor_m3/obsRes.plat$deltime_secs

#' Adjust streamflow rates by reservoir storage change rates. These are now
#' "naturalized" flows.
## ------------------------------------------------------------------------
obsStr$q_cms_adj <- obsStr$q_cms
obsStr <- plyr::join(obsStr, obsRes.plat[c("Station","POSIXct","delstor_cms")], 
                     by=c("Station", "POSIXct"), match="first")
obsStr$q_cms_adj[!(is.na(obsStr$delstor_cms))] <- obsStr$q_cms[!(is.na(obsStr$delstor_cms))] + 
  obsStr$delstor_cms[!(is.na(obsStr$delstor_cms))]


#' Until we can automate this from the data download side, we will have to manually set the gage 
#' drainage areas. We will set it up as an attribute to the observation dataframe.
## ------------------------------------------------------------------------
attr(obsStr, "area_sqmi") <- c(ALATERCO=107, CONMOGCO=282, RIOSFKCO=216, 
                               RIOWAGCO=780, SAGSAGCO=595, TRITURCO=45,
                               SANORTCO=110, RIODELCO=1320,
                               LOSORTCO=167, CONPLACO=40)
attr(obsStr, "gage_name") <- c(ALATERCO="ALAMOSA RIVER ABOVE TERRACE RESERVOIR",
                               CONMOGCO="CONEJOS RIVER NEAR MOGOTE",
                               RIOSFKCO="SOUTH FORK RIO GRANDE RIVER AT SOUTH FORK",
                               RIOWAGCO="RIO GRANDE RIVER AT WAGON WHEEL GAP",
                               SAGSAGCO="SAGUACHE CREEK NEAR SAGUACHE",
                               TRITURCO="TRINCHERA CREEK ABOVE TURNER'S RANCH",
                               SANORTCO="SAN ANTONIO RIVER AT ORTIZ",
                               RIODELCO="RIO GRANDE RIVER NEAR DEL NORTE",
                               LOSORTCO="LOS PINOS RIVER NEAR ORTIZ, CO",
                               CONPLACO="CONEJOS RIVER BELOW PLATORO RESERVOIR")

#' To access this attribute:
## ------------------------------------------------------------------------
#attributes(obsStr)$area_sqmi

#' We can double-check to make sure we have all of the gages we expect.
## ------------------------------------------------------------------------
#unique(obsStr$Station)

###################################################################################################
#' # Run aggregations - DAILY
#' 
#' We can aggregate the observed data to a daily timestep to match the model.
## ------------------------------------------------------------------------
obsStr$Date <- as.Date(trunc(as.POSIXct(format(obsStr$POSIXct, tz="UTC"), tz="UTC"), "days"))
obsStr.dy <- plyr::ddply(obsStr, plyr::.(Station, Date), 
                         plyr::summarise, mean_qcms=mean(q_cms, na.rm=TRUE), 
                         mean_qcms_adj=mean(q_cms_adj, na.rm=TRUE),
                         .parallel=TRUE)
# Unit conversion: m^3/s -> m^3/dy -> ft^3/dy -> ac-ft/dy
obsStr.dy$qvol_acft <- obsStr.dy$mean_qcms * 86400 / (0.3048^3) / 43560
obsStr.dy$qvol_acft_adj <- obsStr.dy$mean_qcms_adj * 86400 / (0.3048^3) / 43560

#' Let's add a POSIXct column for ease of calculations and plotting. We will make it plus
#' one date for better consistency with the model output.
## ------------------------------------------------------------------------
obsStr.dy$POSIXct <- as.POSIXct(paste0(obsStr.dy$Date+1," 00:00", 
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")

#' And we can calculate a cumulative volume for each water year
## ------------------------------------------------------------------------
obsStr.dy$wy <- CalcWaterYear(obsStr.dy$POSIXct)
wyList <- unique(obsStr.dy$wy)
gageList <- unique(obsStr.dy$Station)
obsStr.dy$cumqvol_acft <- 0
obsStr.dy$cumqvol_mm <- 0
obsStr.dy$cumqvol_acft_adj <- 0
obsStr.dy$cumqvol_mm_adj <- 0
for (i in 1:length(gageList)) {
  print(paste0("index ",i, " gage ", gageList[[i]]))
  tmpgage <- subset(obsStr.dy, obsStr.dy$Station==gageList[i])
  for (j in 1:length(wyList)) {
    tmpgagewy <- subset(tmpgage, tmpgage$wy==wyList[j])
    obsStr.dy$cumqvol_acft[obsStr.dy$Station==gageList[i] & obsStr.dy$wy==wyList[j]] <- CumsumNa(tmpgagewy$qvol_acft)
    obsStr.dy$cumqvol_mm[obsStr.dy$Station==gageList[i] & obsStr.dy$wy==wyList[j]] <- CumsumNa(tmpgagewy$qvol_acft) / 
      attr(obsStr, "area_sqmi")[[gageList[i]]] /
      (5280^2) * 43560 * 25.4 * 12
    obsStr.dy$cumqvol_acft_adj[obsStr.dy$Station==gageList[i] & obsStr.dy$wy==wyList[j]] <- CumsumNa(tmpgagewy$qvol_acft_adj)
    obsStr.dy$cumqvol_mm_adj[obsStr.dy$Station==gageList[i] & obsStr.dy$wy==wyList[j]] <- CumsumNa(tmpgagewy$qvol_acft_adj) / 
      attr(obsStr, "area_sqmi")[[gageList[i]]] /
      (5280^2) * 43560 * 25.4 * 12
    rm(tmpgagewy)
  }
  rm(tmpgage)
}

#' Let's also aggregate the reservoir data to a daily timestep so we can adjust the daily
#' outflow predictions. We will assign a date stamp of the following day to properly adjust
#' the modelled flows.
## ------------------------------------------------------------------------
obsRes.plat$Date <- as.Date(trunc(as.POSIXct(format(obsRes.plat$POSIXct, tz="UTC"), tz="UTC"), "days"))
obsRes.plat.dy <- plyr::ddply(obsRes.plat, plyr::.(Date), 
                         plyr::summarise, mean_delstorcms=mean(delstor_cms, na.rm=TRUE), 
                         .parallel=TRUE)
# Unit conversion: m^3/s -> m^3/dy -> ft^3/dy -> ac-ft/dy
obsRes.plat.dy$delvol_acft <- obsRes.plat.dy$mean_delstorcms * 86400 / (0.3048^3) / 43560

#' Let's add a POSIXct column for ease of calculations and plotting.
## ------------------------------------------------------------------------
obsRes.plat.dy$POSIXct <- as.POSIXct(paste0(obsRes.plat.dy$Date+1," 00:00", 
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")

########################################################################################################
#' # Run aggregations - MONTHLY
#' 
#' We can aggregate the observed data to a monthly timestep in acre-ft.
## ------------------------------------------------------------------------
# CO DWR gages
obsStr$mo <- as.integer(format(obsStr$POSIXct, "%m", tz="UTC"))
obsStr$yr <- as.integer(format(obsStr$POSIXct, "%Y", tz="UTC"))
obsStr.mo <- plyr::ddply(obsStr, plyr::.(Station, yr, mo), 
                         plyr::summarise, mean_qcms=mean(q_cms, na.rm=TRUE), 
                         mean_qcms_adj=mean(q_cms_adj, na.rm=TRUE), 
                         .parallel=TRUE)

# Unit conversion: m^3/s -> m^3/mo -> ft^3/mo -> ac-ft/mo
obsStr.mo$qvol_acft <- obsStr.mo$mean_qcms * 86400 *
  CalcMonthDays(obsStr.mo$mo, obsStr.mo$yr) /
  0.3048^3 / 43560
obsStr.mo$qvol_acft_adj <- obsStr.mo$mean_qcms_adj * 86400 *
  CalcMonthDays(obsStr.mo$mo, obsStr.mo$yr) /
  0.3048^3 / 43560

#' Let's add a POSIXct column for ease of plotting. We'll associate monthly vaues with the 1st 
#' day of each month for plotting.
## ------------------------------------------------------------------------
obsStr.mo$POSIXct <- as.POSIXct(paste0(obsStr.mo$yr,"-",obsStr.mo$mo,"-01", 
                                       format="%Y-%m-%d", tz="UTC"), tz="UTC")

########################################################################################################
#' Save results
save(obsStr, obsStr.dy, obsStr.mo, obsRes.plat, obsRes.plat.dy, file=rimgPath)
