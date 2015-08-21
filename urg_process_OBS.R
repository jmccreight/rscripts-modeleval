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
library("rwrfhydro")
library(doMC)
registerDoMC(16)

#' Set the data paths for the Upper Rio Grande (test case general, model output, streamflow  
#' observations, basin masks).
## ------------------------------------------------------------------------
# Mask dataset so we can get station IDs
maskPath <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_masks_NEW.Rdata'
# Path to existing streamflow data workspace
dataPath <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG.Rdata'
# Where to save the R workspace
rimgPath <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG2.Rdata'

# Download options
getStr <- FALSE
getRes <- TRUE


###################################################################################################
#' # Initialize Reservoir Info

# Reservoir to basin lookup
basin2resList <- list(	"CONPLACO"=c("PLARESCO"),
			"CONMOGCO"=c("PLARESCO"),
			"RIOWAGCO"=c("RIORESCO","CONRESCO"),
			"RIODELCO"=c("RIORESCO","CONRESCO"))

# Pan evap estimates
# adjustment factor for pan evap to reservoir
panfact <- 0.7
# reservoir areas in sq m
area_PLARESCO <- 600 * 4046.86 # in sq m (1 ac = 4046.86 m2)
area_RIORESCO <- 600 * 4046.86 # in sq m
area_CONRESCO <- 600 * 4046.86 # in sq m
# pan evap vals in average inches per month
obsPanEvap <- data.frame(mo=seq(1:12),
	days=c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
	PLATORO=c(1.41,  1.25,  2.81,  4.78,  5.86,  8.10,  6.57,  5.24,  5.52,  3.33,  1.35,  1.06),
	WAGONWHEEL=c(1.41,  1.25,  2.81,  4.78,  6.69,  7.90,  7.15,  5.81,  5.30,  2.61,  1.35,  1.06))
obsPanEvap$PLATORO_mmps <- panfact * obsPanEvap$PLATORO * 25.4 / obsPanEvap$days / 86400
obsPanEvap$WAGONWHEEL_mmps <- panfact * obsPanEvap$WAGONWHEEL * 25.4 / obsPanEvap$days / 86400
estResEvap <- data.frame(ResStation=c(rep("PLARESCO", 12), rep("RIORESCO", 12), rep("CONRESCO", 12)),
	area=c(rep(area_PLARESCO, 12), rep(area_RIORESCO, 12), rep(area_CONRESCO, 12)),
	mo=rep(seq(1:12),3),
	RESEVAP_mmps=c(obsPanEvap$PLATORO_mmps, obsPanEvap$WAGONWHEEL_mmps, obsPanEvap$WAGONWHEEL_mmps),
	stringsAsFactors=FALSE)
estResEvap$RESEVAP_cms <- estResEvap$RESEVAP_mmps * estResEvap$area / 1000


###################################################################################################
#' # Download streamflow data from CO DWR website
#' 
#' Download and process data from the CO DWR website. First, we build a 
#' list of stations.
## ------------------------------------------------------------------------
load(maskPath)
if (!is.null(dataPath)) load(dataPath)
stnList <- names(gage2basinList)

#' Then, we grab stage height and discharge data for the 2015 water year 
#' (to date) using GetCoDwrData.
## ------------------------------------------------------------------------
if (getStr) {
	obsStr <- GetCoDwrData(siteIDs=stnList, 
                       paramCodes=c("GAGE_HT", "DISCHRG"), 
                       timeInt="raw", 
                       startDate="01/01/04", endDate="07/31/15")
	}
#' We also download the reservoir storage data.
## ------------------------------------------------------------------------
if (getRes) {
	obsRes <- GetCoDwrData(siteIDs=c("PLARESCO", "RIORESCO", "CONRESCO"), 
                            paramCodes=c("STORAGE", "SURF_AC"), 
                            timeInt="raw", 
                            startDate="01/01/04", endDate="07/31/15")
}
names(obsRes)[names(obsRes)=="Station"] <- "ResStation"
# Process (noisy) reservoir data
obsResNew <- data.frame()
for (i in unique(obsRes$ResStation)) {
	tmp <- subset(obsRes, obsRes$ResStation==i)
	tmp <- tmp[order(tmp$POSIXct),]
	# Remove erroneous values
	tmp$stor_acft_out1 <- FillOutliers(tmp$STORAGE..AF., 1000)
	# Smooth the reservoir volume time series to remove some of the noise.
	# We use a 24-hr (48-hr when 30-min) smoother
	tmp$stor_acft <- CalcRunningMean(tmp$stor_acft_out1, 96)
	# Calculate timestep
	tmp$deltime_secs <- 0
	tmp$deltime_secs[2:nrow(tmp)] <- 
  		as.integer(difftime(tmp$POSIXct[2:nrow(tmp)], 
                      tmp$POSIXct[1:(nrow(tmp)-1)], units="secs"))
	# Recalculate storage in m3
	tmp$stor_m3 <- tmp$stor_acft * 43560 * (0.3048^3)
	# And then as a flowrate
	tmp$delstor_m3 <- 0
	tmp$delstor_m3[2:nrow(tmp)] <- diff(tmp$stor_m3)
	tmp$delstor_cms <- tmp$delstor_m3/tmp$deltime_secs
	# Res evap adjustment
	tmp$mo <- as.integer(format(tmp$POSIXct, format="%m"))
	tmp <- plyr::join(tmp, estResEvap[,c("ResStation","mo","RESEVAP_cms")], by=c("ResStation", "mo"))
	tmp$delstornet_cms <- tmp$delstor_cms + tmp$RESEVAP_cms
        obsResNew <- rbind(obsResNew, tmp)
        }
#obsResNew[obsResNew==Inf] <- NA
#obsResNew[obsResNew==(-Inf)] <- NA
#obsResNew[obsResNew==NaN] <- NA
obsRes <- obsResNew
rm(obsResNew)

#' Adjust streamflow rates by reservoir storage change rates. These are now
#' "naturalized" flows.
## ------------------------------------------------------------------------
obsStr$q_cms_adj <- obsStr$q_cms
obsStr$delstornet_cms <- NULL
for (i in names(basin2resList)) {
	tmp <- subset(obsStr, obsStr$Station==i)
	for (j in basin2resList[[i]]) {
		tmp <- plyr::join(tmp, subset(obsRes, obsRes$ResStation==j)[c("POSIXct", "delstornet_cms")], by="POSIXct", match="first")
		tmp$q_cms_adj[!(is.na(tmp$delstornet_cms))] <- tmp$q_cms_adj[!(is.na(tmp$delstornet_cms))] + 
  								tmp$delstornet_cms[!(is.na(tmp$delstornet_cms))]
		tmp$delstornet_cms <- NULL
		}
	obsStr <- rbind(subset(obsStr, obsStr$Station != i), tmp)
}

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
obsRes$Date <- as.Date(trunc(as.POSIXct(format(obsRes$POSIXct, tz="UTC"), tz="UTC"), "days"))
obsRes.dy <- plyr::ddply(obsRes, plyr::.(Date), 
                         plyr::summarise, mean_delstornetcms=mean(delstornet_cms, na.rm=TRUE), 
                         .parallel=TRUE)
# Unit conversion: m^3/s -> m^3/dy -> ft^3/dy -> ac-ft/dy
obsRes.dy$delvolnet_acft <- obsRes.dy$mean_delstornetcms * 86400 / (0.3048^3) / 43560

#' Let's add a POSIXct column for ease of calculations and plotting.
## ------------------------------------------------------------------------
obsRes.dy$POSIXct <- as.POSIXct(paste0(obsRes.dy$Date+1," 00:00", 
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
save(obsStr, obsStr.dy, obsStr.mo, obsRes, obsRes.dy, file=rimgPath)
