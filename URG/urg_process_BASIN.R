###################################################################################################
# Setup
# 
# Load the rwrfhydro package. 
## ------------------------------------------------------------------------
library("rwrfhydro")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_masks_NEW.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG.Rdata")

# If you want to use R's multi-core capability (make sure  doMC is installed) specify the number 
# of cores.
## ------------------------------------------------------------------------
ncores <- 16
library(doMC)
registerDoMC(ncores)

# Where the raw data lives
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_snowmod_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS.Rdata'
#inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_NLDAS2dwnsc_SIMGM_BATSalb_fullrtng_BAS.Rdata'

inImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_ALL_RAW.Rdata'

# Where to save the processed data
#outImg <- 'urg_wy2015_NLDAS2dwnsc_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_snowmod_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS_PROCESSED.Rdata'
#outImg <- 'urg_wy2015_NLDAS2dwnsc_SIMGM_BATSalb_fullrtng_BAS_PROCESSED.Rdata'

outImg <- '/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_wy2015_BAS_PROCESSED.Rdata'

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
enddate <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

# Range dates for main stats
stdate_stats <- NULL
enddate_stats <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")

# Range dates for seasonal stats (e.g., spring)
stdate_stats_sub <- as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
enddate_stats_sub <- as.POSIXct("2015-08-27 00:00", format="%Y-%m-%d %H:%M", tz="UTC")



###################################################################################################
# Run

load(inImg)

# Setup lookups
stid2gage <- data.frame(st_id=names(stid2gageList), STAID=unlist(stid2gageList), stringsAsFactors=FALSE)
stid2gage$st_id <- as.integer(stid2gage$st_id)
#basin2gage <- read.table(basin2gageTbl, sep=",", header=TRUE)
#names(basin2gage)<-c("basin","STAID1","STAID2")

## ------------------------------------------------------------------------
# Setup processing functions

ProcessFrxstout <- function(modFrxstout, stid2gage, stdate=NULL, enddate=NULL) {
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    modFrxstout <- subset(modFrxstout, modFrxstout$POSIXct >= stdate & modFrxstout$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    modFrxstout <- subset(modFrxstout, modFrxstout$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    modFrxstout <- subset(modFrxstout, modFrxstout$POSIXct <= enddate)
  }
  # Bring in basin IDs
  modFrxstout <- plyr::join(modFrxstout, stid2gage, by="st_id")
  # Calculate accumulated flow
  modFrxstout$q_mm <- NA
  for (i in 1:nrow(modFrxstout)) {
    modFrxstout$q_mm[i] <- ifelse(is.na(modFrxstout$STAID[i]), NA, 
                                modFrxstout$q_cms[i]/
                                  (mskhyd.areaList[[modFrxstout$STAID[i]]]
                                   /100*1000*1000)*1000*(3600*24))
    }
  modFrxstout <- modFrxstout[order(modFrxstout$st_id, modFrxstout$POSIXct),]
  modFrxstout$ACCFLOW <- NA
  for (j in unique(modFrxstout$STAID)[!is.na(unique(modFrxstout$STAID))]) {
    tmp <- subset(modFrxstout, modFrxstout$STAID==j)
    qaccum <- cumsum(tmp$q_mm)
    modFrxstout$ACCFLOW[modFrxstout$STAID==j & !is.na(modFrxstout$STAID)] <- qaccum
  }
  modFrxstout
}

ProcessGwout <- function(modGwout, basin2gage, stdate=NULL, enddate=NULL) {
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    modGwout <- subset(modGwout, modGwout$POSIXct >= stdate & modGwout$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    modGwout <- subset(modGwout, modGwout$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    modGwout <- subset(modGwout, modGwout$POSIXct <= enddate)
  }
  # Bring in basin IDs
  #modGwout <- plyr::join(modGwout, basin2gage, by="basin", type="left", match="first")
  modGwout
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
  # Bring in basin IDs
  names(modLdasout)[which(names(modLdasout)=="STAID")]<-"STAID"
  modLdasout
}

CalcStrStats <- function(modDf, obsDf, stid2gageList, 
                      stdate=NULL, enddate=NULL, 
                      outfile=NULL) {
  sites <- names(stid2gageList)
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
      out$STAID <- stid2gageList[[sites[n]]]
      out
      }
    }
  #results<-plyr::join(results, stid2gageList, by="site_id")
  results[results=="Inf"]<-NA
  results[results=="-Inf"]<-NA
  if (!is.null(outfile)) {
    # Save outputs
    # Change NAs to large negative for QGIS so data type is not affected
    results2<-results
    results2[is.na(results2)]<-(-1e+30)
    write.table(results2, file=outfile, sep="\t", row.names=FALSE)
    }
  results
}


## ------------------------------------------------------------------------
# Run Processing
saveList <- c()
stats_str_all <- data.frame()
for (j in 1:length(objSuffixList)) {
  # Get files
  objSuffix <- objSuffixList[j]
  modLdasout <- get(paste0("modLdasout", objSuffix, "_BAS"))
  modGwout <- get(paste0("modGwout", objSuffix))
  modFrxstout <- get(paste0("modFrxstout", objSuffix))
  # Process
  modLdasout <- ProcessLdasout(modLdasout, stdate=stdate, enddate=stopDates[j])
  modGwout <- ProcessGwout(modGwout, basin2gage, stdate=stdate, enddate=stopDates[j])
  modFrxstout <- ProcessFrxstout(modFrxstout, stid2gage, stdate=stdate, enddate=stopDates[j])
  # Stats
  results <- CalcStrStats(modFrxstout, obsStr, stid2gageList, stdate=stdate_stats, enddate=enddate_stats, 
		outfile=paste0("stats_str", objSuffix, ".txt"))
  results$run <- objSuffix
  results$seas <- "Full"
  stats_str_all <- rbind(stats_str_all, results)
  results <- CalcStrStats(modFrxstout, obsStr, stid2gageList, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                outfile=paste0("stats_str_sub", objSuffix, ".txt"))
  results$run <- objSuffix
  results$seas <- "Sub"
  stats_str_all <- rbind(stats_str_all, results)
  # Save
  assign(paste0("modLdasout", objSuffix, "_BAS"), modLdasout)
  assign(paste0("modFrxstout", objSuffix), modFrxstout)
  assign(paste0("modGwout", objSuffix), modGwout)
  assign(paste0("stats_str", objSuffix), results)
  saveList <- c(saveList, paste0("modLdasout", objSuffix, "_BAS"), 
			paste0("modFrxstout", objSuffix), 
			paste0("modGwout", objSuffix),
			paste0("stats_str", objSuffix))
}

saveList <- c(saveList, "stats_str_all")

# Change NAs to large negative for QGIS so data type is not affected
stats_str_tmp <- stats_str_all
stats_str_tmp[is.na(stats_str_tmp)]<-(-1e+30)
write.table(stats_str_tmp, file="stats_str_all.txt", sep="\t", row.names=FALSE)


## ------------------------------------------------------------------------
# Cleanup
save(list=saveList, file=outImg)
#rm(modLdasout, modFrxstout, modGwout, sites, enddate, inImg, j, ncores, objSuffix, objSuffixList, outImg, stdate)
#save.image(outImg)

proc.time()
quit("no")
