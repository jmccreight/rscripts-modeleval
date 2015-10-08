###################################################################################################
# Run

load(modOutImg)

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
modLdasout_BAS <- data.frame()
modGwout_BAS <- data.frame()
modFrxstout_BAS <- data.frame()

for (j in 1:length(objTagList)) {
  # Get files
  objSuffix <- objTagList[j]
  modLdasout_tmp <- get(paste0("modLdasout_", objSuffix, "_BAS"))
  modGwout_tmp <- get(paste0("modGwout_", objSuffix))
  modFrxstout_tmp <- get(paste0("modFrxstout_", objSuffix))
  # Process
  modLdasout_tmp <- ProcessLdasout(modLdasout_tmp, stdate=stdate, enddate=stopDates[j])
  modGwout_tmp <- ProcessGwout(modGwout_tmp, basin2gage, stdate=stdate, enddate=stopDates[j])
  modFrxstout_tmp <- ProcessFrxstout(modFrxstout_tmp, stid2gage, stdate=stdate, enddate=stopDates[j])
  # Stats
  results <- CalcStrStats(modFrxstout_tmp, obsStr, stid2gageList, stdate=stdate_stats, enddate=enddate_stats, 
		outfile=paste0("stats_str", objSuffix, ".txt"))
  results$run <- objSuffix
  results$seas <- "Full"
  stats_str_all <- rbind(stats_str_all, results)
  results <- CalcStrStats(modFrxstout_tmp, obsStr, stid2gageList, stdate=stdate_stats_sub, enddate=enddate_stats_sub,
                outfile=paste0("stats_str_sub", objSuffix, ".txt"))
  results$run <- objSuffix
  results$seas <- "Sub"
  stats_str_all <- rbind(stats_str_all, results)
  
  ##add current tag to each data frame
  modLdasout_tmp$modTag <- objSuffix
  modGwout_tmp$modTag <- objSuffix
  modFrxstout_tmp$modTag <- objSuffix

  modLdasout_BAS <- rbind(modLdasout_BAS,modLdasout_tmp)
  modGwout_BAS <- rbind(modGwout_BAS,modGwout_tmp)
  modFrxstout_BAS <- rbind(modFrxstout_BAS,modFrxstout_tmp)

}

##saveList <- c(saveList, "stats_str_all")
saveList <- c("modLdasout_BAS","modGwout_BAS","modFrxstout_BAS","stats_str_all")

# Change NAs to large negative for QGIS so data type is not affected
stats_str_tmp <- stats_str_all
stats_str_tmp[is.na(stats_str_tmp)]<-(-1e+30)
write.table(stats_str_tmp, file="stats_str_all.txt", sep="\t", row.names=FALSE)


## ------------------------------------------------------------------------
# Cleanup`
save(list=saveList, file=basinOutImg)
#rm(modLdasout, modFrxstout, modGwout, sites, enddate, modOutImg, j, ncores, objSuffix, objTagList, basinOutImg, stdate)
#save.image(basinOutImg)

##proc.time()
##quit("no")
