###################################################
##             MODEL OUTPUT READS                ##
###################################################


################## General Setup ##################

# Setup temp file
saveList <- c()
tmpRimg <- tempfile(fileext=".Rdata")
message(paste0("Temp output file:", tmpRimg))

# Setup lookups
stid2gage <- data.frame(st_id=names(stid2gageList), STAID=unlist(stid2gageList), stringsAsFactors=FALSE)
stid2gage$st_id <- as.integer(stid2gage$st_id)

# Get needed geo info
ncid <- ncdf4::nc_open(geoFile)
geoDX <- ncdf4::ncatt_get(ncid,varid=0,'DX')$value
ncdf4::nc_close(ncid)
hydDX <- geoDX/aggfact

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# INDEX PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readBasinLdasout | readBasinLdasin) {
 	# Basin means
 	basgeoIndex_Lev0 <- list()
 	basgeoIndex_Lev1 <- list()
 	basgeoIndex_Lev2 <- list()
 	basgeoIndex_Lev3 <- list()
 	basgeoIndex_Lev4 <- list()
 	for (i in 1:length(mskgeo.nameList)) {
   		basgeoIndex_Lev0[[as.character(mskgeo.nameList[[i]])]] <- list(start=c(mskgeo.minInds$x[i], mskgeo.minInds$y[i], 1),
                                               end=c(mskgeo.maxInds$x[i], mskgeo.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskgeo.List[[i]]))
   		basgeoIndex_Lev1[[as.character(mskgeo.nameList[[i]])]] <- list(start=c(mskgeo.minInds$x[i], 1, mskgeo.minInds$y[i], 1),
                                               end=c(mskgeo.maxInds$x[i], 1, mskgeo.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskgeo.List[[i]]))
   		basgeoIndex_Lev2[[as.character(mskgeo.nameList[[i]])]] <- list(start=c(mskgeo.minInds$x[i], 2, mskgeo.minInds$y[i], 1),
                                               end=c(mskgeo.maxInds$x[i], 2, mskgeo.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskgeo.List[[i]]))
   		basgeoIndex_Lev3[[as.character(mskgeo.nameList[[i]])]] <- list(start=c(mskgeo.minInds$x[i], 3, mskgeo.minInds$y[i], 1),
                                               end=c(mskgeo.maxInds$x[i], 3, mskgeo.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskgeo.List[[i]]))
   		basgeoIndex_Lev4[[as.character(mskgeo.nameList[[i]])]] <- list(start=c(mskgeo.minInds$x[i], 4, mskgeo.minInds$y[i], 1),
                                               end=c(mskgeo.maxInds$x[i], 4, mskgeo.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskgeo.List[[i]]))
 		}
 	}

if (readBasinRtout) {
         # Basin means
         bashydIndex_Lev0 <- list()
         bashydIndex_Lev1 <- list()
         bashydIndex_Lev2 <- list()
         bashydIndex_Lev3 <- list()
         bashydIndex_Lev4 <- list()
         for (i in 1:length(mskhyd.nameList)) {
                 bashydIndex_Lev0[[as.character(mskhyd.nameList[[i]])]] <- list(start=c(mskhyd.minInds$x[i], mskhyd.minInds$y[i], 1),
                                               end=c(mskhyd.maxInds$x[i], mskhyd.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskhyd.List[[i]]))
                 bashydIndex_Lev1[[as.character(mskhyd.nameList[[i]])]] <- list(start=c(mskhyd.minInds$x[i], 1, mskhyd.minInds$y[i], 1),
                                               end=c(mskhyd.maxInds$x[i], 1, mskhyd.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskhyd.List[[i]]))
                 bashydIndex_Lev2[[as.character(mskhyd.nameList[[i]])]] <- list(start=c(mskhyd.minInds$x[i], 2, mskhyd.minInds$y[i], 1),
                                               end=c(mskhyd.maxInds$x[i], 2, mskhyd.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskhyd.List[[i]]))
                 bashydIndex_Lev3[[as.character(mskhyd.nameList[[i]])]] <- list(start=c(mskhyd.minInds$x[i], 3, mskhyd.minInds$y[i], 1),
                                               end=c(mskhyd.maxInds$x[i], 3, mskhyd.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskhyd.List[[i]]))
                 bashydIndex_Lev4[[as.character(mskhyd.nameList[[i]])]] <- list(start=c(mskhyd.minInds$x[i], 4, mskhyd.minInds$y[i], 1),
                                               end=c(mskhyd.maxInds$x[i], 4, mskhyd.maxInds$y[i], 1), stat='basin_avg', arg=list(mskvar=mskhyd.List[[i]]))
                 }
         }

if (readAmfLdasout | readAmfLdasin) {
        # Ameriflux stations
        amfIndex_Lev0 <- list()
        amfIndex_Lev1 <- list()
        amfIndex_Lev2 <- list()
        amfIndex_Lev3 <- list()
        amfIndex_Lev4 <- list()
        # Ameriflux sites
        for (i in 1:length(ptgeo.amf$id)) {
                if (!is.na(ptgeo.amf$we[i]) & !is.na(ptgeo.amf$sn[i])) {
                        amfIndex_Lev0[[as.character(ptgeo.amf$id[i])]] <- list(start=c(ptgeo.amf$we[i], ptgeo.amf$sn[i], 1),
                                               end=c(ptgeo.amf$we[i], ptgeo.amf$sn[i], 1), stat="mean")
                        amfIndex_Lev1[[as.character(ptgeo.amf$id[i])]] <- list(start=c(ptgeo.amf$we[i], 1, ptgeo.amf$sn[i], 1),
                                          end=c(ptgeo.amf$we[i], 1, ptgeo.amf$sn[i], 1), stat="mean")
                        amfIndex_Lev2[[as.character(ptgeo.amf$id[i])]] <- list(start=c(ptgeo.amf$we[i], 2, ptgeo.amf$sn[i], 1),
                                               end=c(ptgeo.amf$we[i], 2, ptgeo.amf$sn[i], 1), stat="mean")
                        amfIndex_Lev3[[as.character(ptgeo.amf$id[i])]] <- list(start=c(ptgeo.amf$we[i], 3, ptgeo.amf$sn[i], 1),
                                               end=c(ptgeo.amf$we[i], 3, ptgeo.amf$sn[i], 1), stat="mean")
                        snoIndex_Lev4[[as.character(ptgeo.amf$id[i])]] <- list(start=c(ptgeo.amf$we[i], 4, ptgeo.amf$sn[i], 1),
                                               end=c(ptgeo.amf$we[i], 4, ptgeo.amf$sn[i], 1), stat="mean")
                        }
                }
        }

if (readSnoLdasout | readSnoLdasin) {
 	# SNOTEL points
 	snoIndex_Lev0 <- list()
 	snoIndex_Lev1 <- list()
 	snoIndex_Lev2 <- list()
 	snoIndex_Lev3 <- list()
 	snoIndex_Lev4 <- list()
 	for (i in 1:length(ptgeo.sno$id)) {
   		if (!is.na(ptgeo.sno$we[i]) & !is.na(ptgeo.sno$sn[i])) {
   			snoIndex_Lev0[[as.character(ptgeo.sno$id[i])]] <- list(start=c(ptgeo.sno$we[i], ptgeo.sno$sn[i], 1),
                                               end=c(ptgeo.sno$we[i], ptgeo.sno$sn[i], 1), stat="mean")
   			snoIndex_Lev1[[as.character(ptgeo.sno$id[i])]] <- list(start=c(ptgeo.sno$we[i], 1, ptgeo.sno$sn[i], 1),
                                          end=c(ptgeo.sno$we[i], 1, ptgeo.sno$sn[i], 1), stat="mean")
   			snoIndex_Lev2[[as.character(ptgeo.sno$id[i])]] <- list(start=c(ptgeo.sno$we[i], 2, ptgeo.sno$sn[i], 1),
                                               end=c(ptgeo.sno$we[i], 2, ptgeo.sno$sn[i], 1), stat="mean")
   			snoIndex_Lev3[[as.character(ptgeo.sno$id[i])]] <- list(start=c(ptgeo.sno$we[i], 3, ptgeo.sno$sn[i], 1),
                                               end=c(ptgeo.sno$we[i], 3, ptgeo.sno$sn[i], 1), stat="mean")
   			snoIndex_Lev4[[as.character(ptgeo.sno$id[i])]] <- list(start=c(ptgeo.sno$we[i], 4, ptgeo.sno$sn[i], 1),
                                               end=c(ptgeo.sno$we[i], 4, ptgeo.sno$sn[i], 1), stat="mean")
   			}
 		}
	}

if (readMetLdasout | readMetLdasin) {
	# Other met stations
        metIndex_Lev0 <- list()
        metIndex_Lev1 <- list()
        metIndex_Lev2 <- list()
        metIndex_Lev3 <- list()
        metIndex_Lev4 <- list()
        for (i in 1:length(ptgeo.met$id)) {
                if (!is.na(ptgeo.met$we[i]) & !is.na(ptgeo.met$sn[i])) {
                        metIndex_Lev0[[as.character(ptgeo.met$id[i])]] <- list(start=c(ptgeo.met$we[i], ptgeo.met$sn[i], 1),
                                               end=c(ptgeo.met$we[i], ptgeo.met$sn[i], 1), stat="mean")
                        metIndex_Lev1[[as.character(ptgeo.met$id[i])]] <- list(start=c(ptgeo.met$we[i], 1, ptgeo.met$sn[i], 1),
                                          end=c(ptgeo.met$we[i], 1, ptgeo.met$sn[i], 1), stat="mean")
                        metIndex_Lev2[[as.character(ptgeo.met$id[i])]] <- list(start=c(ptgeo.met$we[i], 2, ptgeo.met$sn[i], 1),
                                               end=c(ptgeo.met$we[i], 2, ptgeo.met$sn[i], 1), stat="mean")
                        metIndex_Lev3[[as.character(ptgeo.met$id[i])]] <- list(start=c(ptgeo.met$we[i], 3, ptgeo.met$sn[i], 1),
                                               end=c(ptgeo.met$we[i], 3, ptgeo.met$sn[i], 1), stat="mean")
                        metIndex_Lev4[[as.character(ptgeo.met$id[i])]] <- list(start=c(ptgeo.met$we[i], 4, ptgeo.met$sn[i], 1),
                                               end=c(ptgeo.met$we[i], 4, ptgeo.met$sn[i], 1), stat="mean")
                        }
                }
	}
	
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# RTOUT PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readBasinRtout) {
 
         ## ------------------------------------------------------------------------
         # Setup RTOUT files
 
         filesList <- list.files(path=modoutPath, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=TRUE)
         rtoutFilesList <- list( rtout = filesList)
 
         ## ------------------------------------------------------------------------
         # Setup variables
 
         varNames <- c('QSTRMVOLRT','SFCHEADSUBRT','QBDRYRT')
         rtoutVars <- as.list( varNames )
         names(rtoutVars) <- varNames
         rtoutVariableList <- list( rtout = rtoutVars )
 
         ## ------------------------------------------------------------------------
         # Setup indexes
 
         level0 <- bashydIndex_Lev0
         rtoutInd <- list( level0, level0, level0 )
         names(rtoutInd) <- names(rtoutVars)
         rtoutIndexList <- list( rtout = rtoutInd )
 
         ## ------------------------------------------------------------------------
         # Get data and flatten files

	 ## Loop through model run output directories
	 modRtout_BAS_tmp <- data.frame()
	 for (i in 1:length(modPathList)) {
        	modoutPath <- modPathList[i]
        	modoutTag <- modTagList[i]
         	# Setup RTOUT files
         	filesList <- list.files(path=modoutPath, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=TRUE)
		if (!is.null(readStart) | !is.null(readEnd)) {		
			filesList <- subDates(filesList, readStart, readEnd, rt2dt)
		}
         	rtoutFilesList <- list( rtout = filesList)
         	# Run basin means
		rtoutDF <- GetMultiNcdf(indexList=rtoutIndexList,
                	         variableList=rtoutVariableList,
                        	 filesList=rtoutFilesList, parallel=parallelFlag )
         	modRtout <- ReshapeMultiNcdf(rtoutDF)
         	modRtout <- modRtout[order(modRtout$statArg, modRtout$POSIXct),]
	 	modRtout$tag <- modoutTag
	 	modRtout_BAS_tmp <- plyr::rbind.fill(modRtout_BAS_tmp, modRtout)
         	rm(rtoutDF, modRtout, filesList, rtoutFilesList)
 	 	gc()
	 }
	 rm(varNames, rtoutVars, rtoutVariableList, level0, rtoutInd, rtoutIndexList)
	 saveList <- c(saveList, "modRtout_BAS_tmp")
 	 save(list=saveList, file=tmpRimg)

 } # end rtout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# LDASOUT PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readBasinLdasout | readAmfLdasout | readSnoLdasout | readMetLdasout) {
 
 	## ------------------------------------------------------------------------
 	# Setup variables

        if (varsLdasoutNFIE) {
                # SUBSET
                varNames <- c('SWFORC', 'LWFORC', 'ALBEDO', 'GRDFLX',
                        'LH', 'HFX', 'ETRAN', 'UGDRNOFF',
                        'SFCRNOFF', 'CANLIQ', 'CANICE', 'ACCPRCP',
			'ACCECAN', 'ACCEDIR', 'ACCETRAN', 'TRAD',
                         rep('SOIL_M',4),
                        'SNOWH', 'SNEQV')
                varLabels <- c('SWFORC', 'LWFORC', 'ALBEDO', 'GRDFLX',
                        'LH', 'HFX', 'ETRAN', 'UGDRNOFF',
                        'SFCRNOFF', 'CANLIQ', 'CANICE', 'ACCPRCP',
                        'ACCECAN', 'ACCEDIR', 'ACCETRAN', 'TRAD',
                         paste0('SOIL_M',1:4),
                        'SNOWH', 'SNEQV')
                ldasoutVars <- as.list( varNames )
                names(ldasoutVars) <- varLabels
                ldasoutVariableList <- list( ldasout = ldasoutVars )
                # INDEXES
                genIndex_LdasoutSub <- function(pref, ldasoutVars) {
                        level0 <- get(paste0(pref, "Index_Lev0"))
                        level1 <- get(paste0(pref, "Index_Lev1"))
                        level2 <- get(paste0(pref, "Index_Lev2"))
                        level3 <- get(paste0(pref, "Index_Lev3"))
                        level4 <- get(paste0(pref, "Index_Lev4"))
                        ldasoutInd <- list( level0, level0, level0, level0,
                                        level0, level0, level0, level0,
                                        level0, level0, level0, level0,
					level0, level0, level0, level0,
                                        level1, level2, level3, level4,
                                        level0, level0 )
                        names(ldasoutInd) <- names(ldasoutVars)
                        ldasoutIndexList <- list( ldasout = ldasoutInd )
                        ldasoutIndexList
                }
                # Run through reads
                if (readBasinLdasout) {
                        ldasoutBasIndexList <- genIndex_LdasoutSub("basgeo", ldasoutVars)
                        } # end runBasin
                if (readAmfLdasout) {
                        ldasoutAmfIndexList <- genIndex_LdasoutSub("amf", ldasoutVars)
                        } # end runAmf
                if (readSnoLdasout) {
                        ldasoutSnoIndexList <- genIndex_LdasoutSub("sno", ldasoutVars)
                        } # end runSno
                if (readMetLdasout) {
                        ldasoutMetIndexList <- genIndex_LdasoutSub("met", ldasoutVars)
                        } # end runMet

        } else { 
 
 	if (varsLdasoutSUB) {
 		# SUBSET
 		varNames <- c('ACCECAN', 'ACCEDIR', 'ETRAN', 'ACCPRCP', 
 			'CANICE', 'CANLIQ',
 			'SFCRNOFF','UGDRNOFF',
 			 rep('SOIL_M',4),
 			'GRDFLX', 'HFX', 'LH', 'LWFORC', 'SWFORC', 
 			'ALBEDO','TRAD',
 			'SNEQV', 'SNOWH')
 		varLabels <- c('ACCECAN', 'ACCEDIR', 'ETRAN', 'ACCPRCP',
                 	'CANICE', 'CANLIQ',
                 	'SFCRNOFF','UGDRNOFF',
                  	paste0('SOIL_M',1:4),
                 	'GRDFLX', 'HFX', 'LH', 'LWFORC', 'SWFORC', 
                 	'ALBEDO', 'TRAD',
                 	'SNEQV', 'SNOWH')
 		ldasoutVars <- as.list( varNames )
 		names(ldasoutVars) <- varLabels
 		ldasoutVariableList <- list( ldasout = ldasoutVars )
 		# INDEXES
		genIndex_LdasoutSub <- function(pref, ldasoutVars) {
			level0 <- get(paste0(pref, "Index_Lev0"))
			level1 <- get(paste0(pref, "Index_Lev1"))
			level2 <- get(paste0(pref, "Index_Lev2"))
			level3 <- get(paste0(pref, "Index_Lev3"))
			level4 <- get(paste0(pref, "Index_Lev4"))
			ldasoutInd <- list( level0, level0, level0, level0,
                                        level0, level0,
                                        level0, level0,
                                        level1, level2, level3, level4,
                                        level0, level0, level0, level0, level0,
                                        level0, level0, 
                                        level0, level0 )
                        names(ldasoutInd) <- names(ldasoutVars)
                        ldasoutIndexList <- list( ldasout = ldasoutInd )
			ldasoutIndexList
		}
		# Run through reads
 		if (readBasinLdasout) {
 			ldasoutBasIndexList <- genIndex_LdasoutSub("basgeo", ldasoutVars)
 			} # end runBasin
                if (readAmfLdasout) {
                        ldasoutAmfIndexList <- genIndex_LdasoutSub("amf", ldasoutVars)
                        } # end runAmf
                if (readSnoLdasout) {
                        ldasoutSnoIndexList <- genIndex_LdasoutSub("sno", ldasoutVars)
                        } # end runSno
 		if (readMetLdasout) {
                 	ldasoutMetIndexList <- genIndex_LdasoutSub("met", ldasoutVars)
 			} # end runMet
 
 	} else {
 		# ALL
 		varNames <- c('ACCECAN', 'ACCEDIR', 'ACCETRAN', 'ACCPRCP', 'ACSNOM', 'ACSNOW', 'ALBEDO', 'APAR', 'CANICE', 'CANLIQ',
                        'CH', 'CHB', 'CHB2', 'CHLEAF', 'CHUC', 'CHV', 'CHV2', 'CM', 'COSZ', 'EAH', 
                        'ECAN', 'EDIR', 'EMISS', 'ETRAN', 'EVB', 'EVC', 'EVG', 'FASTCP', 'FIRA', 'FSA', 
                        'FSNO', 'FVEG', 'FWET', 'GHB', 'GHV', 'GPP', 'GRDFLX', 'HFX', 'IRB', 'IRC', 
                        'IRG', 'ISLTYP', 'ISNOW', 'IVGTYP', 'LAI', 'LFMASS', 'LH', 'LWFORC', 'NEE', 'NPP', 
                        'PSN', 'Q2MB', 'Q2MV', 'QSNOW', 'RAINRATE', 'RTMASS', 'SAG', 'SAI', 'SAV', 'SFCRNOFF', 
                        'SHB', 'SHC', 'SHG', 'SNEQV',
                         rep('SNICE',3), 
                         rep('SNLIQ',3), 
                         'SNOWH', 
                         rep('SNOW_T',3),
                         rep('SOIL_M',4), 
                         rep('SOIL_T',4), 
                         rep('SOIL_W',4), 
                         'STBLCP', 'STMASS', 'SWFORC', 'T2MB', 'T2MV', 'TAH', 'TG', 'TGB', 'TGV', 'TR',
                         'TRAD', 'TV', 'UGDRNOFF', 'WA', 'WOOD', 'WT', 
                         rep('ZSNSO_SN',3), 
                         'ZWT')
 		varLabels <- c('ACCECAN', 'ACCEDIR', 'ACCETRAN', 'ACCPRCP', 'ACSNOM', 'ACSNOW', 'ALBEDO', 'APAR', 'CANICE', 'CANLIQ',
                         'CH', 'CHB', 'CHB2', 'CHLEAF', 'CHUC', 'CHV', 'CHV2', 'CM', 'COSZ', 'EAH', 
                         'ECAN', 'EDIR', 'EMISS', 'ETRAN', 'EVB', 'EVC', 'EVG', 'FASTCP', 'FIRA', 'FSA', 
                         'FSNO', 'FVEG', 'FWET', 'GHB', 'GHV', 'GPP', 'GRDFLX', 'HFX', 'IRB', 'IRC', 
                         'IRG', 'ISLTYP', 'ISNOW', 'IVGTYP', 'LAI', 'LFMASS', 'LH', 'LWFORC', 'NEE', 'NPP', 
                         'PSN', 'Q2MB', 'Q2MV', 'QSNOW', 'RAINRATE', 'RTMASS', 'SAG', 'SAI', 'SAV', 'SFCRNOFF', 
                         'SHB', 'SHC', 'SHG', 'SNEQV',
                         paste0('SNICE',1:3), 
                         paste0('SNLIQ',1:3), 
                         'SNOWH', 
                         paste0('SNOW_T',1:3), 
                         paste0('SOIL_M',1:4), 
                         paste0('SOIL_T',1:4), 
                         paste0('SOIL_W',1:4),
                         'STBLCP', 'STMASS', 'SWFORC', 'T2MB', 'T2MV', 'TAH', 'TG', 'TGB', 'TGV', 'TR',
                         'TRAD', 'TV', 'UGDRNOFF', 'WA', 'WOOD', 'WT',
                         paste0('ZSNSO_SN',1:3),
                         'ZWT')
         	ldasoutVars <- as.list( varNames )
         	names(ldasoutVars) <- varLabels
         	ldasoutVariableList <- list( ldasout = ldasoutVars )
 		# INDEXES
	        genIndex_LdasoutAll <- function(pref, ldasoutVars) {
                        level0 <- get(paste0(pref, "Index_Lev0"))
			level0 <- get(paste0(pref, "Index_Lev0"))
                        level0 <- get(paste0(pref, "Index_Lev0"))
                        level0 <- get(paste0(pref, "Index_Lev0"))
                        level0 <- get(paste0(pref, "Index_Lev0"))
                        ldasoutInd <- list( level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0,
                                     level1, level2, level3,
                                     level1, level2, level3,
                                     level0,
                                     level1, level2, level3,
                                     level1, level2, level3, level4,
                                     level1, level2, level3, level4,
                                     level1, level2, level3, level4,
                                     level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
                                     level0, level0, level0, level0, level0, level0,
                                     level1, level2, level3,
                                     level0 )
                        names(ldasoutInd) <- names(ldasoutVars)
                        ldasoutIndexList <- list( ldasout = ldasoutInd )
			ldasoutIndexList
		}
		# Run through reads
         	if (readBasinLdasout) {
		        ldasoutBasIndexList <- genIndex_LdasoutAll("basgeo", ldasoutVars)
 			} # end runBasin
         	if (readAmfLdasout) {
                 	ldasoutAmfIndexList <- genIndex_LdasoutAll("amf", ldasoutVars)
                 	} # end runAmf
                if (readSnoLdasout) {
                        ldasoutSnoIndexList <- genIndex_LdasoutAll("sno", ldasoutVars)
                        } # end runSno
                if (readMetLdasout) {
                        ldasoutMetIndexList <- genIndex_LdasoutAll("met", ldasoutVars)
                        } # end runMet

 	} # end ifelse vars subset
	} # end ifelse NFIE
 
 	## ------------------------------------------------------------------------
 	# Get data and flatten files
 	## ------------------------------------------------------------------------
	# GENERAL FUNCTION
	genReads_Ldasout <- function(modPathList=modPathList, modTagList=modTagList,
					ldasoutIndexList,
					ldasoutVariableList=ldasoutVariableList, 
					ldasoutFilesList=ldasoutFilesList, 
					parallelFlag=parallelFlag,
					readStart=readStart, readEnd=readEnd, ldas2dt=ldas2dt) {
                modLdasout_tmp <- data.frame()
                for (i in 1:length(modPathList)) {
                        modoutPath <- modPathList[i]
                        modoutTag <- modTagList[i]
                        # Setup LDASOUT files
                        filesList <- list.files(path=modoutPath, pattern=glob2rx('*.LDASOUT_DOMAIN*'), full.names=TRUE)
			if (!is.null(readStart) | !is.null(readEnd)) {
                        	filesList <- subDates(filesList, readStart, readEnd, ldas2dt)
                	}
                        ldasoutFilesList <- list( ldasout = filesList)
                        # Run basin means
                        ldasoutDF <- GetMultiNcdf(indexList=ldasoutIndexList,
                                           variableList=ldasoutVariableList,
                                           filesList=ldasoutFilesList,
                                           parallel=parallelFlag )
                        modLdasout <- ReshapeMultiNcdf(ldasoutDF)
                        modLdasout <- CalcNoahmpFluxes(modLdasout, "statArg")
                        modLdasout$tag <- modoutTag
                        modLdasout_tmp <- plyr::rbind.fill(modLdasout_tmp, modLdasout)
                }
		modLdasout_tmp
	}
	# BASIN
 	if (readBasinLdasout) {
		modLdasout_BAS_tmp <- genReads_Ldasout(ldasoutIndexList=ldasoutBasIndexList)
		saveList <- c(saveList, "modLdasout_BAS_tmp")
         	save(list=saveList, file=tmpRimg)
 	}	
 	# AMERIFLUX
 	if (readAmfLdasout) {
                modLdasout_AMF_tmp <- genReads_Ldasout(ldasoutIndexList=ldasoutAmfIndexList)
                saveList <- c(saveList, "modLdasout_AMF_tmp")
                save(list=saveList, file=tmpRimg)
	}
	# SNOTEL
        if (readSnoLdasout) {
                modLdasout_SNO_tmp <- genReads_Ldasout(ldasoutIndexList=ldasoutSnoIndexList)
                saveList <- c(saveList, "modLdasout_SNO_tmp")
                save(list=saveList, file=tmpRimg)
        }
	# MET
        if (readMetLdasout) {
                modLdasout_MET_tmp <- genReads_Ldasout(ldasoutIndexList=ldasoutMetIndexList)
                saveList <- c(saveList, "modLdasout_MET_tmp")
                save(list=saveList, file=tmpRimg)
        }


 } # end ldasout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# GW PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readGwout) {

	## Loop through model run output directories
	modGwout_tmp <- data.frame()
	for (i in 1:length(modPathList)) {
        	modoutPath <- modPathList[i]
        	modoutTag <- modTagList[i]
		# Read GW
 		modGwout <- ReadGwOut(paste0(modoutPath, '/GW_outflow.txt'))
		# Filter out non-unique dates. Take values from latest run if dups.
 		modGwout <- modGwout[nrow(modGwout):1,]
 		modGwout$uni <- paste(modGwout$basin, modGwout$timest, sep=",")
 		modGwout <- modGwout[!(duplicated(modGwout$uni)),]
 		modGwout$uni <- NULL
 		modGwout <- modGwout[nrow(modGwout):1,]
		modGwout$tag <- modoutTag
                modGwout_tmp <- plyr::rbind.fill(modGwout_tmp, modGwout)
		rm(modGwout)
		gc()
	}
        saveList <- c(saveList, "modGwout_tmp")
        save(list=saveList, file=tmpRimg)

} # end gwout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# FRXST PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readFrxstout) {

        ## Loop through model run output directories
        modFrxstout_tmp <- data.frame()
        for (i in 1:length(modPathList)) {
                modoutPath <- modPathList[i]
                modoutTag <- modTagList[i]
                # Read STR
        	modFrxstout <- ReadFrxstPts(paste0(modoutPath, '/frxst_pts_out.txt'))
        	# Filter out non-unique dates. Take values from latest run if dups.
        	modFrxstout <- modFrxstout[nrow(modFrxstout):1,]
        	modFrxstout$uni <- paste(modFrxstout$st_id, modFrxstout$timest, sep=",")
        	modFrxstout <- modFrxstout[!(duplicated(modFrxstout$uni)),]
        	modFrxstout$uni <- NULL
        	modFrxstout <- modFrxstout[nrow(modFrxstout):1,]
                modFrxstout <- modFrxstout[nrow(modFrxstout):1,]
		# Bring in basin IDs
  		modFrxstout <- plyr::join(modFrxstout, stid2gage, by="st_id")
  		# Calculate accumulated flow
  		modFrxstout$q_mm <- NA
  		for (j in 1:nrow(modFrxstout)) {
			ts <- modFrxstout$secs[j] - ifelse(j==1, 0, modFrxstout$secs[j-1])
    			modFrxstout$q_mm[j] <- ifelse(is.na(modFrxstout$STAID[j]), NA,
                                	modFrxstout$q_cms[j]/
                                  	(mskhyd.areaList[[modFrxstout$STAID[j]]]
                                   	*hydDX*hydDX)*1000*ts)
    		}
  		modFrxstout <- modFrxstout[order(modFrxstout$st_id, modFrxstout$POSIXct),]
  		modFrxstout$ACCFLOW <- NA
  		for (j in unique(modFrxstout$STAID)[!is.na(unique(modFrxstout$STAID))]) {
    			tmp <- subset(modFrxstout, modFrxstout$STAID==j)
    			qaccum <- cumsum(tmp$q_mm)
    			modFrxstout$ACCFLOW[modFrxstout$STAID==j & !is.na(modFrxstout$STAID)] <- qaccum
  		}
		# Add model run tag and bind
                modFrxstout$tag <- modoutTag
                modFrxstout_tmp <- plyr::rbind.fill(modFrxstout_tmp, modFrxstout)
                rm(modFrxstout)
                gc()
        }
        saveList <- c(saveList, "modFrxstout_tmp")
        save(list=saveList, file=tmpRimg)

} # end frxstout processing



## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# LDASIN PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (readBasinLdasin | readAmfLdasin | readSnoLdasin | readMetLdasin) {

        ## ------------------------------------------------------------------------
        # Setup variables

        varNames <- c('T2D', 'Q2D', 'U2D', 'V2D', 'PSFC', 'SWDOWN', 'LWDOWN','RAINRATE')
        ldasinVars <- as.list( varNames )
        names(ldasinVars) <- varNames
        ldasinVariableList <- list( ldasin = ldasinVars )
       
        # INDEXES
	# GENERAL FUNCTION
        genIndex_Ldasin <- function(pref, ldasinVars=ldasinVars) {
        	level0 <- get(paste0(pref, "Index_Lev0"))
		ldasinInd <- list( level0, level0, level0, level0, level0, level0, level0, level0 )
                names(ldasinInd) <- names(ldasinVars)
                ldasinIndexList <- list( ldasin = ldasinInd )        
                ldasinIndexList
                }
        # Run through reads
        if (readBasinLdasin) {
                ldasinBasIndexList <- genIndex_Ldasin("basgeo")
                } # end runBasin
        if (readAmfLdasin) {
                ldasinAmfIndexList <- genIndex_Ldasin("amf")
                } # end runAmf
        if (readSnoLdasin) {
                ldasinSnoIndexList <- genIndex_Ldasin("sno")
                } # end runSno
        if (readMetLdasin) {
                ldasinMetIndexList <- genIndex_Ldasin("met")
                } # end runMet

        ## ------------------------------------------------------------------------
        # Get data and flatten files
        ## ------------------------------------------------------------------------
        # GENERAL FUNCTION
        genReads_Ldasin <- function(forcPathList=forcPathList, forcTagList=forcTagList,
                                        ldasinIndexList,
                                        ldasinVariableList=ldasinVariableList,
                                        ldasinFilesList=ldasinFilesList,
                                        parallelFlag=parallelFlag,
					readStart=readStart, readEnd=readEnd, ldas2dt=ldas2dt) {
                modLdasin_tmp <- data.frame()
                for (i in 1:length(forcPathList)) {
                        forcPath <- forcPathList[i]
                        forcTag <- forcTagList[i]
                        # Setup LDASIN files
                        filesList <- list.files(path=forcPath, pattern=glob2rx('*.LDASIN_DOMAIN*'), full.names=TRUE)
                        if (!is.null(readStart) | !is.null(readEnd)) {
                                filesList <- subDates(filesList, readStart, readEnd, ldas2dt)
                        }                
        		ldasinFilesList <- list( ldasin = filesList)
                        # Run basin means
                        ldasinDF <- GetMultiNcdf(indexList=ldasinIndexList,
                                           variableList=ldasinVariableList,
                                           filesList=ldasinFilesList,
                                           parallel=parallelFlag )
                        modLdasin <- ReshapeMultiNcdf(ldasinDF)
			modLdasin <- modLdasin[order(modLdasin$statArg, modLdasin$POSIXct),]
                        modLdasin$tag <- forcTag
                        modLdasin_tmp <- plyr::rbind.fill(modLdasin_tmp, modLdasin)
                }
                modLdasin_tmp
        }
        # BASIN
        if (readBasinLdasin) {
                modLdasin_BAS_tmp <- genReads_Ldasin(ldasinIndexList=ldasinBasIndexList)
                saveList <- c(saveList, "modLdasin_BAS_tmp")
                save(list=saveList, file=tmpRimg)
        }
        # AMERIFLUX
        if (readAmfLdasin) {
                modLdasin_AMF_tmp <- genReads_Ldasin(ldasinIndexList=ldasinAmfIndexList)
                saveList <- c(saveList, "modLdasin_AMF_tmp")
                save(list=saveList, file=tmpRimg)
        }
        # SNOTEL
        if (readSnoLdasin) {
                modLdasin_SNO_tmp <- genReads_Ldasin(ldasinIndexList=ldasinSnoIndexList)
                saveList <- c(saveList, "modLdasin_SNO_tmp")
                save(list=saveList, file=tmpRimg)
        }
        # MET
        if (readMetLdasin) {
                modLdasin_MET_tmp <- genReads_Ldasin(ldasinIndexList=ldasinMetIndexList)
                saveList <- c(saveList, "modLdasin_MET_tmp")
                save(list=saveList, file=tmpRimg)
        }

 } # end LDASIN

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# SAVE DATA
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

saveList <- c()

if (append2exist) {
	load(modReadFileIn)
	if (readBasinLdasout) {
		if (exists("modLdasout_BAS")) {
			modLdasout_BAS <- plyr::rbind.fill(modLdasout_BAS, modLdasout_BAS_tmp)
		} else {
			modLdasout_BAS <- modLdasout_BAS_tmp
		}
		saveList <- c(saveList, "modLdasout_BAS")
	}
        if (readAmfdasout) {
                if (exists("modLdasout_AMF")) {
                        modLdasout_AMF <- plyr::rbind.fill(modLdasout_AMF, modLdasout_AMF_tmp)
                } else {
                        modLdasout_AMF <- modLdasout_AMF_tmp
                }
		saveList <- c(saveList, "modLdasout_AMF")
        }
        if (readSnoLdasout) {
                if (exists("modLdasout_SNO")) {
                        modLdasout_SNO <- plyr::rbind.fill(modLdasout_SNO, modLdasout_SNO_tmp)
                } else {
                        modLdasout_SNO <- modLdasout_SNO_tmp
                }
		saveList <- c(saveList, "modLdasout_SNO")
        }
        if (readMetLdasout) {
                if (exists("modLdasout_MET")) {
                        modLdasout_MET <- plyr::rbind.fill(modLdasout_MET, modLdasout_MET_tmp)
                } else {
                        modLdasout_MET <- modLdasout_MET_tmp
                }
		saveList <- c(saveList, "modLdasout_MET")
        }
        if (readBasinLdasin) {
                if (exists("modLdasin_BAS")) {
                        modLdasin_BAS <- plyr::rbind.fill(modLdasin_BAS, modLdasin_BAS_tmp)
                } else {
                        modLdasin_BAS <- modLdasin_BAS_tmp
                }
		saveList <- c(saveList, "modLdasin_BAS")
        }
        if (readAmfdasin) {
                if (exists("modLdasin_AMF")) {
                        modLdasin_AMF <- plyr::rbind.fill(modLdasin_AMF, modLdasin_AMF_tmp)
                } else {
                        modLdasin_AMF <- modLdasin_AMF_tmp
                }
		saveList <- c(saveList, "modLdasin_AMF")
        }
        if (readSnoLdasin) {
                if (exists("modLdasin_SNO")) {
                        modLdasin_SNO <- plyr::rbind.fill(modLdasin_SNO, modLdasin_SNO_tmp)
                } else {
                        modLdasin_SNO <- modLdasin_SNO_tmp
                }
		saveList <- c(saveList, "modLdasin_SNO")
        }
        if (readMetLdasin) {
                if (exists("modLdasin_MET")) {
                        modLdasin_MET <- plyr::rbind.fill(modLdasin_MET, modLdasin_MET_tmp)
                } else {
                        modLdasin_MET <- modLdasin_MET_tmp
                }
		saveList <- c(saveList, "modLdasin_MET")
        }
	if (readBasinRtout) {
		if (exists("modRtout_BAS")) {
                        modRtout_BAS <- plyr::rbind.fill(modRtout_BAS, modRtout_BAS_tmp)
                } else {
                        modRtout_BAS <- modRtout_BAS_tmp
                }
		saveList <- c(saveList, "modRtout_BAS")
	}
        if (readGwout) {
                if (exists("modGwout")) {
                        modGwout <- plyr::rbind.fill(modGwout, modGwout_tmp)
                } else {
                        modGwout <- modGwout_tmp
                }
		saveList <- c(saveList, "modGwout")
        }
        if (readFrxstout) {
                if (exists("modFrxstout")) {
                        modFrxstout <- plyr::rbind.fill(modFrxstout, modFrxstout_tmp)
                } else {
                        modFrxstout <- modFrxstout_tmp
                }
		saveList <- c(saveList, "modFrxstout")
        }

} else {
        if (readBasinLdasout) {
                        modLdasout_BAS <- modLdasout_BAS_tmp
			saveList <- c(saveList, "modLdasout_BAS")
        }
        if (readAmfdasout) {
                        modLdasout_AMF <- modLdasout_AMF_tmp
			saveList <- c(saveList, "modLdasout_AMF")
        }
        if (readSnoLdasout) {
                        modLdasout_SNO <- modLdasout_SNO_tmp
			saveList <- c(saveList, "modLdasout_SNO")
        }
        if (readMetLdasout) {
                        modLdasout_MET <- modLdasout_MET_tmp
			saveList <- c(saveList, "modLdasout_MET")
        }
        if (readBasinLdasin) {
                        modLdasin_BAS <- modLdasin_BAS_tmp
			saveList <- c(saveList, "modLdasin_BAS")
        }
        if (readAmfdasin) {
                        modLdasin_AMF <- modLdasin_AMF_tmp
			saveList <- c(saveList, "modLdasin_AMF")
        }
        if (readSnoLdasin) {
                        modLdasin_SNO <- modLdasin_SNO_tmp
			saveList <- c(saveList, "modLdasin_SNO")
        }
        if (readMetLdasin) {
                        modLdasin_MET <- modLdasin_MET_tmp
			saveList <- c(saveList, "modLdasin_MET")
        }
        if (readBasinRtout) {
                        modRtout_BAS <- modRtout_BAS_tmp
			saveList <- c(saveList, "modRtout_BAS")
        }
        if (readGwout) {
                        modGwout <- modGwout_tmp
			saveList <- c(saveList, "modGwout")
        }
        if (readFrxstout) {
                        modFrxstout <- modFrxstout_tmp
			saveList <- c(saveList, "modFrxstout")
        }
}

save(saveList, file=modReadFileOut)
#file.remove(tmpRimg)

quit("no")
proc.time()

