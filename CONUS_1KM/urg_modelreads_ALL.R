###################################################################################################
# Run

for(j in 1:length(objTagList)){

# Setup basin mean function
basin_avg <- function(myvar, mskvar, minValid=-1e+29) {
   myvar[which(myvar<minValid)]<-NA
   sum(mskvar*myvar, na.rm=TRUE)/sum(mskvar, na.rm=TRUE)
 }

saveList <- c()

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# INDEX PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runBasinLdasout) {
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

if (runBasinRtout) {
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


if (runSnoLdasout | runSnoLdasin) {
 	# SNOTEL points
 	snoIndex_Lev0 <- list()
 	snoIndex_Lev1 <- list()
 	snoIndex_Lev2 <- list()
 	snoIndex_Lev3 <- list()
 	snoIndex_Lev4 <- list()
 	for (i in 1:length(sno.URG.sites$site_id)) {
   		if (!is.na(sno.URG.sites$ew[i]) & !is.na(sno.URG.sites$sn[i])) {
   			snoIndex_Lev0[[as.character(sno.URG.sites$site_id[i])]] <- list(start=c(sno.URG.sites$ew[i], sno.URG.sites$sn[i], 1),
                                               end=c(sno.URG.sites$ew[i], sno.URG.sites$sn[i], 1), stat="mean")
   			snoIndex_Lev1[[as.character(sno.URG.sites$site_id[i])]] <- list(start=c(sno.URG.sites$ew[i], 1, sno.URG.sites$sn[i], 1),
                                          end=c(sno.URG.sites$ew[i], 1, sno.URG.sites$sn[i], 1), stat="mean")
   			snoIndex_Lev2[[as.character(sno.URG.sites$site_id[i])]] <- list(start=c(sno.URG.sites$ew[i], 2, sno.URG.sites$sn[i], 1),
                                               end=c(sno.URG.sites$ew[i], 2, sno.URG.sites$sn[i], 1), stat="mean")
   			snoIndex_Lev3[[as.character(sno.URG.sites$site_id[i])]] <- list(start=c(sno.URG.sites$ew[i], 3, sno.URG.sites$sn[i], 1),
                                               end=c(sno.URG.sites$ew[i], 3, sno.URG.sites$sn[i], 1), stat="mean")
   			snoIndex_Lev4[[as.character(sno.URG.sites$site_id[i])]] <- list(start=c(sno.URG.sites$ew[i], 4, sno.URG.sites$sn[i], 1),
                                               end=c(sno.URG.sites$ew[i], 4, sno.URG.sites$sn[i], 1), stat="mean")
   			}
 		}
	# Other met stations
        for (i in 1:length(met.URG.sites$site_id)) {
                if (!is.na(met.URG.sites$ew[i]) & !is.na(met.URG.sites$sn[i])) {
                        snoIndex_Lev0[[as.character(met.URG.sites$site_id[i])]] <- list(start=c(met.URG.sites$ew[i], met.URG.sites$sn[i], 1),
                                               end=c(met.URG.sites$ew[i], met.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev1[[as.character(met.URG.sites$site_id[i])]] <- list(start=c(met.URG.sites$ew[i], 1, met.URG.sites$sn[i], 1),
                                          end=c(met.URG.sites$ew[i], 1, met.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev2[[as.character(met.URG.sites$site_id[i])]] <- list(start=c(met.URG.sites$ew[i], 2, met.URG.sites$sn[i], 1),
                                               end=c(met.URG.sites$ew[i], 2, met.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev3[[as.character(met.URG.sites$site_id[i])]] <- list(start=c(met.URG.sites$ew[i], 3, met.URG.sites$sn[i], 1),
                                               end=c(met.URG.sites$ew[i], 3, met.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev4[[as.character(met.URG.sites$site_id[i])]] <- list(start=c(met.URG.sites$ew[i], 4, met.URG.sites$sn[i], 1),
                                               end=c(met.URG.sites$ew[i], 4, met.URG.sites$sn[i], 1), stat="mean")
                        }
                }
	# Ameriflux proxy sites
	for (i in 1:length(amf.URG.sites$site_id)) {
                if (!is.na(amf.URG.sites$ew[i]) & !is.na(amf.URG.sites$sn[i])) {
                        snoIndex_Lev0[[as.character(amf.URG.sites$site_id[i])]] <- list(start=c(amf.URG.sites$ew[i], amf.URG.sites$sn[i], 1),
                                               end=c(amf.URG.sites$ew[i], amf.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev1[[as.character(amf.URG.sites$site_id[i])]] <- list(start=c(amf.URG.sites$ew[i], 1, amf.URG.sites$sn[i], 1),
                                          end=c(amf.URG.sites$ew[i], 1, amf.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev2[[as.character(amf.URG.sites$site_id[i])]] <- list(start=c(amf.URG.sites$ew[i], 2, amf.URG.sites$sn[i], 1),
                                               end=c(amf.URG.sites$ew[i], 2, amf.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev3[[as.character(amf.URG.sites$site_id[i])]] <- list(start=c(amf.URG.sites$ew[i], 3, amf.URG.sites$sn[i], 1),
                                               end=c(amf.URG.sites$ew[i], 3, amf.URG.sites$sn[i], 1), stat="mean")
                        snoIndex_Lev4[[as.character(amf.URG.sites$site_id[i])]] <- list(start=c(amf.URG.sites$ew[i], 4, amf.URG.sites$sn[i], 1),
                                               end=c(amf.URG.sites$ew[i], 4, amf.URG.sites$sn[i], 1), stat="mean")
                        }
                }
        }


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# RTOUT PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runBasinRtout) {
 
         ## ------------------------------------------------------------------------
         # Setup RTOUT files
 
         filesList <- list.files(path=modoutPath, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=TRUE)
 	if (exists("modoutPath2")) {
 		filesList <- c(filesList, list.files(path=modoutPath2, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=TRUE))
 		nameList <- c(list.files(path=modoutPath, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=FALSE),
 				list.files(path=modoutPath2, pattern=glob2rx('*.RTOUT_DOMAIN*'), full.names=FALSE))
 		filesDf <- data.frame(x=rev(filesList), y=rev(nameList), stringsAsFactors=FALSE)
 		filesList <- as.list(rev(filesDf[!(duplicated(filesDf$y)),1]))
 		}
 	#filesList <- filesList[1:100]
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
 
         rtoutDF <- GetMultiNcdf(indexList=rtoutIndexList,
                         variableList=rtoutVariableList,
                         filesList=rtoutFilesList, parallel=TRUE )
         modRtout <- ReshapeMultiNcdf(rtoutDF)
         modRtout <- modRtout[order(modRtout$statArg, modRtout$POSIXct),]
 	#save(modRtout, file="test.Rdata")
 	assign(paste0("modRtout_", objTagList[j], "_BAS"), modRtout)
         saveList <- c(saveList, paste0("modRtout_", objTagList[j], "_BAS"))
         rm(rtoutDF, modRtout)
 	rm(filesList, rtoutFilesList, varNames, rtoutVars, rtoutVariableList, level0, rtoutInd, rtoutIndexList)
 	gc()
 	save(list=saveList, file=modOutImg)
 } # end rtout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# LDASOUT PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runBasinLdasout | runSnoLdasout) {
 
 	## ------------------------------------------------------------------------
 	# Setup LDASOUT files
 
 	filesList <- list.files(path=modoutPath, pattern=glob2rx('*.LDASOUT_DOMAIN*'), full.names=TRUE)
         if (exists("modoutPath2")) {
                 filesList <- c(filesList, list.files(path=modoutPath2, pattern=glob2rx('*.LDASOUT_DOMAIN*'), full.names=TRUE))
                 nameList <- c(list.files(path=modoutPath, pattern=glob2rx('*.LDASOUT_DOMAIN*'), full.names=FALSE),
                                 list.files(path=modoutPath2, pattern=glob2rx('*.LDASOUT_DOMAIN*'), full.names=FALSE))
                 filesDf <- data.frame(x=rev(filesList), y=rev(nameList), stringsAsFactors=FALSE)
                 filesList <- as.list(rev(filesDf[!(duplicated(filesDf$y)),1]))
                 }
 	#filesList <- filesList[1:100]
 	ldasoutFilesList <- list( ldasout = filesList)
 
 	## ------------------------------------------------------------------------
 	# Setup variables
 
 	if (varsLdasoutSUB) {
 		# SUBSET
 		varNames <- c('ACCECAN', 'ACCEDIR', 'ACCETRAN', 'ACCPRCP', 
 			'CANICE', 'CANLIQ',
 			'SFCRNOFF','UGDRNOFF',
 			 rep('SOIL_M',4),
 		 	rep('SOIL_T',4),
			rep('SOIL_W', 4),
 			'FIRA', 'FSA', 'GRDFLX', 'HFX', 'LH', 'LWFORC', 'SWFORC', 
 			'ACSNOM', 'ACSNOW', 'ALBEDO',
 			'SNEQV', 'FSNO', 'ISNOW', 'QSNOW', 'SNOWH',
 		 	rep('SNICE',3),
 			 rep('SNLIQ',3), 
 		 	rep('SNOW_T',3),
 		 	rep('ZSNSO_SN',3),
			'CHUC','LAI')
 		varLabels <- c('ACCECAN', 'ACCEDIR', 'ACCETRAN', 'ACCPRCP',
                 	'CANICE', 'CANLIQ',
                 	'SFCRNOFF','UGDRNOFF',
                  	paste0('SOIL_M',1:4),
                  	paste0('SOIL_T',1:4),
			paste0('SOIL_W',1:4),
                 	'FIRA', 'FSA', 'GRDFLX', 'HFX', 'LH', 'LWFORC', 'SWFORC', 
                 	'ACSNOM', 'ACSNOW', 'ALBEDO',
                 	'SNEQV', 'FSNO', 'ISNOW', 'QSNOW', 'SNOWH',
                  	paste0('SNICE',1:3),
                  	paste0('SNLIQ',1:3), 
                  	paste0('SNOW_T',1:3),
                  	paste0('ZSNSO_SN',1:3),
			'CHUC','LAI')
 		ldasoutVars <- as.list( varNames )
 		names(ldasoutVars) <- varLabels
 		ldasoutVariableList <- list( ldasout = ldasoutVars )
 		# INDEXES
 		if (runBasinLdasout) {
 			level0 <- basgeoIndex_Lev0
 			level1 <- basgeoIndex_Lev1
 			level2 <- basgeoIndex_Lev2
 			level3 <- basgeoIndex_Lev3
 			level4 <- basgeoIndex_Lev4
 			ldasoutBasInd <- list( level0, level0, level0, level0,
                         		level0, level0,
                         		level0, level0,
                         		level1, level2, level3, level4,
                         		level1, level2, level3, level4,
					level1, level2, level3, level4,
                         		level0, level0, level0, level0, level0, level0, level0,
                         		level0, level0, level0,
                         		level0, level0, level0, level0, level0,
                         		level1, level2, level3,
                         		level1, level2, level3,
                         		level1, level2, level3,
                         		level1, level2, level3,
					level0, level0 )
 			names(ldasoutBasInd) <- names(ldasoutVars)
 			ldasoutBasIndexList <- list( ldasout = ldasoutBasInd )
 			} # end runBasin
 		if (runSnoLdasout) {
 			level0 <- snoIndex_Lev0
 			level1 <- snoIndex_Lev1
 			level2 <- snoIndex_Lev2
 			level3 <- snoIndex_Lev3
 			level4 <- snoIndex_Lev4
 			ldasoutSnoInd <- list( level0, level0, level0, level0,
                 		        level0, level0,
                         		level0, level0,
                         		level1, level2, level3, level4,
                       	   		level1, level2, level3, level4,
					level1, level2, level3, level4,
                     		        level0, level0, level0, level0, level0, level0, level0,
                        		level0, level0, level0,
                     			level0, level0, level0, level0, level0,
                       			level1, level2, level3,
                			level1, level2, level3,
                         		level1, level2, level3,
                         		level1, level2, level3,
					level0, level0 )
 	        	names(ldasoutSnoInd) <- names(ldasoutVars)
                 	ldasoutSnoIndexList <- list( ldasout = ldasoutSnoInd )
 			} # end runSno
 
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
         	if (runBasinLdasout) {
                 	level0 <- basgeoIndex_Lev0
                 	level1 <- basgeoIndex_Lev1
                 	level2 <- basgeoIndex_Lev2
                 	level3 <- basgeoIndex_Lev3
                 	level4 <- basgeoIndex_Lev4
 			ldasoutBasInd <- list( level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
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
                 	names(ldasoutBasInd) <- names(ldasoutVars)
                 	ldasoutBasIndexList <- list( ldasout = ldasoutBasInd )
 			} # end runBasin
         	if (runSnoLdasout) {
                 	level0 <- snoIndex_Lev0
                		level1 <- snoIndex_Lev1
                 	level2 <- snoIndex_Lev2
                 	level3 <- snoIndex_Lev3
                 	level4 <- snoIndex_Lev4
                 	ldasoutSnoInd <- list( level0, level0, level0, level0, level0, level0, level0, level0, level0, level0,
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
                 	names(ldasoutSnoInd) <- names(ldasoutVars)
                 	ldasoutSnoIndexList <- list( ldasout = ldasoutSnoInd )
                 	} # end runSno
 	} # end ifelse vars subset
 
 	## ------------------------------------------------------------------------
 	# Get data and flatten files
 	## ------------------------------------------------------------------------
 	if (runBasinLdasout) {
 		# BASIN
 		ldasoutDF <- GetMultiNcdf(indexList=ldasoutBasIndexList,
                           variableList=ldasoutVariableList,
                           filesList=ldasoutFilesList,
                           parallel=TRUE )
 		modLdasout <- ReshapeMultiNcdf(ldasoutDF)
 		modLdasout <- CalcNoahmpFluxes(modLdasout, "statArg")
 		assign(paste0("modLdasout_", objTagList[j], "_BAS"), modLdasout)
 		saveList <- c(saveList, paste0("modLdasout_", objTagList[j], "_BAS"))
 		rm(ldasoutDF, modLdasout)
 		gc()
 		save(list=saveList, file=modOutImg)
 	}	
 
 	if (runSnoLdasout) {
 		# SNOTEL
 		ldasoutDF <- GetMultiNcdf(indexList=ldasoutSnoIndexList, 
                           variableList=ldasoutVariableList, 
                           filesList=ldasoutFilesList, 
                           parallel=TRUE )
 		modLdasout <- ReshapeMultiNcdf(ldasoutDF)
 		modLdasout <- CalcNoahmpFluxes(modLdasout, "statArg")
 		assign(paste0("modLdasout_", objTagList[j], "_SNO"), modLdasout)
 		saveList <- c(saveList, paste0("modLdasout_", objTagList[j], "_SNO"))
 		rm(ldasoutDF, modLdasout)
 		gc()
 		save(list=saveList, file=modOutImg)
 		}

 	} # end ldasout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# GW PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runGwout) {
 	modGwout <- ReadGwOut(paste0(modoutPath, '/GW_outflow.txt'))
 	if (exists("modoutPath2")) {
 		modGwout <- rbind(modGwout, ReadGwOut(paste0(modoutPath2, '/GW_outflow_PART2.txt')))
		} 
	# Filter out non-unique dates. Take values from latest run if dups.
 	modGwout <- modGwout[nrow(modGwout):1,]
 	modGwout$uni <- paste(modGwout$basin, modGwout$timest, sep=",")
 	modGwout <- modGwout[!(duplicated(modGwout$uni)),]
 	modGwout$uni <- NULL
 	modGwout <- modGwout[nrow(modGwout):1,]
 	assign(paste0("modGwout_", objTagList[j]), modGwout)
        saveList <- c(saveList, paste0("modGwout_", objTagList[j]))
        rm(modGwout)
} # end gwout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# FRXST PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runFrxstout) {
        modFrxstout <- ReadFrxstPts(paste0(modoutPath, '/frxst_pts_out.txt'))
 	if (exists("modoutPath2")) {
 		modFrxstout <- rbind(modFrxstout, ReadFrxstPts(paste0(modoutPath2, '/frxst_pts_out_PART2.txt')))
		}       
        # Filter out non-unique dates. Take values from latest run if dups.
        modFrxstout <- modFrxstout[nrow(modFrxstout):1,]
        modFrxstout$uni <- paste(modFrxstout$st_id, modFrxstout$timest, sep=",")
        modFrxstout <- modFrxstout[!(duplicated(modFrxstout$uni)),]
        modFrxstout$uni <- NULL
        modFrxstout <- modFrxstout[nrow(modFrxstout):1,]
        assign(paste0("modFrxstout_", objTagList[j]), modFrxstout)
        saveList <- c(saveList, paste0("modFrxstout_", objTagList[j]))
        rm(modFrxstout)
} # end frxstout processing


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# LDASIN PROCESSING
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

if (runSnoLdasin) {
 
 	## ------------------------------------------------------------------------
 	# Setup LDASIN files
 
 	filesList <- list.files(path=forcPath, pattern=glob2rx('2014*.LDASIN_DOMAIN*'), full.names=TRUE)
 	filesList <- c(filesList, list.files(path=forcPath, pattern=glob2rx('2015*.LDASIN_DOMAIN*'), full.names=TRUE))
	#filesList <- filesList[1:100]
 	ldasinFilesList <- list( ldasin = filesList)
 
 
 	## ------------------------------------------------------------------------
 	# Setup variables
 
 	varNames <- c('T2D', 'Q2D', 'U2D', 'V2D', 'PSFC', 'SWDOWN', 'LWDOWN','RAINRATE')
 	ldasinVars <- as.list( varNames )
 	names(ldasinVars) <- varNames
 	ldasinVariableList <- list( ldasin = ldasinVars )
 
	# SNOTEL
	if (runSnoLdasin) {
 		## ------------------------------------------------------------------------
 		# Setup indexes
 
 		level0 <- snoIndex_Lev0
 		ldasinInd <- list( level0, level0, level0, level0, level0, level0, level0, level0 )
 		names(ldasinInd) <- names(ldasinVars)
 		ldasinIndexList <- list( ldasin = ldasinInd )
 
 		## ------------------------------------------------------------------------
 		# Get data and flatten files
 
 		ldasinDF <- GetMultiNcdf(indexList=ldasinIndexList,
                         variableList=ldasinVariableList,
                         filesList=ldasinFilesList, parallel=TRUE )
 		modLdasin <- ReshapeMultiNcdf(ldasinDF)
 		modLdasin <- modLdasin[order(modLdasin$statArg, modLdasin$POSIXct),]
 		assign(paste0("modLdasin_", objTagList[j], "_SNO"), modLdasin)
 		saveList <- c(saveList, paste0("modLdasin_", objTagList[j], "_SNO"))
 		rm(ldasinDF, modLdasin)
 		rm(filesList, ldasinFilesList, varNames, ldasinVars, ldasinVariableList, level0, ldasinInd, ldasinIndexList)
 		gc()
 		#save(list=saveList, file=modOutImg)
		}


 }

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# SAVE DATA
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

save(list=saveList, file=modOutImg)

##quit("no")
##proc.time()
}
