mskgeoPath <- '../DOMAIN/MASKS'
fileList <- list.files(mskgeoPath, pattern=glob2rx("*.tif"), full.names=TRUE)
nameList <- list.files(mskgeoPath, pattern=glob2rx("*.tif"))

file2gageList <- list("alamosa_1k"="ALATERCO", 
               "conejos_1k"="CONMOGCO", 
               "s_frk_1k"="RIOSFKCO", 
               "rio_wagw_1k"="RIOWAGCO", 
               "saguache_1k"="SAGSAGCO", 
               "trinchera_1k"="TRITURCO",
               "pinos_1k"="08248000",
               "plat_inflo_1k"="PLATERO_INFLOW",
               "rio_deln_1k"="RIODELCO",
               "sn_anton_1k"="SANORTCO")

stats_lstnight_basins <- data.frame()
for (i in 1:length(fileList)) {
	nm <- unlist(strsplit(nameList[i],"[.]"))[1]
	msk <- raster(fileList[i])
	msk[msk==-9999] <- NA
	tmp <- mask(x=lst_night.b, mask=msk)
	stats <- CalcStatsRS(tmp)
	stats$basin_id <- file2gageList[[nm]]
	stats_lstnight_basins <- rbind(stats_lstnight_basins, stats)
}