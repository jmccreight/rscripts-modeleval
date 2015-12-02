###################################################
##                GENERATE PLOTS                 ##
###################################################


################## General Setup ##################

dir.create(writePlotDir, showWarnings = FALSE)
source("PlotBasin.R")
source("PlotSnotel.R")
source("PlotMaps.R")
library(ggplot2)

lineColors <- c("dodgerblue", "darkorange1", "olivedrab", "chocolate", "darkmagenta")
lineTyp <- 1
lineWd <- 3

# Get needed geo info
ncid <- ncdf4::nc_open(geoFile)
geoDX <- ncdf4::ncatt_get(ncid,varid=0,'DX')$value
ncdf4::nc_close(ncid)
hydDX <- geoDX/aggfact

if (!exists("gageList") & exists("rtLinks")) {
	gageList<-subset(rtLinks[c("link","site_no")], rtLinks$gages!="")
	#names(link2gage)[names(link2gage)=="gages"]<-"site_no"
}

if (writeHtml) {
	library(knitr)
	library(pander)
	library(xtable)
	if (accflowPlot | hydroPlot | flowlsmPlot) {
		writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', con=paste0(writePlotDir,"/plots_hydro.Rmd"))
		cat('# MODEL OUTPUT: HYDROLOGY\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
	}
	if (accprecipPlot) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', con=paste0(writePlotDir,"/plots_climate.Rmd"))
                cat('# MODEL OUTPUT: CLIMATE\n', file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
	}
	if (flowswePlot | swePlot) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', con=paste0(writePlotDir,"/plots_snow.Rmd"))
                cat('# MODEL OUTPUT: SNOW\n', file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
	}
        if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', con=paste0(writePlotDir,"/plots_stats.Rmd"))
                cat('# MODEL OUTPUT: STATS\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
        }
}

## -----------------------------------------------------------------------
# Generate Plots

if (accprecipPlot | flowswePlot | flowlsmPlot) {
	modLdasout_BAS <- list(native=subset(modLdasout[["native"]], modLdasout[["native"]]$fileGroup=="ldasout.basgeo"),
                       snoday=subset(modLdasout[["snoday"]], modLdasout[["snoday"]]$fileGroup=="ldasout.basgeo"),
                       utcday=subset(modLdasout[["utcday"]], modLdasout[["utcday"]]$fileGroup=="ldasout.basgeo"))
}

if (swePlot) {
	modLdasout_SNO <- list(native=subset(modLdasout[["native"]], modLdasout[["native"]]$fileGroup=="ldasout.sno"),
                       snoday=subset(modLdasout[["snoday"]], modLdasout[["snoday"]]$fileGroup=="ldasout.sno"),
                       utcday=subset(modLdasout[["utcday"]], modLdasout[["utcday"]]$fileGroup=="ldasout.sno"))
}


## -----------------------------------------------------------------------
# Generate Plots

# Accumulated Flow
if (accflowPlot) {
message("Generating accumulated flow plots...")
# Setup
accflowList <- list()
if (is.null(accflowTags)) accflowTags <- unique(modFrxstout$tag)
for (i in 1:length(accflowTags)) {
        accflowList[[i]] <- subset(modFrxstout, modFrxstout$tag==accflowTags[i])
}
accflowColors <- lineColors[1:length(accflowList)]
accflowTypes <- rep(lineTyp, length(accflowList))
accflowWidths <- rep(lineWd, length(accflowList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/accstrflow_", n, ".png"), width=2100, height=1350, res=225)
        PlotAccFlow(n, modDfs=accflowList,
                        obs=obsStrData.dy,
                        stdate=accflowStartDate,
                        enddate=accflowEndDate,
                        labMods=accflowTags,
                        labObs="Observed",
                        lnCols=accflowColors,
                        lnTyps=accflowTypes, lnWds=accflowWidths,
                        labTitle=paste0("Accumulated Flow: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")"), obsCol="cumqvol_mm")
        dev.off()
}
if (writeHtml) {
	cat('## Accumulated Flow Plots\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
	for (n in names(gage2basinList)) {
		cat(paste0("```{r accflow_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
		plottxt <- knitr::knit_expand(text='PlotAccFlow("{{n}}", modDfs=accflowList,
                        obs=obsStrData.dy,
                        stdate=accflowStartDate,
                        enddate=accflowEndDate,
                        labMods=accflowTags,
                        labObs="Observed",
                        lnCols=accflowColors,
                        lnTyps=accflowTypes, lnWds=accflowWidths,
                        labTitle=paste0("Accumulated Flow: ", "{{n}}", " (", obsStrMeta$site_name[obsStrMeta$site_no=="{{n}}"], ")"), obsCol="cumqvol_mm")\n')
		cat(plottxt, file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
		cat('```\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
	}
}
}

# Hydrographs
if (hydroPlot) {
message("Generating hydrograph plots...")
# Setup
hydroList <- list()
if (is.null(hydroTags)) {
	if (exists("modFrxstout")) {
		hydroTags <- unique(modFrxstout$tag)
		gageNames <- names(gage2basinList)
		idCol <- "site_no"
	} else if (exists("modChrtout")) {
		hydroTags <- unique(modChrtout$tag)
		if (exists("gageList")) {
			gageNames <- unique(gageList$link)
		} else {
			gageNames <- unique(obsStrData$link)
		}
		idCol <- "link"
	}
}
for (i in 1:length(hydroTags)) {
	if (exists("modFrxstout")) {
        	hydroList[[i]] <- subset(modFrxstout, modFrxstout$tag==hydroTags[i])
	} else if (exists("modChrtout")) {
		hydroList[[i]] <- subset(modChrtout, modChrtout$tag==hydroTags[i])
	} else {
		stop()
	}
}
hydroColors <- lineColors[1:length(hydroList)]
hydroTypes <- rep(lineTyp, length(hydroList))
hydroWidths <- rep(lineWd, length(hydroList))
# Loop plots
for (n in gageNames) {
        png(paste0(writePlotDir, "/hydrogr_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlow(n, modDfs=hydroList,
                        obs=obsStrData,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=paste0("Streamflow: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")"),
                        stdate=hydroStartDate, enddate=hydroEndDate, obsCol="q_cms", idCol=idCol)
        dev.off()
}
if (writeHtml) {
        cat('## Hydrographs\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
        for (n in gageNames) {
                cat(paste0("```{r hydro_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotFlow("{{n}}", modDfs=hydroList,
                        obs=obsStrData,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=paste0("Streamflow: ", "{{n}}", " (", obsStrMeta$site_name[obsStrMeta$site_no=="{{n}}"], ")"),
                        stdate=hydroStartDate, enddate=hydroEndDate, obsCol="q_cms")\n')
                cat(plottxt, file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
        }
}
}

# Accumulated Precip
if (accprecipPlot) {
message("Generating accumulated precip plots...")
# Setup
accprecipList <- list()
if (is.null(accprecipTags)) accprecipTags <- unique(modLdasout_BAS[["native"]]$tag)
for (i in 1:length(accprecipTags)) {
        accprecipList[[i]] <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==accprecipTags[i])
}
accprecipColors <- lineColors[1:length(accprecipList)]
accprecipTypes <- rep(lineTyp, length(accprecipList))
accprecipWidths <- rep(lineWd, length(accprecipList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/accprecip_", n, ".png"), width=2100, height=1350, res=225)
        PlotAccPrecip(n, modDfs=accprecipList,
                        stdate=accprecipStartDate,
                        enddate=accprecipEndDate,
                        labMods=accprecipTags,
                        lnCols=accprecipColors,
                        lnTyps=accprecipTypes, lnWds=accprecipWidths,
                        labTitle=paste0("Accumulated Precip: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")"))
        dev.off()
}
if (writeHtml) {
        cat('## Accumulated Basin-Mean Precip Plots\n', file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r accprecip_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
				file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotAccPrecip("{{n}}", modDfs=accprecipList,
                        stdate=accprecipStartDate,
                        enddate=accprecipEndDate,
                        labMods=accprecipTags,
                        lnCols=accprecipColors,
                        lnTyps=accprecipTypes, lnWds=accprecipWidths,
                        labTitle=paste0("Accumulated Precip: ", "{{n}}", " (", obsStrMeta$site_name[obsStrMeta$site_no=="{{n}}"], ")"))\n')
                cat(plottxt, file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
        }
}
}

# Flow and basin mean SWE
if (flowswePlot) {
message("Generating flow + basin SWE plots...")
# Setup
flowsweStrList <- list()
flowsweLsmList <- list()
if (is.null(flowsweTags)) flowsweTags <- unique(modLdasout_BAS[["native"]]$tag)
for (i in 1:length(flowsweTags)) {
        flowsweStrList[[i]] <- subset(modFrxstout, modFrxstout$tag==flowsweTags[i])
        flowsweLsmList[[i]] <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowsweTags[i])
}
flowsweColors <- lineColors[1:length(flowsweStrList)]
flowsweTypes <- rep(lineTyp, length(flowsweStrList))
flowsweWidths <- rep(lineWd, length(flowsweStrList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/flowswe_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlowSwe(n, modDfs=flowsweStrList,
                        lsmDfs=flowsweLsmList,
                        obs=obsStrData.dy,
                        labMods=flowsweTags,
                        lnCols=flowsweColors,
                        lnWds=flowsweWidths,
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n),
                        stdate=flowsweStartDate, enddate=flowsweEndDate)
        dev.off()
}
if (writeHtml) {
        cat('## Streamflow & Basin-mean SWE Plots\n', file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r flowswe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotFlowSwe("{{n}}", modDfs=flowsweStrList,
                        lsmDfs=flowsweLsmList,
                        obs=obsStrData.dy,
                        labMods=flowsweTags,
                        lnCols=flowsweColors,
                        lnWds=flowsweWidths,
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", "{{n}}"),
                        stdate=flowsweStartDate, enddate=flowsweEndDate)\n')
                cat(plottxt, file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
        }
}
}

# Flow and basin mean LSM runoff
if (flowlsmPlot) {
message("Generating flow + basin runoff plots...")
# Setup
if (is.null(flowlsmTags)) {
        if (exists("modFrxstout")) {
                flowlsmTags <- unique(modFrxstout$tag)
                gageNames <- names(gage2basinList)
                idCol <- "site_no"
        } else if (exists("modChrtout")) {
                flowlsmTags <- unique(modChrtout$tag)
                if (!is.null(gageList)) {
                        gageNames <- unique(gageList$link)
                } else {
                        gageNames <- unique(obsStrData$link)
                }
                idCol <- "link"
        }
}
flowlsmColors <- lineColors[1:2]
flowlsmTypes <- rep(lineTyp, 2)
flowlsmWidths <- rep(lineWd, 2)
# Loop plots
for (i in 1:length(flowlsmTags)) {
        if (exists("modFrxstout")) {
                strDf <- subset(modFrxstout, modFrxstout$tag==flowlsmTags[i])
        } else if (exists("modChrtout")) {
                strDf <- subset(modChrtout, modChrtout$tag==flowlsmTags[i])
        } else {
                stop()
        }
        lsmDf <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowlsmTags[i])
	ts <- as.integer(difftime(lsmDf$POSIXct[2],lsmDf$POSIXct[1], units="secs"))
	for (n in gageNames) {
		ngage <- ifelse(idCol=="link", as.integer(subset(gageList$site_no, gageList$link==n)), n)
		ngageChar <- ifelse(nchar(as.character(ngage))==7, paste0("0", as.character(ngage)), as.character(ngage))
        	png(paste0(writePlotDir, "/flowlsm_", n, ".png"), width=2100, height=1350, res=225)
        	PlotFlowLsm(n, modDf=strDf, lsmDf=lsmDf, 
                        obs=obsStrData,
                        labMods=flowlsmTags,
                        labObs="Observed",
                        lnCols=flowlsmColors,
                        lnWds=flowlsmWidths,
                        labTitle=paste0("Streamflow: ", ngageChar, " (", obsStrMeta$site_name[obsStrMeta$site_no==ngageChar], ")"),
                        stdate=flowlsmStartDate, enddate=flowlsmEndDate, obsCol="q_cms", idCol=idCol,
			tsSecs=ts, areaSqKm=mskgeo.areaList[[as.character(ngage)]]*geoDX/1000, ngage=ngage)
        	dev.off()
	}
}
if (writeHtml) {
for (i in 1:length(flowlsmTags)) {
        if (exists("modFrxstout")) {
                strDf <- subset(modFrxstout, modFrxstout$tag==flowlsmTags[i])
        } else if (exists("modChrtout")) {
                strDf <- subset(modChrtout, modChrtout$tag==flowlsmTags[i])
        } else {
                stop()
        }
        lsmDf <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowlsmTags[i])
        ts <- as.integer(difftime(lsmDf$POSIXct[2],lsmDf$POSIXct[1], units="secs"))
        cat(paste0('## Streamflow & Basin-mean LSM Runoff Plots:', i, '\n'), file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
        for (n in gageNames) {
                cat(paste0("```{r flowlsm_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"),
                        file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='ngage <- ifelse(idCol=="link", 
					as.integer(subset(gageList$site_no, gageList$link=="{{n}}")), "{{n}}");
			ngageChar <- ifelse(nchar(as.character(ngage))==7, paste0("0", as.character(ngage)), as.character(ngage));
			PlotFlowLsm("{{n}}", modDf=strDf,
			lsmDf=lsmDf,
                        obs=obsStrData,
                        labMods=flowlsmTags,
                        labObs="Observed",
                        lnCols=flowlsmColors,
                        lnWds=flowlsmWidths,
                        labTitle=paste0("Streamflow: ", ngageChar, " (", obsStrMeta$site_name[obsStrMeta$site_no==ngageChar], ")"),
                        stdate=flowlsmStartDate, enddate=flowlsmEndDate, obsCol="q_cms", idCol=idCol,
			tsSecs=ts, areaSqKm=mskgeo.areaList[[as.character(ngage)]]*geoDX/1000, ngage=ngage)\n')
                cat(plottxt, file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/plots_hydro.Rmd"), append=TRUE)
        }
}
}
}
   
# SWE
if (swePlot) {
message("Generating SWE plots...")
# Setup
sweList <- list()
if (is.null(sweTags)) sweTags <- unique(modLdasout_SNO[["native"]]$tag)
for (i in 1:length(sweTags)) {
        sweList[[i]] <- subset(modLdasout_SNO[["native"]], modLdasout_SNO[["native"]]$tag==sweTags[i])
}
sweColors <- lineColors[1:length(sweList)]
sweTypes <- rep(lineTyp, length(sweList))
sweWidths <- rep(lineWd, length(sweList))
# Loop plots
sites <- unique(modLdasout_SNO[["native"]]$statArg)
for (n in sites) {
  png(paste0(writePlotDir, "/swe_", n, ".png"), width=2100, height=1350, res=225)
  PlotSwe(n, modDfs=sweList,
                obs=obsSnoData, obsmeta=obsSnoMeta,
                labMods=sweTags,
                lnCols=sweColors,
                lnWds=sweWidths,
                precCol.obs="CumPrec_mm", precCol.mod="ACCPRCP",
                sweCol.obs="SWE_mm", sweCol.mod="SNEQV", fact=1, snowh=FALSE,
                labTitle="Accumulated Precipitation and SWE",
                stdate=sweStartDate, enddate=sweEndDate)
  dev.off()
}
if (writeHtml) {
        cat('## Station SWE Plots\n', file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
        for (n in sites) {
                cat(paste0("```{r swe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotSwe("{{n}}", modDfs=sweList,
                	obs=obsSnoData, obsmeta=obsSnoMeta,
                	labMods=sweTags,
                	lnCols=sweColors,
                	lnWds=sweWidths,
                	precCol.obs="CumPrec_mm", precCol.mod="ACCPRCP",
                	sweCol.obs="SWE_mm", sweCol.mod="SNEQV", fact=1, snowh=FALSE,
                	labTitle="Accumulated Precipitation and SWE",
                	stdate=sweStartDate, enddate=sweEndDate)\n')
                	cat(plottxt, file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
                	cat('```\n', file=paste0(writePlotDir,"/plots_snow.Rmd"), append=TRUE)
        }
}
}

# MET
if (metPlot) {
	message("Generating MET plots...")
	# Setup
	modLdasin_MET <- subset(modLdasin[["utcday"]], modLdasin[["utcday"]]$fileGroup=="ldasin.met")
	if (is.null(metTags)) metTags <- unique(modLdasin_MET$tag)
	metSites <- unique(modLdasin_MET$statArg)
	for (i in metTags) {
		modLdasin_MET_TAG <- subset(modLdasin_MET, modLdasin_MET$tag==i)
		# Loop Sites
		for (n in metSites) {
  			print(n)
  			# Temperature
  			png(paste0(writePlotDir, "/met_temp_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("Tmean_K", "Tmax_K", "Tmin_K"),
                        	modVars=c("T2D_mean", "T2D_max", "T2D_min"),
                        	lnLabs=c("Mean Temp (C)", "Max Temp (C)", "Min Temp (C)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Temperature"),
                        	xLab="", adj=(-273.15), 
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
			if (writeHtml) {
			        cat('## MET Station Plots\n', file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
                		cat(paste0("```{r met_", i, "_", n, ", fig.width = 13.5, fig.height = 21, out.width='1350', out.height='2100', echo=FALSE}\n"),
                        		file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
                		plottxt <- knitr::knit_expand(text='PlotMet(obs=obsMetData.dy,
                                	mod=modLdasin_MET_TAG,
                                	site="{{n}}",
                                	obsVars=c("Tmean_K", "Tmax_K", "Tmin_K"),
                                	modVars=c("T2D_mean", "T2D_max", "T2D_min"),
                                	lnLabs=c("Mean Temp (C)", "Max Temp (C)", "Min Temp (C)"),
                                	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id=="{{n}}"], ":\nDaily Temperature"),
                                	xLab="", adj=(-273.15), 
                                	stdate=metStartDate, enddate=metEndDate)\n')
                        	cat(plottxt, file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
                        	cat('```\n', file=paste0(writePlotDir,"/plots_climate.Rmd"), append=TRUE)
        		}
  			# SW Radiation
  			png(paste0(writePlotDir, "/met_swrad_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("SWRad_mean", "SWRad_max", "SWRad_min"),
                        	modVars=c("SWDOWN_mean", "SWDOWN_max", "SWDOWN_min"),
                        	lnLabs=c("Mean Rad (W/m2)", "Max Rad (W/m2)", "Min Rad (W/m2)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Shortwave Radiation"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Wind
  			png(paste0(writePlotDir, "/met_wind_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("Wind_mean", "Wind_max", "Wind_min"),
                        	modVars=c("Wind_mean", "Wind_max", "Wind_min"),
                        	lnLabs=c("Mean Speed (m/s)", "Max Speed (m/s)", "Min Speed (m/s)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Wind Speed"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Humidity
  			png(paste0(writePlotDir, "/met_relhum_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("RH_mean", "RH_max", "RH_min"),
                        	modVars=c("RelHum_mean", "RelHum_max", "RelHum_min"),
                        	lnLabs=c("Mean RH (0-1)", "Max RH (0-1)", "Min RH (0-1)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nRelative Humidity"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Pressure
  			png(paste0(writePlotDir, "/met_press_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("SurfPressmean_Pa", "SurfPressmax_Pa", "SurfPressmin_Pa"),
                        	modVars=c("PSFC_mean", "PSFC_max", "PSFC_min"),
                        	lnLabs=c("Mean Press (kPa)", "Max Press (kPa)", "Min Press (kPa)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nSurface Pressure"),
                        	xLab="", mult=0.001,
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
		}
	}
}

# MAPS

# Initialize for maps
if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap ) {
	library(ggplot2)
	library(ggmap)
	library(gridExtra)
	if (reachRting) {
		modStrout <- modChrtout
	} else {
		modStrout <- modFrxstout
	}
	# Setup date ranges
	stdate_stats_PRINT <- ifelse(is.null(stdate_stats), min(modStrout$POSIXct), stdate_stats)
	enddate_stats_PRINT <- ifelse(is.null(enddate_stats), max(modStrout$POSIXct), enddate_stats)
	stdate_stats_sub_PRINT <- ifelse(is.null(stdate_stats_sub), min(modStrout$POSIXct), stdate_stats_sub)
	enddate_stats_sub_PRINT <- ifelse(is.null(enddate_stats_sub), max(modStrout$POSIXct), enddate_stats_sub)
	statsDateList <- list("Full" = paste0(format(as.POSIXct(stdate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"), 
			" to ", format(as.POSIXct(enddate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")), 
			"Sub" = paste0(format(as.POSIXct(stdate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"), 
                        " to ", format(as.POSIXct(enddate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")))
	# Setup map
	geoMap <- SetupMap(geoFile)
}

# STRFLOW Bias Maps
if (strBiasMap) {
	message("Generating STRFLOW Bias error map...")
	# Setup
	if (is.null(strBiasTags)) strBiasTags <- unique(stats_str$tag)
	if (is.null(strBiasSeas)) strBiasSeas <- unique(stats_str$seas)
	if (writeHtml) {
        	cat('## Streamflow Bias Maps\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
		strBias.ggList <- list()
		strBias.tblList <- list()
	}
	for (i in strBiasTags) {
        	for (j in strBiasSeas) {
                	gg <- PlotMapErrors(geoMap, stats_str,
                        	statsTag=i, statsVar=NULL, statsSeas=j,
                        	plotTitle="Modeled Streamflow Bias at CODWR Gages",
				plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        	sizeVar="t_mae", colorVar="t_bias",
                        	sizeLab="Mean Abs Error (cms)", colorLab="Bias (%)",
				minThreshSize=0, maxThreshSize=100,
				minThreshCol=(-100), maxThreshCol=100,
				minPtsize=2, maxPtsize=8,
				colBreaks=c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"), 
                        	valBreaks=c((-1000),(-50),(-10),10,50,1000))
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/str_bias_map_", i, "_", j, ".png"),
                        	plot=gg, units="in", width=8, height=6, dpi=100)
			if (writeHtml) {
				strBias.ggList <- c(strBias.ggList, list(gg))
				if (statsMapTables) {
					tbltmp <- subset(stats_str, stats_str$tag==i & stats_str$seas==j)
					tbltmp <- tbltmp[order(tbltmp$site_no),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_no, lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
					strBias.tblList <- c(strBias.tblList, list(tbltmp))
				}
			}
		}
	}
	#i <- 1
	if (writeHtml) {
		for (i in 1:length(strBias.ggList)) {
		#while (i <= length(ggList)) {
			#if ( (i+1) <= length(ggList)) {
			# Map
			cat(paste0("```{r strbiasmap_", i, ", fig.width = 12, fig.height = 9, echo=FALSE}\n"), 
				file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			#plottxt <- knitr::knit_expand(text='grid.arrange(ggList[[{{i}}]], ggList[[{{i}}+1]], ncol=2)\n')
			plottxt <- knitr::knit_expand(text='strBias.ggList[[{{i}}]]\n')
			#} else {
			#	cat(paste0("```{r strbiasmap_", i, ", fig.width = 10, fig.height = 5, out.width='1000', out.height='500', echo=FALSE}\n"),   
                	#                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                	#        plottxt <- knitr::knit_expand(text='grid.arrange(ggList[[{{i}}]], ncol=1)\n')
			#}
			cat(plottxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			# Table
			if (statsMapTables) {
				cat(paste0("```{r strbiastbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
				#tbltxt <- knitr::knit_expand(text='pandoc.table(strBias.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(strBias.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
				cat(tbltxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			}
			#i <- i + 2
		}
           }
}

# STRFLOW Correlation Maps
if (strCorrMap) {
	message("Generating STRFLOW Corr error map...")
	# Setup
	if (is.null(strCorrTags)) strCorrTags <- unique(stats_str$tag)
	if (is.null(strCorrSeas)) strCorrSeas <- unique(stats_str$seas)
	if (writeHtml) {
        	cat('## Streamflow Correlation Maps\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
        	strCorr.ggList <- list()
		strCorr.tblList <- list()
	}
	for (i in strCorrTags) {
        	for (j in strCorrSeas) {
                	gg <- PlotMapErrors(geoMap, stats_str,
                        	statsTag=i, statsVar=NULL, statsSeas=j,
                        	plotTitle="Modeled Streamflow Correlation at CODWR Gages",
				plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        	sizeVar="dy_cor", colorVar="dy_cor",
                        	sizeLab="Correlation", colorLab="Correlation",
				colorLow="orange", colorMid="yellow", colorHigh="cyan4",
				minThreshSize=0, maxThreshSize=1,
                                minThreshCol=0, maxThreshCol=1,
				minPtsize=0.5, maxPtsize=6,
                                colBreaks=c("#f7f7f7", "#ffffcc", "#c2e699", "#78c679", "#238443"),
                                valBreaks=c(-1, 0.2, 0.4, 0.6, 0.8, 1.0))
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/str_corr_map_", i, "_", j, ".png"),
                        	plot=gg, units="in", width=8, height=6, dpi=100)
                	if (writeHtml) {
                        	strCorr.ggList <- c(strCorr.ggList, list(gg))
				if (statsMapTables) {
                                	tbltmp <- subset(stats_str, stats_str$tag==i & stats_str$seas==j)
					tbltmp <- tbltmp[order(tbltmp$site_no),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_no, lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
						bias=tbltmp$t_bias, mae=tbltmp$t_mae,
						corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor, 
						nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	strCorr.tblList <- c(strCorr.tblList, list(tbltmp))
                		}
			}
          	}
   	}
	if (writeHtml) {
        	for (i in 1:length(strCorr.ggList)) {
			# Map
                	cat(paste0("```{r strcorrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='strCorr.ggList[[{{i}}]]\n')
                	cat(plottxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			# Table
			if (statsMapTables) {
                        	cat(paste0("```{r strcorrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                	file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(strCorr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(strCorr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			}
                }
        }
}

# SNOTEL SWE Maps
if (snosweErrMap) {
	message("Generating SNOTEL SWE error map...")
	# Setup
	if (is.null(snosweErrTags)) snosweErrTags <- unique(stats_ldasout_sno$tag)
	if (is.null(snosweErrSeas)) snosweErrSeas <- unique(stats_ldasout_sno$seas)
        if (writeHtml) {
                cat('## SNOTEL SWE Error Maps\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                snosweErr.ggList <- list()
		snosweErr.tblList <- list()
        }
	for (i in snosweErrTags) {
		for (j in snosweErrSeas) {
			gg <- PlotMapErrors(geoMap, stats_ldasout_sno,
                        	statsTag=i, statsVar="SWE", statsSeas=j,
                        	plotTitle="Modeled SWE Errors at SNOTEL Stations",
				plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        	sizeVar="t_mae", colorVar="t_msd",
                        	sizeLab="Mean Absolute Error (mm)", colorLab="Mean Signed Deviation (mm)")
			ggplot2::ggsave(filename=paste0(writePlotDir, "/sno_sweerr_map_", i, "_", j, ".png"), 
				plot=gg, units="in", width=9, height=6, dpi=100)
                        if (writeHtml) {
                                snosweErr.ggList <- c(snosweErr.ggList, list(gg))
                                tbltmp <- subset(stats_ldasout_sno, stats_ldasout_sno$tag==i & stats_ldasout_sno$seas==j & stats_ldasout_sno$var=="SWE")
				tbltmp <- tbltmp[order(tbltmp$site_id),]
                                tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$site_name, 
						lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                snosweErr.tblList <- c(snosweErr.tblList, list(tbltmp))
                        }
		}
	}
        if (writeHtml) {
                for (i in 1:length(snosweErr.ggList)) {
			# Map
                        cat(paste0("```{r snosweerrmap_", i, ", fig.width = 9, fig.height = 6, out.width='900', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snosweErr.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        # Table
                        cat(paste0("```{r snosweerrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        #tbltxt <- knitr::knit_expand(text='pandoc.table(snosweErr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
			tbltxt <- knitr::knit_expand(text='print(xtable(snosweErr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        cat(tbltxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                }
        }
}

# SNOTEL Precip Maps
if (snoprecipErrMap) {
	message("Generating SNOTEL precip error map...")
	# Setup
	if (is.null(snoprecipErrTags)) snoprecipErrTags <- unique(stats_ldasout_sno$tag)
	if (is.null(snoprecipErrSeas)) snoprecipErrSeas <- unique(stats_ldasout_sno$seas)
        if (writeHtml) {
                cat('## SNOTEL Precipitation Error Maps\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                snoprecipErr.ggList <- list()
		snoprecipErr.tblList <- list()
        }
	for (i in snoprecipErrTags) {
        	for (j in snoprecipErrSeas) {
                	gg <- PlotMapErrors(geoMap, stats_ldasout_sno,
                        	statsTag=i, statsVar="Precip", statsSeas=j,
                        	plotTitle="Modeled Precipitation Errors at SNOTEL Stations",
				plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        	sizeVar="t_mae", colorVar="t_msd",
                        	sizeLab="Mean Absolute Error (mm)", colorLab="Mean Signed Deviation (mm)")
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/sno_preciperr_map_", i, "_", j, ".png"),
                        	plot=gg, units="in", width=9, height=6, dpi=100)
                        if (writeHtml) {
                                snoprecipErr.ggList <- c(snoprecipErr.ggList, list(gg))
                                tbltmp <- subset(stats_ldasout_sno, stats_ldasout_sno$tag==i & stats_ldasout_sno$seas==j & stats_ldasout_sno$var=="Precip")
				tbltmp <- tbltmp[order(tbltmp$site_id),]
                                tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$site_name, 
                                                lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n,
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                snoprecipErr.tblList <- c(snoprecipErr.tblList, list(tbltmp))
                        }
                }
        }
        if (writeHtml) {
                for (i in 1:length(snoprecipErr.ggList)) {
                        cat(paste0("```{r snopreciperrmap_", i, ", fig.width = 9, fig.height = 6, out.width='900', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snoprecipErr.ggList[[{{i}}]]\n')
			cat(plottxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        # Table
                        cat(paste0("```{r snopreciperrtbl_", i, ", fig.width = 9, fig.height = 6, out.width='900', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                        #tbltxt <- knitr::knit_expand(text='pandoc.table(snoprecipErr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
			tbltxt <- knitr::knit_expand(text='print(xtable(snoprecipErr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        cat(tbltxt, file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
			cat('```\n', file=paste0(writePlotDir,"/plots_stats.Rmd"), append=TRUE)
                }
        }
}




# Output HTML
if (writeHtml) {
	if (accflowPlot | hydroPlot | flowlsmPlot) {
		knit2html(paste0(writePlotDir,"/plots_hydro.Rmd"), paste0(writePlotDir,"/plots_hydro.html"))
		file.remove("plots_hydro.md")
	}
	if (accprecipPlot) {	
		knit2html(paste0(writePlotDir,"/plots_climate.Rmd"), paste0(writePlotDir,"/plots_climate.html"))
		file.remove("plots_climate.md")
	}
	if (flowswePlot | swePlot) {	
		knit2html(paste0(writePlotDir,"/plots_snow.Rmd"), paste0(writePlotDir,"/plots_snow.html"))
		file.remove("plots_snow.md")
	}
	if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap) {	
		knit2html(paste0(writePlotDir,"/plots_stats.Rmd"), paste0(writePlotDir,"/plots_stats.html"))
		file.remove("plots_stats.md")
	}
	unlink("figure", recursive=TRUE)
}

## ------------------------------------------------------------------------
# Cleanup

proc.time()    

