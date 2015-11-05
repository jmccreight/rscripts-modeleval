###################################################
##                GENERATE PLOTS                 ##
###################################################


################## General Setup ##################

dir.create(writePlotDir, showWarnings = FALSE)
source("PlotBasins.R")
source("PlotSnotel.R")
source("PlotMaps.R")

lineColors <- c("dodgerblue", "darkorange1", "olivedrab", "chocolate", "darkmagenta")
lineTyp <- 1
lineWd <- 3

if (writeHtml) {
	library(knitr)
	writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', con="example.Rmd")
	cat('# MODEL OUTPUT\n', file="example.Rmd", append=TRUE)
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
                        labTitle=paste0("Accumulated Flow: ", n, " (", obsStrMeta$statname[obsStrMeta$site_no==n], ")"), obsCol="cumqvol_mm")
        dev.off()
}
if (writeHtml) {
	cat('## Accumulated Flow Plots\n', file="example.Rmd", append=TRUE)
	for (n in names(gage2basinList)) {
		cat(paste0("```{r accflow_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), file="example.Rmd", append=TRUE)
		plottxt <- knitr::knit_expand(text='PlotAccFlow("{{n}}", modDfs=accflowList,
                        obs=obsStrData.dy,
                        stdate=accflowStartDate,
                        enddate=accflowEndDate,
                        labMods=accflowTags,
                        labObs="Observed",
                        lnCols=accflowColors,
                        lnTyps=accflowTypes, lnWds=accflowWidths,
                        labTitle=paste0("Accumulated Flow: ", "{{n}}", " (", obsStrMeta$statname[obsStrMeta$site_no=="{{n}}"], ")"), obsCol="cumqvol_mm")\n')
		cat(plottxt, file="example.Rmd", append=TRUE)
		cat('```\n', file="example.Rmd", append=TRUE)
	}
}
}

# Hydrographs
if (hydroPlot) {
message("Generating hydrograph plots...")
# Setup
hydroList <- list()
if (is.null(hydroTags)) hydroTags <- unique(modFrxstout$tag)
for (i in 1:length(hydroTags)) {
        hydroList[[i]] <- subset(modFrxstout, modFrxstout$tag==hydroTags[i])
}
hydroColors <- lineColors[1:length(hydroList)]
hydroTypes <- rep(lineTyp, length(hydroList))
hydroWidths <- rep(lineWd, length(hydroList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/hydrogr_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlow(n, modDfs=hydroList,
                        obs=obsStrData.dy,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=paste0("Streamflow: ", n, " (", obsStrMeta$statname[obsStrMeta$site_no==n], ")"),
                        stdate=hydroStartDate, enddate=hydroEndDate, obsCol="mean_qcms")
        dev.off()
}
if (writeHtml) {
        cat('## Hydrographs\n', file="example.Rmd", append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r hydro_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), file="example.Rmd", append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotFlow("{{n}}", modDfs=hydroList,
                        obs=obsStrData.dy,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=paste0("Streamflow: ", "{{n}}", " (", obsStrMeta$statname[obsStrMeta$site_no=="{{n}}"], ")"),
                        stdate=hydroStartDate, enddate=hydroEndDate, obsCol="mean_qcms")\n')
                cat(plottxt, file="example.Rmd", append=TRUE)
                cat('```\n', file="example.Rmd", append=TRUE)
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
                        labTitle=paste0("Accumulated Precipitation: ", n))
        dev.off()
}
if (writeHtml) {
        cat('## Accumulated Precip Plots\n', file="example.Rmd", append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r accprecip_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), file="example.Rmd", append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotAccPrecip("{{n}}", modDfs=accprecipList,
                        stdate=accprecipStartDate,
                        enddate=accprecipEndDate,
                        labMods=accprecipTags,
                        lnCols=accprecipColors,
                        lnTyps=accprecipTypes, lnWds=accprecipWidths,
                        labTitle=paste0("Accumulated Precipitation: ", "{{n}}"))\n')
                cat(plottxt, file="example.Rmd", append=TRUE)
                cat('```\n', file="example.Rmd", append=TRUE)
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
        cat('## Streamflow & Basin-mean SWE Plots\n', file="example.Rmd", append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r flowswe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), file="example.Rmd", append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotFlowSwe("{{n}}", modDfs=flowsweStrList,
                        lsmDfs=flowsweLsmList,
                        obs=obsStrData.dy,
                        labMods=flowsweTags,
                        lnCols=flowsweColors,
                        lnWds=flowsweWidths,
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", "{{n}}"),
                        stdate=flowsweStartDate, enddate=flowsweEndDate)\n')
                cat(plottxt, file="example.Rmd", append=TRUE)
                cat('```\n', file="example.Rmd", append=TRUE)
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
        cat('## Station SWE Plots\n', file="example.Rmd", append=TRUE)
        for (n in sites) {
                cat(paste0("```{r swe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), file="example.Rmd", append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotSwe("{{n}}", modDfs=sweList,
                	obs=obsSnoData, obsmeta=obsSnoMeta,
                	labMods=sweTags,
                	lnCols=sweColors,
                	lnWds=sweWidths,
                	precCol.obs="CumPrec_mm", precCol.mod="ACCPRCP",
                	sweCol.obs="SWE_mm", sweCol.mod="SNEQV", fact=1, snowh=FALSE,
                	labTitle="Accumulated Precipitation and SWE",
                	stdate=sweStartDate, enddate=sweEndDate)\n')
                	cat(plottxt, file="example.Rmd", append=TRUE)
                	cat('```\n', file="example.Rmd", append=TRUE)
        }
}
}

# Initialize for maps
if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap ) {
	library(ggplot2)
	library(ggmap)
	# Setup date ranges
	stdate_stats_PRINT <- ifelse(is.null(stdate_stats), min(modFrxstout$POSIXct), stdate_stats)
	enddate_stats_PRINT <- ifelse(is.null(enddate_stats), max(modFrxstout$POSIXct), enddate_stats)
	stdate_stats_sub_PRINT <- ifelse(is.null(stdate_stats_sub), min(modFrxstout$POSIXct), stdate_stats_sub)
	enddate_stats_sub_PRINT <- ifelse(is.null(enddate_stats_sub), max(modFrxstout$POSIXct), enddate_stats_sub)
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
        cat('## Streamflow Bias Maps\n', file="example.Rmd", append=TRUE)
	ggList <- list()
}
counter <- 1
for (i in strBiasTags) {
        for (j in strBiasSeas) {
                gg <- PlotMapErrors(geoMap, stats_str,
                        statsTag=i, statsVar=NULL, statsSeas=j,
                        plotTitle="Modeled Streamflow Bias at CODWR Gages",
			plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        sizeVar="t_bias", colorVar="t_bias",
                        sizeLab="Bias (%)", colorLab="Bias (%)")
                ggplot2::ggsave(filename=paste0(writePlotDir, "/str_bias_map_", i, "_", j, ".png"),
                        plot=gg, units="in", width=8, height=6, dpi=100)
		if (writeHtml) {
			ggList <- list(ggList, gg)
		}
	}
}
counter <- 1
if (writeHtml) {
	while (counter <= length(ggList)) {
			if (!(counter %% 2 == 0)) {
			cat(paste0("```{r strbiasmap_", i, "_", j, ", fig.width = 16, fig.height = 6, out.width='1600', out.height='600', echo=FALSE, out.extra='style=\"float:left\"'}\n"), 
				file="example.Rmd", append=TRUE)
			if (counter %% 2 == 0) {

#			plottxt <- knitr::knit_expand(text='PlotMapErrors(geoMap, stats_str,
#                        	statsTag="{{i}}", statsVar=NULL, statsSeas="{{j}}",
#                        	plotTitle="Modeled Streamflow Bias at CODWR Gages",
#                        	plotSubTitle=paste0("{{i}}", ", ", statsDateList[["{{j}}"]]),
#                        	sizeVar="t_bias", colorVar="t_bias",
#                        	sizeLab="Bias (%)", colorLab="Bias (%)")\n')
			cat(plottxt, file="example.Rmd", append=TRUE)
			if (counter %% 2 == 0) cat('```\n', file="example.Rmd", append=TRUE)
			counter <- counter + 1
			}
                }
        }
}

# STRFLOW Correlation Maps
if (strCorrMap) {
message("Generating STRFLOW Corr error map...")
# Setup
if (is.null(strCorrTags)) strCorrTags <- unique(stats_str$tag)
if (is.null(strCorrSeas)) strCorrSeas <- unique(stats_str$seas)
for (i in strCorrTags) {
        for (j in strCorrSeas) {
                gg <- PlotMapErrors(geoMap, stats_str,
                        statsTag=i, statsVar=NULL, statsSeas=j,
                        plotTitle="Modeled Streamflow Correlation at CODWR Gages",
			plotSubTitle=paste0(i, ", ", statsDateList[[j]]),
                        sizeVar="dy_cor", colorVar="dy_cor",
                        sizeLab="Correlation", colorLab="Correlation",
			colorLow="orange", colorMid="yellow", colorHigh="cyan4")
                ggplot2::ggsave(filename=paste0(writePlotDir, "/str_corr_map_", i, "_", j, ".png"),
                        plot=gg, units="in", width=8, height=6, dpi=100)
                }
        }
}

# SNOTEL SWE Maps
if (snosweErrMap) {
message("Generating SNOTEL SWE error map...")
# Setup
if (is.null(snosweErrTags)) snosweErrTags <- unique(stats_ldasout_sno$tag)
if (is.null(snosweErrSeas)) snosweErrSeas <- unique(stats_ldasout_sno$seas)
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
		}
	}
}

# SNOTEL Precip Maps
if (snoprecipErrMap) {
message("Generating SNOTEL precip error map...")
# Setup
if (is.null(snoprecipErrTags)) snoprecipErrTags <- unique(stats_ldasout_sno$tag)
if (is.null(snoprecipErrSeas)) snoprecipErrSeas <- unique(stats_ldasout_sno$seas)
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
                }
        }
}

if (writeHtml) {
	knit2html("example.Rmd", "example.html")
}

## ------------------------------------------------------------------------
# Cleanup

proc.time()    

