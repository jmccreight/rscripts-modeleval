###################################################################################################
# Setup
# 
# Load the rwrfhydro package. 
## ------------------------------------------------------------------------
library("rwrfhydro")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_masks_NEW.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/STRFLOW/strflow_URG.Rdata")

# Where the processed data lives
rimgPath <- 'urg_wy2015_ALL_PROCESSED.Rdata'

# Where to write files
outPath <- 'plots'

# Range dates to restrict plots
stdate <- NULL
#enddate <- NULL
enddate <- as.POSIXct("2015-07-14 00:00", format="%Y-%m-%d %H:%M", tz="UTC")


###################################################################################################
# Run

load(rimgPath)
dir.create(outPath, showWarnings = FALSE)

## ------------------------------------------------------------------------
# Plot Functions

source("PlotBasins.R")


### PLOTS

# Cumulative flow plots comparing all
for (n in names(gage2basinList)) {
	png(paste0(outPath, "/accstrflow_ALL", "_", n, ".png"), width=800, height=500)
	PlotAccFlow(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_snowmod_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_NSSL_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_nlcd11_fullrtng),
			obs=obsStr.dy,
			stdate=as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
			enddate=enddate, 
			labMods=c("NLDAS-2, Base", "NLDAS-2, Snow Modifications", "NLDAS-2, Snow Mods + Mike Recs", 
				  "NSSL, Base", "NSSL, Snow Modifications", "NSSL, Snow Mods + Mike Recs", "NSSL, Snow Mods + Mike Recs + NLCD11"),
			lnCols=c("darkmagenta", "darkmagenta", "darkmagenta", "cyan4", "cyan4", "cyan4", "chocolate"),
			lnTyps=c(1,2,3,1,2,3,3), lnWds=c(1.5,1.5,2,1.5,1.5,2,2),
			labTitle=paste0("Accumulated Flow: ", n, ", April-July 2015"), obsCol="cumqvol_mm")
        dev.off()
}

# Streamflow + SWE plots comparing NLDAS & NSSL
for (n in names(gage2basinList)) {
	png(paste0(outPath, "/strflow_swe_NLDAS_NSSL", "_", n, ".png"), width=800, height=500)
	PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("Model1: NLDAS-2, Base", "Model2: NSSL, Base"),
			lnCols=c("chocolate", "dodgerblue"),
			lnWds=c(2,2),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2014"),
                        stdate=NULL, enddate=enddate)
	dev.off()
}


# Streamflow + SWE plots comparing NLDAS w/ & w/o snow mods & Mike recs
for (n in names(gage2basinList)) {
        png(paste0(outPath, "/strflow_swe_NLDAS_modcomp", "_", n, ".png"), width=800, height=500)
        PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_snowmod_fullrtng,
			modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_snowmod_fullrtng_BAS,
			modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("Model1: NLDAS-2, Base", "Model2: NLDAS-2, Snow Modifications",
				"Model3: NLDAS-2, Snow Mods with Mike's Recs"),
			lnCols=c("chocolate", "dodgerblue", "chartreuse3"),
                        lnWds=c(2,2,2),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2014"),
                        stdate=stdate, enddate=enddate)
        dev.off()
}

# Streamflow + SWE plots comparing NSSL w/ & w/o snow mods
for (n in names(gage2basinList)) {
        png(paste0(outPath, "/strflow_swe_NSSL_modcomp", "_", n, ".png"), width=800, height=500)
        PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_NSSL_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_nlcd11_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_nlcd11_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("Model1: NSSL, Base", "Model2: NSSL, Snow Modifications",
                                        "Model3: NSSL, Snow Modifications, Mike Recs",
                                        "Model4: NSSL, Snow Modifications, Mike Recs, NLCD11"),
                        lnCols=c("chocolate", "dodgerblue", "chartreuse3", "purple"),
                        lnWds=c(2,2,2,2),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2014"),
                        stdate=stdate, enddate=enddate)
        dev.off()
}


# Streamflow stats

# Bias
stats_str_all$runsort <- factor(stats_str_all$run, as.character(stats_str_all$run))
gg <- ggplot(data=stats_str_all, aes(x=STAID, y=t_bias)) + geom_violin() +
	geom_point(size=3, aes(color=factor(runsort))) +
	scale_color_manual(values=c("darkblue","brown4","dodgerblue","coral2","cadetblue2","sandybrown","olivedrab","pink"), name="Model Run") + 
	ylim(-50,400) +
	labs(x="STATION", y="STREAMFLOW BIAS (%)") +
	theme(axis.text.x=element_text(size=8, vjust=0.5, angle=50)) +
	ggtitle("Streamflow Bias (Apr 1 through ~Jun 30, 2015)") +
	geom_hline(yintercept=0)
ggsave(filename=paste0(outPath, "/stats_str_bias.png"), plot=gg,
       units="in", width=12, height=6, dpi=100)

# Correlation
gg <- ggplot(data=stats_str_all, aes(x=STAID, y=dy_cor)) + geom_violin() +
        geom_point(size=3, aes(color=factor(runsort))) +
        scale_color_manual(values=c("darkblue","brown4","dodgerblue","coral2","cadetblue2","sandybrown","olivedrab","pink"), name="Model Run") +
        #ylim(-50,400) +
        labs(x="STATION", y="STREAMFLOW DAILY CORRELATION (%)") +
        theme(axis.text.x=element_text(size=9, vjust=0.5, angle=50)) +
        ggtitle("Streamflow Daily Correlation (Apr 1 through ~Jun 30, 2015)") +
	geom_hline(yintercept=0)
ggsave(filename=paste0(outPath, "/stats_str_corr.png"), plot=gg,
       units="in", width=12, height=6, dpi=100)


# PRESENTATION PLOTS

# Cumulative flow plots
for (n in names(gage2basinList)) {
        png(paste0(outPath, "/accstrflow_adj_PRES", "_", n, ".png"), width=2100, height=1350, res=225)
       	if (n %in% c("CONMOGCO")) {
		labObs <- "Observed (Naturalized)"
	} else {
		labObs <- "Observed"} 
       	PlotAccFlow(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng,
                                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng),
                        obs=obsStr.dy,
                        stdate=as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                        enddate=enddate,
                        labMods=c("NLDAS-2", "NLDAS-2 + NSSL Precipitation"),
			labObs=labObs,
                        lnCols=c("dodgerblue", "darkorange1"),
                        lnTyps=c(1,1), lnWds=c(3,3),
                        labTitle=paste0("Accumulated Flow: ", n, ", April-July 2015"), obsCol="cumqvol_mm_adj")
        dev.off()
}

# Streamflow + SWE plots
for (n in names(gage2basinList)) {
        if (n %in% c("CONMOGCO", "CONPLACO", "RIOWAGCO", "RIODELCO")) {
                labObs <- "Observed (Naturalized)"
        } else {
                labObs <- "Observed"}
        png(paste0(outPath, "/strflow_swe_NLDAS_NSSL_adj_PRES", "_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("NLDAS-2", "NLDAS-2 + NSSL Precipitation"),
                        labObs=labObs,
                        lnCols=c("dodgerblue", "darkorange1"),
                        lnWds=c(3,3),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2014"),
                        stdate=NULL, enddate=enddate, obsCol="mean_qcms_adj")
        dev.off()
}


for (n in names(gage2basinList)) {
        if (n %in% c("CONMOGCO", "CONPLACO", "RIOWAGCO", "RIODELCO")) {
                labObs <- "Observed (Naturalized)"
        } else {
                labObs <- "Observed"}
        png(paste0(outPath, "/strflow_swe_NSSL_adj_PRES", "_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("NSSL Precipitation", "NSSL w/Resistance Mods"),
                        labObs=labObs,
                        lnCols=c("darkorange1", "olivedrab"),
                        lnWds=c(3,3),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2015"),
                        stdate=NULL, enddate=enddate, obsCol="mean_qcms_adj")
        dev.off()
}



################################

# Acc Streamflow - 150908
for (n in names(gage2basinList)) {
        png(paste0(outPath, "/accstrflow_NLDAS_NSSL_", n, ".png"), width=2100, height=1350, res=225)
        if (n %in% c("CONMOGCO", "CONPLACO", "RIOWAGCO", "RIODELCO")) {
                labObs <- "Observed (Naturalized)"
        } else {
                labObs <- "Observed"}
        PlotAccFlow(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_snowresist50_fullrtng,
					modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng),
                        obs=obsStr.dy,
                        stdate=as.POSIXct("2015-04-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                        enddate=enddate,
                        labMods=c("NLDAS-2", "NLDAS-2 + NSSL Precipitation"),
                        labObs=labObs,
                        lnCols=c("dodgerblue", "darkorange1"),
                        lnTyps=c(1,1), lnWds=c(3,3),
                        labTitle=paste0("Accumulated Flow: ", n, ", April-August 2015"), obsCol="cumqvol_mm_adj")
        dev.off()
}

# Streamflow - 150908
for (n in names(gage2basinList)) {
        if (n %in% c("CONMOGCO", "CONPLACO", "RIOWAGCO", "RIODELCO")) {
                labObs <- "Observed (Naturalized)"
        } else {
                labObs <- "Observed"}
        png(paste0(outPath, "/strflow_NLDAS_NSSL_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlow(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_snowresist50_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng),
                        obs=obsStr.dy,
                        labMods=c("NLDAS-2", "NLDAS-2 + NSSL Precipitation"),
                        labObs=labObs,
                        lnCols=c("dodgerblue", "darkorange1"),
                        lnWds=c(3,3),
                        labTitle=paste0("Streamflow: ", n, ", WY2015"),
                        stdate=NULL, enddate=enddate, obsCol="mean_qcms_adj")
        dev.off()
}

# Acc Precip - 150908
for (n in names(gage2basinList)) {
        png(paste0(outPath, "/accprecip_NLDAS_NSSL_", n, ".png"), width=2100, height=1350, res=225)
        PlotAccPrecip(n, modDfs=list(modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_snowresist50_fullrtng_BAS,
                                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_snowresist50_fullrtng_BAS),
                        stdate=as.POSIXct("2014-10-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC"),
                        enddate=enddate,
                        labMods=c("NLDAS-2 Precipitation", "NSSL Precipitation"),
                        lnCols=c("dodgerblue", "darkorange1"),
                        lnTyps=c(1,1), lnWds=c(3,3),
                        labTitle=paste0("Accumulated Precipitation: ", n, ", WY2015"))
        dev.off()
}





### EXIT

quit("no")

