###################################################################################################
# Setup
# 
# Load the rwrfhydro package. 
## ------------------------------------------------------------------------
library("rwrfhydro")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/ANALYSIS/urg_masks_NEW.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/SNOTEL/snotel_URG.Rdata")
load("/glade/p/ral/RHAP/adugger/Upper_RioGrande/OBS/MET/met_URG.Rdata")

# Where the processed data lives
rimgPath <- 'urg_wy2015_SNOCLIM_PROCESSED.Rdata'

# Where to write files
outPath <- 'plots'

# Range dates to restrict plots
stdate <- NULL
enddate <- as.POSIXct("2015-07-14 00:00", format="%Y-%m-%d %H:%M", tz="UTC")


###################################################################################################
# Run

load(rimgPath)
dir.create(outPath, showWarnings = FALSE)

source("PlotSWE.R")


## ------------------------------------------------------------------------
### PLOTS

sites <- unique(sno.URG.sites$site_id)
sites <- c(sites, unique(met.URG.sites$site_id))

# Accumulated precip and SWE plots: NLDAS & Snow Mods & Mike Recs
for (n in 1:length(sites)) {
  print(sites[n])
  png(paste0(outPath, "/cumprec_swe_nldas_snowmods_", sites[n], ".png"), width=800, height=500)
  PlotSwe(sites[n], modDfs=list(modLdasout_wy2015_NLDAS2dwnsc_fullrtng_SNO, 
				modLdasout_wy2015_NLDAS2dwnsc_snowmod_fullrtng_SNO,
				modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_SNO),
		obs=sno.URG.data, obsmeta=sno.URG.sites, 
          	labMods=c("Model1: NLDAS Precip, Base", 
			"Model2: NLDAS Precip, NoahMP Snow Mods",
			"Model3: NLDAS Precip, NoahMP Snow Mods, Mike Recs"),
		lnCols=c('chocolate', 'dodgerblue', 'chartreuse3'),
		lnWds=c(2,2,2,2))
  dev.off()
}

# Accumulated precip and SWE plots: NSSL & Snow Mods
for (n in 1:length(sites)) {
  print(sites[n])
  png(paste0(outPath, "/cumprec_swe_nssl_snowmods_", sites[n], ".png"), width=800, height=500)
  PlotSwe(sites[n], modDfs=list(modLdasout_wy2015_NLDAS2dwnsc_NSSL_fullrtng_SNO, 
				modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_fullrtng_SNO,
 				modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_SNO,
                                modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_nlcd11_fullrtng_SNO),
                obs=sno.URG.data, obsmeta=sno.URG.sites,
                labMods=c("Model1: NSSL Precip, Base",
                        "Model2: NSSL Precip, NoahMP Snow Mods",
                        "Model3: NSSL Precip, NoahMP Snow Mods, Mike Recs",
                        "Model4: NSSL Precip, NoahMP Snow Mods, Mike Recs, NLCD 2011"),
                lnCols=c('chocolate', 'dodgerblue', 'chartreuse3', 'purple'),
		lnWds=c(2,2,2,2,2))
  dev.off()
}

# Boxplot
stats_sno_all$runsort <- factor(stats_sno_all$run, as.character(stats_sno_all$run))
gg <- ggplot(data=stats_sno_all, aes(x=runsort, y=t_bias, color=runsort)) +
        geom_boxplot() +
	theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
	ylab("SWE Bias (%)") + xlab("Model Run") +
        scale_color_manual(values=c("darkblue","brown4","dodgerblue","coral2","cadetblue2","sandybrown","olivedrab","pink"), name="Model Run") +
        ggtitle("SWE Bias (Oct 1, 2014 through ~Jun 30, 2015)") +
        geom_hline(yintercept=0)
ggsave(filename=paste0(outPath, "/stats_sno_bias.png"), plot=gg,
       units="in", width=12, height=6, dpi=100)

# PRESENTATION

# Accumulated precip and SWE plots
sites <- unique(sno.URG.sites$site_id)
for (n in 1:length(sites)) {
  print(sites[n])
  png(paste0(outPath, "/cumprec_swe_pres_", sites[n], ".png"), width=2100, height=1350, res=225)
  PlotSwe(sites[n], modDfs=list(modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_SNO,
				modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_SNO),
                obs=sno.URG.data, obsmeta=sno.URG.sites,
		labMods=c("Model1: NLDAS-2", "Model2: NLDAS-2 + NSSL Precipitation"),
                lnCols=c('dodgerblue', 'darkorange1'),
                lnWds=c(3,3))
  dev.off()
}

sites <- unique(met.URG.sites$site_id)
for (n in 1:length(sites)) {
  print(sites[n])
  png(paste0(outPath, "/cumprec_swe_pres_", sites[n], ".png"), width=2100, height=1350, res=225)
  PlotSwe(sites[n], modDfs=list(modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_SNO,
                                modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_SNO),
                obs=met.URG.data.dy, obsmeta=met.URG.sites,
                labMods=c("Model1: NLDAS-2", "Model2: NSSL Precipitation"),
                lnCols=c('dodgerblue', 'darkorange1'),
                lnWds=c(3,3),
		precCol.obs="PrecAcc_max", precCol.mod="ACCPRCP", 
		sweCol.obs="SnoDep_mean", sweCol.mod="SNOWH", fact=1000, snowh=TRUE,
		labTitle="Accumulated Precipitation and Snow Depth",
		#stdate_prec=as.POSIXct("2014-11-21", format="%Y-%m-%d", tz="UTC"))
		stdate_prec=min(subset(met.URG.data.dy$POSIXct, met.URG.data.dy$site_id==sites[n])))
  dev.off()
}



### EXIT

quit("no")

