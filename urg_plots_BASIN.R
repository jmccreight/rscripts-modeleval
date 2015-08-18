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

PlotAccPrecipFlow <- function(n, str1=modFrxstout_wy2015_NLDAS2dwnsc_fullrtng, 
                        lsm1=modLdasout_wy2015_NLDAS2dwnsc_fullrtng_BAS, 
                        obsstr=obsStr.dy,
                        str2=modFrxstout_wy2015_NLDAS2dwnsc_NSSL_fullrtng, 
                        lsm2=modLdasout_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS) {
  str1 <- subset(str1, str1$STAID==n)
  lsm1 <- subset(lsm1, lsm1$STAID==n)
  obsstr <- subset(obsstr, obsstr$Station==n)
  str2 <- subset(str2, str2$STAID==n)
  lsm2 <- subset(lsm2, lsm2$STAID==n)
  with(lsm1, plot(POSIXct, ACCPRCP, typ='l', col='darkmagenta'))
  with(str1, lines(POSIXct, ACCFLOW, col='deepskyblue'))
  with(lsm1, lines(POSIXct, ACCECAN+ACCEDIR+ACCETRAN, col='darkolivegreen3'))
  with(lsm2, lines(POSIXct, ACCPRCP, col='darkmagenta', lty=2))
  with(str2, lines(POSIXct, ACCFLOW, col='deepskyblue', lty=2))
  with(lsm2, lines(POSIXct, ACCECAN+ACCEDIR+ACCETRAN, col='darkolivegreen3', lty=2))
  with(obsstr, lines(POSIXct, cumqvol_mm, col='blue', lwd=2, lty=1))
}

PlotAccFlow <- function(n, modDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Accumulated Flow",
			lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
			modCol="ACCFLOW", obsCol="cumqvol_mm") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
	str1 <- modDfs
	strcnt <- 1
  } else if (is.list(modDfs)) {
	str1 <- modDfs[[1]]
	strcnt <- length(modDfs)
  } else {
	stop("modDfs must be a dataframe or a list of dataframes")
  }
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  str1 <- subset(str1, str1$STAID==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  obs <- subset(obs, obs$Station==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[nrow(str1),modCol]-str1[1,modCol], obs[nrow(obs),obsCol]-obs[1,obsCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
	for (stri in modDfs) {
		stri <- subset(stri, stri$STAID==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
		ymax <- max(ymax, stri[nrow(stri),modCol]-stri[1,modCol], na.rm=TRUE)
		}
	}
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  # Create plot
  plot(str1$POSIXct, str1[,modCol]-str1[1,modCol], typ='l', ylim=c(0, ymax), 
	xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Accumulated flow (mm)", cex.axis=1.2, cex.lab=1.2)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
		stri <- modDfs[[j]]
                stri <- subset(stri, stri$STAID==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
		lines(stri$POSIXct, stri[,modCol]-stri[1,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
		}
	}
  lines(obs$POSIXct, obs[,obsCol]-obs[1,obsCol], col='black', lwd=2, lty=1)
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt], labObs),
                lty=c(lnTyps[1:strcnt],1), lwd=c(lnWds[1:strcnt],2),
                col=c(lnCols[1:strcnt], 'black'), cex=1.2)
}


PlotFlowSwe <- function(n, modDfs, lsmDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Streamflow with Basin-Mean SWE",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="q_cms", lsmCol="SNEQV", obsCol="mean_qcms") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
        str1 <- modDfs
        strcnt <- 1
  } else if (is.list(modDfs)) {
        str1 <- modDfs[[1]]
        strcnt <- length(modDfs)
  } else {
        stop("modDfs must be a dataframe or a list of dataframes")
  }
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(lsmDfs)) {
        lsm1 <- lsmDfs
        lsmcnt <- 1
  } else if (is.list(lsmDfs)) {
        lsm1 <- lsmDfs[[1]]
        lsmcnt <- length(lsmDfs)
  } else {
        stop("lsmDfs must be a dataframe or a list of dataframes")
  }
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  str1 <- subset(str1, str1$STAID==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  lsm1 <- subset(lsm1, lsm1$statArg==n)
  obs <- subset(obs, obs$Station==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[,modCol], obs[,obsCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (stri in modDfs) {
                stri <- subset(stri, stri$STAID==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                ymax <- max(ymax, stri[,modCol], na.rm=TRUE)
                }
        }
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  # Create plot
  plot(str1$POSIXct, str1[,modCol], typ='l', ylim=c(0, ymax),
        xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Streamflow (m3/s) or SWE (cm)", cex.axis=1.2, cex.lab=1.2)
  lines(lsm1$POSIXct, lsm1[,lsmCol]/10, col=lnCols[1], lty=2, lwd=lnWds[1])
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
		lsmi <- lsmDfs[[j]]
                stri <- subset(stri, stri$STAID==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
		lsmi <- subset(lsmi, lsmi$statArg==n)
                lines(stri$POSIXct, stri[,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                lines(lsmi$POSIXct, lsmi[,lsmCol]/10, col=lnCols[j], lty=2, lwd=lnWds[j])
		}
        }
  lines(obs$POSIXct, obs[,obsCol], col='black', lwd=2, lty=1)
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt], labObs, "", "Streamflow (m3/s)", "SWE (cm)"),
                lty=c(lnTyps[1:strcnt],1,1,1,2), lwd=c(lnWds[1:strcnt],2,1,2,2),
                col=c(lnCols[1:strcnt], 'black','white','grey40','grey40'), cex=1.2,
		bg="white")
}




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
        png(paste0(outPath, "/strflow_swe_NLDAS_NSSL_adj_PRES", "_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlowSwe(n, modDfs=list(modFrxstout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng,
                        modFrxstout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng),
                        lsmDfs=list(modLdasout_wy2015_NLDAS2dwnsc_snowmod_mikerec_fullrtng_BAS,
                        modLdasout_wy2015_NLDAS2dwnsc_NSSL_snowmod_mikerec_fullrtng_BAS),
                        obs=obsStr.dy,
                        labMods=c("NLDAS-2", "NLDAS-2 + NSSL Precipitation"),
			lnCols=c("dodgerblue", "darkorange1"),
                        lnWds=c(3,3),
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n, ", WY2014"),
                        stdate=NULL, enddate=enddate, obsCol="mean_qcms_adj")
        dev.off()
}


### EXIT

quit("no")

