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

PlotAccPrecip <- function(n, modDfs,
                        labMods=NULL,
                        labTitle="Accumulated Precipitation",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="ACCPRCP") {
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
  str1 <- subset(str1, str1$statArg==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[nrow(str1),modCol]-str1[1,modCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (stri in modDfs) {
                stri <- subset(stri, stri$statArg==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
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
        xlab="", ylab="Accumulated precipitation (mm)", cex.axis=1.2, cex.lab=1.2)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
                stri <- subset(stri, stri$statArg==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                lines(stri$POSIXct, stri[,modCol]-stri[1,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                }
        }
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt]),
                lty=c(lnTyps[1:strcnt]), lwd=c(lnWds[1:strcnt]),
                col=c(lnCols[1:strcnt]), cex=1.2)
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


PlotFlow <- function(n, modDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Streamflow with Basin-Mean SWE",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="q_cms", obsCol="mean_qcms") {
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
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
                stri <- subset(stri, stri$STAID==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                lines(stri$POSIXct, stri[,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                }
        }
  lines(obs$POSIXct, obs[,obsCol], col='black', lwd=2, lty=1)
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt], labObs),
                lty=c(lnTyps[1:strcnt],1), lwd=c(lnWds[1:strcnt],2),
                col=c(lnCols[1:strcnt], 'black'), cex=1.2,
                bg="white")
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

