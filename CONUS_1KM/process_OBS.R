library(doMC)
registerDoMC(16)

# Derived vars
obsAmfData$RgNet <- with(obsAmfData, Rg - RgOut)
obsAmfData$RglNet <- with(obsAmfData, Rgl - RglOut)
obsAmfData$EnResid <- with(obsAmfData, 
     (Rg - RgOut) + (Rgl - RglOut) - LE - H - FG)
obsAmfData$H2O_mmps <- with(obsAmfData, ifelse(TA>0, LE/2510/1000, LE/2844/1000))

# Date
obsAmfData$UTC_date <- as.Date(trunc(as.POSIXct(format(obsAmfData$POSIXct, tz="UTC"), tz="UTC"), "days"))

# Daily means
obsAmfData.d <- plyr::ddply(obsAmfData, plyr::.(site_id, UTC_date), 
                         plyr::summarise, 
                            TA_mean=mean(TA, na.rm=TRUE),
                            WS_mean=mean(WS, na.rm=TRUE),
                            NEE_mean=mean(NEE, na.rm=TRUE),
                            FC_mean=mean(FC, na.rm=TRUE),
                            H_mean=mean(H, na.rm=TRUE),
                            SH_mean=mean(SH, na.rm=TRUE),
                            LE_mean=mean(LE, na.rm=TRUE),
                            RH_mean=mean(RH, na.rm=TRUE),
                            PRESS_mean=mean(PRESS, na.rm=TRUE),
                            CO2_mean=mean(CO2, na.rm=TRUE),
                            VPD_mean=mean(VPD, na.rm=TRUE),
                            SWC1_mean=mean(SWC1, na.rm=TRUE),
                            SWC2_mean=mean(SWC2, na.rm=TRUE),
                            Rn_mean=mean(Rn, na.rm=TRUE),
                            Rg_mean=mean(Rg, na.rm=TRUE),
                            RgOut_mean=mean(RgOut, na.rm=TRUE),
                            Rgl_mean=mean(Rgl, na.rm=TRUE),
                            RglOut_mean=mean(RglOut, na.rm=TRUE),
                            H2O_mean=mean(H2O, na.rm=TRUE),
                            RE_mean=mean(RE, na.rm=TRUE),
                            GPP_mean=mean(GPP, na.rm=TRUE),
                            ZL_mean=mean(ZL, na.rm=TRUE),
                            H2Ommps_mean=mean(H2O_mmps, na.rm=TRUE),
                         .parallel=TRUE)
obsAmfData.d$POSIXct <- as.POSIXct(paste0(obsAmfData.d$UTC_date," 12:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")
obsAmfData.d$wy <- CalcWaterYear(obsAmfData.d$POSIXct)
obsAmfData.d$H2O_mmpd <- with(obsAmfData.d, H2Ommps_mean*86400)
