amf.URG.data <- GetAmeriflux(siteIDs=c("US-Vcm","US-Vcp","US-NR1","US-Cop"))
names(amf.URG.data)[names(amf.URG.data)=="site_id"]<-"site_nm"
amf.URG.data$H2O_mmps <- with(amf.URG.data, LE/2510/1000)
amf.URG.data$H2O_mmps2 <- with(amf.URG.data, ifelse(TA>0, LE/2510/1000, LE/2844/1000))
amf.URG.data$Date <- as.Date(trunc(as.POSIXct(format(amf.URG.data$POSIXct, tz="UTC"), tz="UTC"), "days"))
amf.URG.data$yr <- as.integer(format(amf.URG.data$POSIXct, "%Y"))
amf.URG.data$mo <- as.integer(format(amf.URG.data$POSIXct, "%m"))

# Daily means
amf.URG.data.dy <- plyr::ddply(amf.URG.data, plyr::.(site_nm, Date), 
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
                            H2Ommps2_mean=mean(H2O_mmps2, na.rm=TRUE),
                         .parallel=TRUE)
amf.URG.data.dy$POSIXct <- as.POSIXct(paste0(amf.URG.data.dy$Date+1," 00:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")
amf.URG.data.dy$wy <- CalcWaterYear(amf.URG.data.dy$POSIXct)
amf.URG.data.dy$H2O_mmpd <- with(amf.URG.data.dy, H2Ommps_mean*86400)
amf.URG.data.dy$H2O_mmpd2 <- with(amf.URG.data.dy, H2Ommps2_mean*86400)

# Month-year means
amf.URG.data.my <- plyr::ddply(amf.URG.data, plyr::.(site_nm, yr, mo),
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
                            H2Ommps2_mean=mean(H2O_mmps2, na.rm=TRUE),
                         .parallel=TRUE)
amf.URG.data.my$POSIXct <- as.POSIXct(paste0(amf.URG.data.my$yr,"-",amf.URG.data.my$mo,"-15"," 00:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")
amf.URG.data.my$wy <- CalcWaterYear(amf.URG.data.my$POSIXct)
amf.URG.data.my$H2O_mmpd <- with(amf.URG.data.my, H2Ommps_mean*86400)
amf.URG.data.my$H2O_mmpd2 <- with(amf.URG.data.my, H2Ommps2_mean*86400)
