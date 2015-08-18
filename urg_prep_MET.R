library(rwrfhydro)
load("/glade/p/ral/RHAP/alyssah/DG_stations/UpperRioGrande_stations_10012014-07012015.RData")
load("../OBS/MET/met_URG.Rdata")
met.URG.data <- data.frame()
for (i in 1:6) {
    tmp <- get(paste0("URG",i))[[1]]
    tmp <- tmp[order(tmp$POSIXct),]
    if ("Total_Accumulated_Precipitation" %in% names(tmp)) {
        tmp$Prec_mm <- 0
        tmp$Prec_mm[2:nrow(tmp)] <- diff(tmp$Total_Accumulated_Precipitation)
        }
    tmp$site_id <- i
    tmp$ptrun_id <- met.URG.sites$ptrun_id[met.URG.sites$site_id==i]
    met.URG.data <- plyr::rbind.fill(met.URG.data, tmp)
    }
# Calculate hourly aggregations
met.URG.data$Hour <- as.character(trunc(as.POSIXct(format(met.URG.data$POSIXct, tz="UTC"), tz="UTC"), "hours"))
met.URG.data.hr <- plyr::ddply(met.URG.data, plyr::.(site_id, Hour),
                         plyr::summarise,
                         RH_mean=mean(relative_humidity, na.rm=TRUE),
                         T_mean=mean(temperature, na.rm=TRUE),
                         Wind_mean=mean(wind, na.rm=TRUE),
                         WindDir_mean=mean(wdirection, na.rm=TRUE),
                         SoilTemp1_mean=mean(Soil_Temp1, na.rm=TRUE),
                         SM1_mean=mean(Soil_moist1, na.rm=TRUE),
                         SoilTemp2_mean=mean(Soil_Temp2, na.rm=TRUE),
                         SM2_mean=mean(Soil_moist2, na.rm=TRUE),
                         SnoDep_mean=mean(snow_depth, na.rm=TRUE),
                         SWRad_mean=mean(shortwave_radiation, na.rm=TRUE),
                         LeafWet_mean=mean(Leaf_Wetness, na.rm=TRUE),
                         PrecAcc_max=max(Total_Accumulated_Precipitation, na.rm=TRUE),
                         PrecTot=sum(Prec_mm, na.rm=TRUE),
                         PrecInt_mean=mean(Precipitation_Intensity, na.rm=TRUE),
                         SurfPress_mean=mean(Surface_Pressure, na.rm=TRUE),
                         .parallel=FALSE)
# Create daily aggregations
met.URG.data$Date <- as.Date(trunc(as.POSIXct(format(met.URG.data$POSIXct, tz="UTC"), tz="UTC"), "days"))
met.URG.data.dy <- plyr::ddply(met.URG.data, plyr::.(site_id, Date),
                         plyr::summarise,
                         RH_mean=mean(relative_humidity, na.rm=TRUE),
                         T_mean=mean(temperature, na.rm=TRUE),
                         T_min=min(temperature, na.rm=TRUE),
                         T_max=max(temperature, na.rm=TRUE),
                         Wind_mean=mean(wind, na.rm=TRUE),
                         WindDir_mean=mean(wdirection, na.rm=TRUE),
                         SoilTemp1_mean=mean(Soil_Temp1, na.rm=TRUE),
                         SoilTemp1_min=min(Soil_Temp1, na.rm=TRUE),
                         SoilTemp1_max=max(Soil_Temp1, na.rm=TRUE),
                         SM1_mean=mean(Soil_moist1, na.rm=TRUE),
                         SoilTemp2_mean=mean(Soil_Temp2, na.rm=TRUE),
                         SoilTemp2_min=min(Soil_Temp2, na.rm=TRUE),
                         SoilTemp2_max=max(Soil_Temp2, na.rm=TRUE),
                         SM2_mean=mean(Soil_moist2, na.rm=TRUE),
                         SnoDep_mean=mean(snow_depth, na.rm=TRUE),
                         SWRad_mean=mean(shortwave_radiation, na.rm=TRUE),
                         LeafWet_mean=mean(Leaf_Wetness, na.rm=TRUE),
                         PrecAcc_max=max(Total_Accumulated_Precipitation, na.rm=TRUE),
                         PrecTot=sum(Prec_mm, na.rm=TRUE),
                         PrecInt_mean=mean(Precipitation_Intensity, na.rm=TRUE),
                         SurfPress_mean=mean(Surface_Pressure, na.rm=TRUE),
                         .parallel=FALSE)
# Do unit conversions to match LDASIN
met.URG.data$Temp_K <- met.URG.data$temperature + 273.15
met.URG.data$SurfPress_Pa <- met.URG.data$Surface_Pressure * 1000
met.URG.data$SnoDep_m <- met.URG.data$snow_depth / 100
met.URG.data$SoilTemp1_K <- met.URG.data$Soil_Temp1 + 273.15
met.URG.data$SoilTemp2_K <- met.URG.data$Soil_Temp2 + 273.15
# Similar for hourly
met.URG.data.hr$Tmean_K <- met.URG.data.hr$T_mean + 273.15
met.URG.data.hr$SurfPressmean_Pa <- met.URG.data.hr$SurfPress_mean * 1000
met.URG.data.hr$SnoDepmean_m <- met.URG.data.hr$SnoDep_mean / 100
met.URG.data.hr$SoilT1mean_K <- met.URG.data.hr$SoilTemp1_mean + 273.15
met.URG.data.hr$SoilT2mean_K <- met.URG.data.hr$SoilTemp2_mean + 273.15
# Similar for daily
met.URG.data.dy$Tmean_K <- met.URG.data.dy$T_mean + 273.15
met.URG.data.dy$Tmin_K <- met.URG.data.dy$T_min + 273.15
met.URG.data.dy$Tmax_K <- met.URG.data.dy$T_max + 273.15
met.URG.data.dy$SurfPressmean_Pa <- met.URG.data.dy$SurfPress_mean * 1000
met.URG.data.dy$SnoDepmean_m <- met.URG.data.dy$SnoDep_mean / 100
met.URG.data.dy$SoilT1mean_K <- met.URG.data.dy$SoilTemp1_mean + 273.15
met.URG.data.dy$SoilT1min_K <- met.URG.data.dy$SoilTemp1_min + 273.15
met.URG.data.dy$SoilT1max_K <- met.URG.data.dy$SoilTemp1_max + 273.15
met.URG.data.dy$SoilT2mean_K <- met.URG.data.dy$SoilTemp2_mean + 273.15
met.URG.data.dy$SoilT2min_K <- met.URG.data.dy$SoilTemp2_min + 273.15
met.URG.data.dy$SoilT2max_K <- met.URG.data.dy$SoilTemp2_max + 273.15

# Make NAs consistent
met.URG.data.hr[met.URG.data.hr==(-Inf)] <- NA
met.URG.data.hr[met.URG.data.hr==(Inf)] <- NA
met.URG.data.dy[met.URG.data.dy==(-Inf)] <- NA
met.URG.data.dy[met.URG.data.dy==(Inf)] <- NA

# Add POSIXct for ease of plotting. Set to 00:00 on the next day to match SNOTEL
# and daily LDASIN aggregations.
met.URG.data.hr$POSIXct <- as.POSIXct(met.URG.data.hr$Hour, format="%Y-%m-%d %H:%M", tz="UTC")
met.URG.data.dy$POSIXct <- as.POSIXct(paste0(met.URG.data.dy$Date+1," 00:00",
                                       format="%Y-%m-%d %H:%M", tz="UTC"), tz="UTC")


save(met.URG.sites, met.URG.data, met.URG.data.hr, met.URG.data.dy, file="../OBS/MET/met_URG.Rdata")
