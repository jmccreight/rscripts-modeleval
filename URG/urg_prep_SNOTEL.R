library(rwrfhydro)
sno.URG.data <- GetSnotel(c(sno.URG.sites$site_id), series="Daily", startYr=2004, endYr=2015, duration="CY")
sno.inds <- GetGeogridIndex(data.frame(lon=sno.URG.sites$lon, lat=sno.URG.sites$lat), "/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/geo_em.d02.nc")
sno.URG.sites <- cbind(sno.URG.sites, sno.inds)
sno.URG.data$POSIXct <- as.POSIXct(paste0(sno.URG.data$Date, " 00:00"), format="%Y-%m-%d", tz="UTC")
sno.URG.data$Tobs_K <- sno.URG.data$Tobs_C+273.15
sno.URG.data$Tmax_K <- sno.URG.data$Tmax_C+273.15
sno.URG.data$Tmin_K <- sno.URG.data$Tmin_C+273.15
sno.URG.data$Tavg_K <- sno.URG.data$Tavg_C+273.15


met.URG.sites <- read.table("../OBS/DG_obs/dg_sites_latlon.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
met.inds <- GetGeogridIndex(data.frame(lon=met.URG.sites$lon, lat=met.URG.sites$lat), "/glade/p/ral/RHAP/gochis/Col_Upp_Rio_Grande/DOMAIN/geo_em.d02.nc")
met.URG.sites <- cbind(met.URG.sites, met.inds)

save(sno.URG.sites, sno.URG.data, met.URG.sites, file="../OBS/snotel_URG.Rdata")
