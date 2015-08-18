library(rwrfhydro)
sno.URG.data <- GetSnotel(c(sno.URG.sites$site_id), series="Daily", startYr=2004, endYr=2015, duration="CY")
sno.inds <- GetGeogridIndex(data.frame(lon=sno.URG.sites$lon, lat=sno.URG.sites$lat), "../DOMAIN/geo_em.d02.nc")
sno.URG.sites <- cbind(sno.URG.sites, sno.inds)

dg.URG.sites <- read.table("../OBS/DG_obs/dg_sites_latlon.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
dg.inds <- GetGeogridIndex(data.frame(lon=dg.URG.sites$lon, lat=dg.URG.sites$lat), "../DOMAIN/geo_em.d02.nc")
dg.URG.sites <- cbind(dg.URG.sites, dg.inds)

save(sno.URG.sites, sno.URG.data, dg.URG.sites, file="../OBS/snotel_URG.Rdata")