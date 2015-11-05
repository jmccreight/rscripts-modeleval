SetupMap <- function(file) {
        ncid <- ncdf4::nc_open(file)
        lats<-ncdf4::ncvar_get(ncid, "XLAT_M")
        lons<-ncdf4::ncvar_get(ncid, "XLONG_M")
        ncdf4::nc_close(ncid)
        myLocation <- c(min(lons)-0.1, min(lats)-0.1, max(lons)+0.1, max(lats)+0.1)
        myMap <- ggmap::get_map(location=myLocation, source="stamen", maptype="terrain", crop=TRUE)
	myMap
}

PlotMapErrors <- function(myMap, statsObj, 
			statsTag, statsVar, statsSeas, 
			plotTitle="Model Errors", plotSubTitle="",
			sizeVar="t_mae", colorVar="t_msd",
			sizeLab="Mean Absolute Error", colorLab="Mean Signed Deviation",
			colorLow="blue", colorMid="white", colorHigh="red") {
	if (is.null(statsVar)) {
		myData <- subset(statsObj, statsObj$tag==statsTag & statsObj$seas==statsSeas)
	} else {
		myData <- subset(statsObj, statsObj$tag==statsTag & statsObj$var==statsVar & statsObj$seas==statsSeas)
	}
	xCol <- ifelse("lon" %in% names(statsObj), "lon", "st_lon")
	yCol <- ifelse("lat" %in% names(statsObj), "lat", "st_lat")
	gg <- ggmap::ggmap(myMap) + 
		ggplot2::geom_point(aes_string(x=xCol, y=yCol, size=sizeVar, fill=colorVar), data=myData, alpha=0.8, shape=21) + 
		ggplot2::scale_size(sizeLab, range=c(2,10)) + 
		ggplot2::scale_fill_gradient2(colorLab, low=colorLow, mid=colorMid, high=colorHigh, midpoint=0) +
		ggplot2::ggtitle(bquote(atop(.(plotTitle), atop(italic(.(plotSubTitle)), "")))) +
		ggplot2::theme(plot.title = element_text(size=18, face="bold", vjust=-1))
	gg
}
