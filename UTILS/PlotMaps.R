SetupMap <- function(file) {
        ncid <- ncdf4::nc_open(file)
        lats<-ncdf4::ncvar_get(ncid, "XLAT_M")
        lons<-ncdf4::ncvar_get(ncid, "XLONG_M")
	clat <- mean(lats)
	clon <- mean(lons)
        ncdf4::nc_close(ncid)
        myLocation <- c(min(lons), min(lats), max(lons), max(lats)-5.0)
        myMap <- ggmap::get_map(location=myLocation, source="stamen", maptype="terrain", crop=TRUE)
	#myMap <- ggmap::get_googlemap(center=c(lon=clon, lat=clat), zoom=4, maptype="terrain", format="png8", size=c(640,480), scale=2)
	myMap
}

PlotMapErrors <- function(myMap, statsObj, 
			plotTitle="Model Errors", plotSubTitle="",
			sizeVar="t_mae", colorVar="t_msd",
			sizeLab="Mean Absolute Error", colorLab="Mean Signed Deviation",
			colorLow="blue", colorMid="white", colorHigh="red",
			minThreshSize=NULL, maxThreshSize=NULL,
			minThreshCol=NULL, maxThreshCol=NULL,
			minPtsize=1, maxPtsize=8,
			exclVar="t_n", exclThresh=0.8,
			colBreaks, 
			valBreaks) {
	#if (is.null(statsVar)) {
	#	myData <- subset(statsObj, statsObj$tag==statsTag & statsObj$seas==statsSeas)
	#} else {
	#	myData <- subset(statsObj, statsObj$tag==statsTag & statsObj$var==statsVar & statsObj$seas==statsSeas)
	#}
	#maxCnt <- max(myData[,exclVar], na.rm=TRUE)
	myData <- subset(statsObj, statsObj[,exclVar] >= exclThresh)
	#myData <- subset(myData, myData[,sizeVar]>quantile(myData[,sizeVar], 0.05, na.rm=TRUE) & 
	#		myData[,sizeVar]<quantile(myData[,sizeVar], 0.95, na.rm=TRUE))
	if (is.null(minThreshSize)) minThreshSize <- min(myData[,sizeVar], na.rm=TRUE)
	if (is.null(maxThreshSize)) maxThreshSize <- max(myData[,sizeVar], na.rm=TRUE)
	if (is.null(minThreshCol)) minThreshCol <- min(myData[,colorVar], na.rm=TRUE)
	if (is.null(maxThreshCol)) maxThreshCol <- max(myData[,colorVar], na.rm=TRUE)
	xCol <- ifelse("lon" %in% names(statsObj), "lon", "st_lon")
	yCol <- ifelse("lat" %in% names(statsObj), "lat", "st_lat")
	myData$plotcol <- cut(myData[,colorVar], breaks = valBreaks, right = TRUE)
	myData <- subset(myData, !is.na(myData$plotcol))
	valBreaksScaled <- scales::rescale(valBreaks, from=range(myData[,colorVar], na.rm = TRUE, finite = TRUE))
	gg <- ggmap::ggmap(myMap) + 
		ggplot2::geom_point(aes_string(x=xCol, y=yCol, size=sizeVar, fill="plotcol"), data=myData, alpha=0.8, shape=21) + 
		ggplot2::scale_size(sizeLab, range=c(minPtsize, maxPtsize), limits=c(minThreshSize, maxThreshSize)) +
		ggplot2::scale_fill_manual(colorLab, values=colBreaks) +
		#ggplot2::scale_fill_gradient2(colorLab, low=colorLow, mid=colorMid, high=colorHigh, midpoint=0, limits=c(minThreshCol, maxThreshCol)) +
		#ggplot2::scale_fill_gradientn(colorLab, colours = colBreaks, values = scales::rescale(valBreaks)) +
		ggplot2::ggtitle(bquote(atop(.(plotTitle), atop(italic(.(plotSubTitle)), "")))) +
		ggplot2::theme(plot.title = element_text(size=18, face="bold", vjust=-1)) +
		ggplot2::guides(fill = guide_legend(override.aes = list(size=3)))
	freqtbl <- table(myData$plotcol)
	#myDataSub <- subset(myData, !is.na(myData$plotcol))
	gghist <- ggplot2::ggplot(data=myData, aes(plotcol, fill=plotcol)) + 
			ggplot2::geom_histogram() +
			ggplot2::labs(x=colorLab, y="Site Count") +
			ggplot2::ggtitle(paste0("Distribution of ", colorLab)) +
			ggplot2::scale_fill_manual(colorLab, values=colBreaks) 
                	#ggplot2::theme(plot.title = element_text(size=12, face="bold", vjust=1))
	list(gg, freqtbl, gghist)
}
