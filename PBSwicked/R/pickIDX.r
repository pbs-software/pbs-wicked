pickIDX <- function(survey_lab="SYN QCS", species="silvergray rockfish",
   datnam, idxnam, outnam, plot.maps=FALSE, 
   png=FALSE, pngres=400, PIN=c(8,6), lang=c("e"))
{
	vomit <- function() {
		gc(verbose=FALSE)
		while (any(grepl("file:", search())))
			detach(grep("file:", search(),value=TRUE)[1], character.only=TRUE)
	}
	on.exit( vomit() )

	## Need to attach the binary data and index objects
	if (missing(datnam))
		datnam = paste0("data.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (file.exists(paste0(datnam,".rda"))) {
		attach(paste0(datnam,".rda"))
	} else {
		stop("Run the models first using function 'wickSSD()' to get suitable data")
	}
	if (missing(idxnam))
		idxnam = paste0("index.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (file.exists(paste0(idxnam,".rda"))) {
		attach(paste0(idxnam,".rda"))
	} else {
		stop("Run the models first using function 'wickTMB()' to get a suitable index")
	}
	## Plotting using wicked Wickham's routines
	theme_set(ggsidekick::theme_sleek())

	plot_map <- function(dat, column) {
		ggplot(dat, aes(X, Y, fill = {{ column }})) +
		geom_raster() +
		facet_wrap(~year) +
		coord_fixed()
	}
	## 'predictions' created in 'wickTMB()'
	## may need to move these function calls inside language loop if french is needed
	fig.FE.RE <- plot_map(predictions$data, exp(est2)) +
		scale_fill_viridis_c(trans = "sqrt") +
		ggtitle("Prediction (fixed effects + all random effects)")

	fig.FE <- plot_map(predictions$data, exp(est_non_rf2)) +
		ggtitle("Prediction (fixed effects only)") +
		scale_fill_viridis_c(trans = "sqrt")

	fig.S.RE <- plot_map(predictions$data, omega_s2) +
		ggtitle("Spatial random effects only") +
		scale_fill_gradient2()

	fig.ST.RE <- plot_map(predictions$data, epsilon_st2) +
		ggtitle("Spatiotemporal random effects only") +
		scale_fill_gradient2()

	if (png) createFdir(lang, dir=".")

	if (plot.maps) {
		fout.map = paste0("map.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
		## Predictions
		map.effects = c("FE.RE", "FE", "S.RE", "ST.RE")
		for (i in 1:length(map.effects)) {
			ii = map.effects[i]
			for (l in lang) {
				changeLangOpts(L=l)
				fout = switch(l, 'e' = paste0("./english/",fout.map, ".", ii), 'f' = paste0("./french/",fout.map, ".", ii) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(file=paste0(fout,".png"), width=9, height=8.5, units="in", res=pngres)
				}
				mess = paste0("print(fig.", ii, ")")
				eval(parse(text=mess))
				if (png) dev.off()
			}; eop()  ## end lang loop
#browser();return()
		}
		#print(fig.FE)
		#print(fig.S.RE)
		#print(fig.ST.RE)
	}
	## Raw Index
	#rawindex <- gfdata::get_survey_index(species, ssid = ssid)  ## now run in 'wickSSD()'
	raw      <- as.data.frame(cbind(rawindex$year,rawindex$biomass,rawindex$lowerci,rawindex$upperci))
	colnames(raw) <- c('year','est','lwr','upr')
	raw$type <- 'Design'
	index$type <- 'Model'
	index_DM <- full_join(raw, index) ## D=design, M=model

	if (missing(outnam))
		outnam = paste0("plot.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	save("fig.FE.RE", "fig.FE", "fig.S.RE", "fig.ST.RE", "raw", "index_DM", file=paste0(outnam,".rda"))

	## Plotting by geomean
	## -------------------
	## Start with design-based points
	xdes = index_DM$year[index_DM$type=="Design"]-0.1
	ydes = index_DM$est[index_DM$type=="Design"]/exp(mean(log(index_DM$est[index_DM$type=="Design"])))
	ydlo = index_DM$lwr[index_DM$type=="Design"]/exp(mean(log(index_DM$est[index_DM$type=="Design"])))
	ydhi = index_DM$upr[index_DM$type=="Design"]/exp(mean(log(index_DM$est[index_DM$type=="Design"])))
	## Gather geostatistical points
	xgeo = index_DM$year[index_DM$type=="Model"]
	ygeo = index_DM$est[index_DM$type=="Model"]/exp(mean(log(index_DM$est[index_DM$type=="Model"])))
	yglo = index_DM$lwr[index_DM$type=="Model"]/exp(mean(log(index_DM$est[index_DM$type=="Model"])))
	yghi = index_DM$upr[index_DM$type=="Model"]/exp(mean(log(index_DM$est[index_DM$type=="Model"])))
	xlim = c(2002,2025); ylim=c(0,max(ydhi))

	fout.e = outnam
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(file=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
		}
		expandGraph(mfrow=c(1,1), mar=c(3,3,0.5,0.5), mgp=c(1.8,0.5,0))
		plot(xdes, ydes, type="n", xlab="Year", xlim=xlim, ylim=ylim, ylab="Index", las=1, cex.axis=1.2, cex.lab=1.5)
		axis(1, at=xlim[1]:xlim[2], tcl=-0.2, labels=FALSE)
		axis(2, at=seq(ylim[1],ylim[2],0.1), tcl=-0.15, labels=FALSE)
		polygon(x=c(xgeo,rev(xgeo)), y=c(yghi,rev(yglo)), col=lucent("blue",0.1), border="gainsboro") #paleturquoise powderblue
		
		## Add design-based bounds then points
		arrows(xdes, ydlo, xdes, ydhi, angle=0, lwd=2, pch=15, col="black")
		points(xdes, ydes, col="red", bg="orangered", pch=22, cex=1.4)
	
		## Add geostatistical points and confidence bounds
		lines(xgeo, ygeo, lwd=2, col="blue")
		points(xgeo, ygeo, pch=21, cex=1.5, col="blue", bg="dodgerblue")
		lines(xgeo, yglo, lty=2, col="blue")
		lines(xgeo, yghi, lty=2, col="blue")
	
		legend("topright", legend=c("Geostatistical","Design-based"), col=c("blue","black"), lwd=1, pch=c(21,22), pt.bg=c("dodgerblue","orangered"), pt.cex=c(1.5,1.4), seg.len=3, bty="n", inset=0.05)
		addLabel(0.025, 0.975, gsub("\\_"," ",substring(mnam,4)), col="slategray", cex=0.8, adj=c(0,1))
		if (png) dev.off()
	}; eop()  ## end lang loop
#browser();return()
}

pickIDX(plot.maps=T, png=T)

