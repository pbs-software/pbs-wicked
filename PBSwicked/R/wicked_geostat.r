##==============================================================================
## Module : wicked_geostat
## -----------------------
## wickSSD..........Get survey series data (SSD) for a species
## wickTMB..........Run TMB geostatistical models and choose the best one
## pickIDX..........Plot results of fitting the best TMB model
## tmbModels........List of 14 TMB models to fit geostat indices

##-----Supplementary hidden functions-----
##===============================================================================

## wickSSD -----------------------------2026-05-28
## Get survey series data (SSD) for a species
##  using wicked Hadley Wickham code
## Code from NF to geo-standardize synoptic indices
## ------------------------------------------NF|RH
wickSSD <- function(survey_lab="SYN QCS", species="silvergray rockfish",
   datnam, plot=TRUE)
{
	if (missing(datnam))
		datnam = paste0("data.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (!file.exists(paste0(datnam,".rda"))) {
		eval(parse(text="ssids = gfdata::get_ssids()"))
		if(survey_lab=="SYN QCS"){
			ssid<-1
		} else if (survey_lab=="SYN HS"){
			ssid<-3
		} else if (survey_lab=="SYN WCVI"){
			ssid<-4
		} else if (survey_lab=="SYN WCHG"){
			ssid<-16
		}
		##  1 = Queen Charlotte Sound Synoptic Bottom Trawl
		##  3 = Hecate Strait Synoptic Bottom Trawl
		##  4 = West Coast Vancouver Island Synoptic Bottom Trawl
		## 16 = West Coast Haida Gwaii Synoptic Bottom Trawl SYN WCHG (8km mesh)

		## Executed through query 'get-survey-index.sql', which relies on GFBioSQL.dbo.BOOT_HEADER (not always up to date)
		.flush.cat("Running query to get design-based index from BOOT_HEADER and BOOT_DETAIL\n")
		eval(parse(text="rawindex <- gfdata::get_survey_index(species, ssid = ssid)"))  ## do this now rather than later in 'wickIDX()'

		## Get data, note that WCHG removed 2014 data because of < half of tows were completed. 
		.flush.cat("Running query to get survey set data (SSD) -- takes forever\n")
		eval(parse(text="survindex <- gfdata::get_survey_sets(species, ssid = ssid)"))  ## is this the query that takes forever?
		.flush.cat("Finished query to get SSD\n")
		save("ssids", "rawindex", "survindex", file=paste0(datnam,".rda"))  ## save query results temporarily in case following code bugs out
	} else {
		load(paste0(datnam, ".rda"))
	}
	surv = survindex
	if (survey_lab=="SYN WCHG"){
		surv <- surv[-which(surv$year==2014),]
	}
	if (survey_lab=="SYN WCVI" && species=="widow rockfish"){
		surv <- surv[-which(surv$year==2010),]
	}
	surv$present <- NA
	surv$present[surv$density_kgpm2>0]  <- 1
	surv$present[surv$density_kgpm2==0] <- 0
	surv$area_swept  = surv$tow_length_m * surv$doorspread_m 
	surv$area_swept1 = surv$doorspread_m  * (surv$speed_mpm  * surv$duration_min)
	surv$area_swept  = ifelse(is.na(surv$area_swept),surv$area_swept1,surv$area_swept)
	surv$log_area_swept = log(surv$area_swept/1e5)
	surv_ini = surv ## keep initial survey results before subsetting (RH 260213)

	## Inexplicably overwrites 'surv'
	surv <- as.data.frame(cbind(surv$year,surv$catch_weight,surv$log_area_swept,surv$present,surv$depth_m,surv$longitude,surv$latitude))
	colnames(surv) <- c('year','catch_weight','log_area_swept','present','depth','X','Y')

	## Code to fill in missing depths
	eval(parse(text="surv_wutm <- gfplot:::ll2utm(surv, utm_zone = 9)"))
	## The following (originally) seemed to assume that there are missing depths (RH 260213)
	## Treat zero depths as missing
	bad_depths = surv_wutm$depth<=0 | is.na(surv_wutm$depth)
	if (any(bad_depths)) {
		surv_wutm$depth[bad_depths] = NA
		eval(parse(text="get_depth <- gfplot:::interp_survey_bathymetry(surv_wutm)"))
		## Line that actually fills it in
		surv$depth[bad_depths] <- surv_wutm$depth[bad_depths] <- get_depth$data$akima_depth[bad_depths]
	}
	## logging and standardizing depth
	surv$log_depth_c  <- (log(surv$depth)-mean(log(surv$depth)))/sd(log(surv$depth))
	surv$log_depth_c2 <- surv$log_depth_c^2
#browser();return()
	
	## Getting UTMs back in
	eval(parse(text="surv_wUTM <- gfplot:::ll2utm(surv, utm_zone = 9)"))

	if (survey_lab=="SYN WCHG"){
		mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 8)  #10km mesh except for WCHG, which had mesh cutoff of 8km, these cutoffs are based on default for gfsynopsis report
	} else { 
		mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 10)  #10km mesh except for WCHG, which had mesh cutoff of 8km
	}
	## save this mess around (overwrite the RDA if it was queried the first time)
	save("ssids", "rawindex", "survindex", "surv_ini", "surv", "surv_wutm", "surv_wUTM", "mesh", file=paste0(datnam,".rda"))

	## Mesh (RH 251205) cannot seem to change appearance of mesh with par() arguments
	if (plot) {
		##require(fmesher) ## functions needed imported from fmesher in package PBSwicked
		expandGraph()
		plot(mesh)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~wickSSD


## wickTMB -----------------------------2026-05-28
## Run TMB geostatistical models and choose the best one;
##  uses wicked Hadley Wickham code
## Code from NF to geo-standardize synoptic indices
## ------------------------------------------NF|RH
wickTMB <- function(survey_lab="SYN QCS", species="silvergray rockfish",
   tmbModels, datnam, modnam, idxnam, verbose=FALSE, plot=TRUE)
{
	## Need to attach the survey series binary data object
	if (missing(datnam))
		datnam = paste0("data.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (file.exists(paste0(datnam,".rda"))) {
		attach(paste0(datnam,".rda"))
	} else {
		stop("Need object 'surv_wUTM'; extract data for survey and species using 'wickSSD()'")
	}
	surv_wUTM_TMB = surv_wUTM = surv_wUTM
	mesh = mesh
	detach(2)  ## second place because it was just attached

	if (missing(tmbModels)) {
		ici = lenv()
		data("tmbModels", package="PBSwicked", envir=ici)
	}
	if (missing(modnam))
		modnam = paste0("model.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (!file.exists(paste0(modnam,".rda"))) {
		.flush.cat("Started fitting ", length(tmbModels), " models -- takes forever\n", sep="")
		modres    <- lapply(tmbModels[1:14], function(fun) { fun(surv_wUTM, mesh) })  ## Takes about 30 minutes
		save("modres", file=paste0(modnam, ".rda"))
	} else {
		load(paste0(modnam, ".rda"))
	}
	## Remove the non-converged models
	modres = modres[!is.na(modres)]
	## Get the AIC values for the 14 TMB models
	AIC_table <- as.data.frame(matrix(sapply(modres,AIC), ncol=1, dimnames=list(model=names(modres), stat='AIC')))
#browser();return()
	AIC_order <- AIC_table[order(AIC_table$AIC),,drop=FALSE]
	## Remove models with AIC=NA
	good_AIC = !is.na(AIC_order$AIC)
	AIC_order = AIC_order[good_AIC,,drop=FALSE]
	print(AIC_order) ## always want to see this
	nmods  = nrow(AIC_order)

	## Basically troll through models until sanity is achieved
	sane = FALSE; eye=0
	.flush.cat("Checking model sanity\n", sep="")
	while (!sane) {
		eye  = eye + 1
		if (eye>nmods) stop ("ALL MODELS INSANE !!!")
		.flush.cat("...", rownames(AIC_order)[eye], "\n", sep="")
		mod  = substring(rownames(AIC_order)[eye], 1, 2)
		#m    = modres[[as.numeric(mod)]] ## won't work when some models have not converged
		mpos = grep(mod, names(modres))
		m    = modres[[mpos]]
		mnam = names(modres)[mpos]
#browser();return()
		exam = try( suppressMessages(sanity(m, silent=TRUE)), silent=TRUE )
		if (inherits(exam,"try-error"))
			sane = FALSE
		else
			sane = all(unlist(exam))
	}
	if (verbose) {
		m$sd_report
		summary(m)
		tidy(m,conf.int=TRUE)
		tidy(m,effects='fixed')
		tidy(m,effects='ran_pars',conf.int=TRUE)
	}
	set.seed(5)
#browser();return()
	## Look at residuals
	#  If it's a two-step model, 1=presence-absence, 2=positive catch rates
	#surv_wUTM_TMB$resids <- residuals(m, type="mle-mvn", model=2)
	surv_wUTM_TMB$resids <- residuals(m, type="mle-mvn")
	if (plot) {
		expandGraph()
		qqnorm(surv_wUTM_TMB$resids)
		abline(a=0,b=1)
	}
	## Create extrapolation grid for prediction to sum over to make eventual index 
	eval(parse(text="survey_grid <- gfplot::synoptic_grid |>
		filter(survey == survey_lab) |>
		select(X, Y, area = cell_area, depth=depth)"))

	eval(parse(text="nd <- sdmTMB::replicate_df(dat = survey_grid, time_name = \"year\", time_values = unique(surv_wUTM_TMB$year))"))
	## Making the standardized variable for the prediction. Remember this is using the mean log depth of the original dataset
	nd$log_depth_c <- (log(nd$depth)-mean(log(surv_wUTM_TMB$depth)))/sd(log(surv_wUTM_TMB$depth))
	nd$log_depth_c2 <- nd$log_depth_c^2

	## Predicting on new grid
	predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
	## Getting the index
	index <- get_index(predictions, area = nd$area, bias_correct = TRUE)
	## Here is the data you would want for your assessment
	if (verbose){
		print(index)
	}
	if (missing(idxnam))
		idxnam = paste0("index.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	save("index", "predictions", "nd", "survey_grid", "surv_wUTM_TMB", "AIC_order", "m", "mnam", file=paste0(idxnam,".rda"))
#browser();return()
	return(index) ## data for geostatistical index
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~wickTMB


## pickIDX -----------------------------2026-05-28
## Plot the results of fitting the best TMB model;
##   uses wicked Hadley Wickham code
## Code from NF to geo-standardize synoptic indices
## ------------------------------------------NF|RH
pickIDX <- function(survey_lab="SYN QCS", species="silvergray rockfish",
   datnam, idxnam, pixnam, plot.maps=FALSE, tabs=FALSE,
   png=FALSE, pngres=400, PIN=c(8,6), lang=c("f","e"))
{
	vomit <- function() {  ## detach binary RDA files
		gc(verbose=FALSE)
		while (any(grepl("file:", search())))
			detach(grep("file:", search(),value=TRUE)[1], character.only=TRUE)
	}
	on.exit( vomit() )

	## Get three-letter species code
	strSpp = switch(species,
		'silvergray rockfis' = "SGR",
		'widow rockfish'     = "WWR",
		"UNK" )
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
	eval(parse(text="theme_set(ggsidekick::theme_sleek())"))  ## SA construct

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
		}
	}
	## Raw Index
	#rawindex <- gfdata::get_survey_index(species, ssid = ssid)  ## now run in 'wickSSD()'
	raw           <- as.data.frame(cbind(rawindex$year,rawindex$biomass,rawindex$lowerci,rawindex$upperci))
	colnames(raw) <- c('year','est','lwr','upr')
	raw$type      <- 'Design'
	index$type    <- 'Model'
	index_DM      <- full_join(raw, index) ## D=design, M=model

	if (missing(pixnam))
		pixnam = paste0("plot.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])

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

	isgeo = match(xgeo, xdes+0.1)
	xdes = xdes[isgeo]; ydes=ydes[isgeo]; ydlo=ydlo[isgeo]; ydhi=ydhi[isgeo]
	index_GM = data.frame(year=xgeo, ygeo=ygeo, yglo=yglo, yghi=yghi, ydes=ydes, ydlo=ydlo, ydhi=ydhi)

	save("fig.FE.RE", "fig.FE", "fig.S.RE", "fig.ST.RE", "index_DM", "index_GM", file=paste0(pixnam,".rda"))
	## Schmush together for output
	if (tabs) {
		write.csv (index_DM, paste0(sub("plot","indexDM",pixnam),".csv"), row.names=F)
		write.csv (index_GM, paste0(sub("plot","indexGM",pixnam),".csv"), row.names=F)
	}
	xlim = c(2002,2025); ylim=c(0,max(ydhi))
	xlim = range(xgeo) + c(-0.5,0.5); ylim=c(0,max(c(yghi,ydhi), na.rm=TRUE))

	fout.e = pixnam
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(file=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
		}
		expandGraph(mfrow=c(1,1), mar=c(3,3,0.5,1), mgp=c(1.8,0.5,0))
		plot(xdes, ydes, type="n", xlab=linguaFranca("Year",l), xlim=xlim, ylim=ylim, ylab=linguaFranca("Index",l), las=1, cex.axis=1.2, cex.lab=1.5)
		axis(1, at=xgeo[1]:xgeo[length(xgeo)], tcl=-0.2, labels=FALSE)
		by = ifelse(as.logical(floor(max(ylim))%%2),1,2) * 10^floor(log10(max(ylim))-1)
		axis(2, at=seq(ylim[1], ylim[2], by=by), tcl=-0.15, labels=FALSE)
		polygon(x=c(xgeo,rev(xgeo)), y=c(yghi,rev(yglo)), col=lucent("blue",0.1), border="gainsboro") #paleturquoise powderblue
		
		## Add design-based bounds then points
		arrows(xdes, ydlo, xdes, ydhi, angle=0, lwd=2, pch=15, col="black")
		points(xdes, ydes, col="red", bg="orangered", pch=22, cex=1.4)
	
		## Add geostatistical points and confidence bounds
		lines(xgeo, ygeo, lwd=2, col="blue")
		points(xgeo, ygeo, pch=21, cex=1.5, col="blue", bg="dodgerblue")
		lines(xgeo, yglo, lty=2, col="blue")
		lines(xgeo, yghi, lty=2, col="blue")
	
		addLegend(0.05, 0.95, legend=linguaFranca(c("Geostatistical","Design-based"),l), col=c("blue","black"), lwd=1, pch=c(21,22), pt.bg=c("dodgerblue","orangered"), pt.cex=c(1.5,1.4), seg.len=3, bty="n", xjust=0)
		labtxt = paste0(strSpp, " (", survey_lab, ") TMB model ", gsub("\\_"," ", sub("\\_"," : ", mnam)))
		addLabel(0.025, 0.975, txt=linguaFranca(labtxt,l), col="slategray", cex=0.8, adj=c(0,1))
		if (png) dev.off()
	}; eop()  ## end lang loop
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pickIDX


## tmbModels----------------------------2025-12-12
##  NB : Fitting 14 models.
##  The authors of the gfsynopsis report choose the likelihood to
##  use as the one that passes all sanity checks and has lowest AIC.
##  (RH 251205) Crikey! this will take forever!
## ------------------------------------------NF|RH
tmbModels <- list(
'01_Delta_Lognormal' =
	m_delta_lognormal <- function(surv_wUTM, mesh) {
	.flush.cat("   01_Delta_Lognormal\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_lognormal(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'02_Delta_Lognormal_with_depth' = 
	m_delta_lognormal_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   02_Delta_Lognormal_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_lognormal(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'03_Tweedie' = 
	m_tweedie <- function(surv_wUTM, mesh) {
	.flush.cat("   03_Tweedie\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = tweedie(link = "log"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'04_Tweedie_with_depth' = 
	m_tweedie_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   04_Tweedie_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = tweedie(link = "log"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'05_Delta_gamma' = 
	m_delta_gamma <- function(surv_wUTM, mesh) {
	.flush.cat("   05_Delta_gamma\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'06_Delta_gamma_with_depth' = 
	m_delta_gamma_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   06_Delta_gamma_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'07_Delta_poisson_link_gamma' = 
	m_delta_poislinkgamma <- function(surv_wUTM, mesh) {
	.flush.cat("   07_Delta_poisson_link_gamma\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'08_Delta_poisson_link_gamma_with_depth' = 
	m_delta_poislinkgamma_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   08_Delta_poisson_link_gamma_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'09_Delta_poisson_link_lognormal' = 
	m_delta_poislinklognorm <- function(surv_wUTM, mesh) {
	.flush.cat("   09_Delta_poisson_link_lognormal\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'10_Delta_poisson_link_lognormal_with_depth' = 
	m_delta_poislinklognorm_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   10_Delta_poisson_link_lognormal_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'11_Delta_generalized_gamma_distribution' = 
	m_delta_gengamma <- function(surv_wUTM, mesh) {
	.flush.cat("   11_Delta_generalized_gamma_distribution\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gengamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'12_Delta_generalized_gamma_distribution_with_depth' = 
	m_delta_gengamma_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   12_Delta_generalized_gamma_distribution_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gengamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'13_Delta_poisson_link_generalized_gamma_distribution' = 
	m_delta_poislinkgengamma <- function(surv_wUTM, mesh) {
	.flush.cat("   13_Delta_poisson_link_generalized_gamma_distribution\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	},
'14_Delta_poisson_link_generalized_gamma_distribution_with_depth' = 
	m_delta_poislinkgengamma_depth <- function(surv_wUTM, mesh) {
	.flush.cat("   14_Delta_poisson_link_generalized_gamma_distribution_with_depth\n")
	out = try(sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=TRUE,
	  control = sdmTMBcontrol(newton_loops = 1)
	), silent=TRUE)
	if (inherits(out,"try-error")) return(NA) else return(out)
	}
)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tmbModels


