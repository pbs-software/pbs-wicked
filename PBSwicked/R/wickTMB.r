wickTMB <- function(survey_lab="SYN QCS", species="silvergray rockfish",
   tmbModels, datnam, modnam, outnam, verbose=TRUE, plot=TRUE)
{
	if (missing(modnam))
		modnam = paste0("model.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (!file.exists(paste0(modnam,".rda"))) {
		modres    <- lapply(tmbModels, function(fun) { fun(surv_wUTM) })  ## Takes about 30 minutes
		save("modres", file=paste0(modnam, ".rda"))
	} else {
		load(paste0(modnam, ".rda"))
	}
	## Get the AIC values for the 14 TMB models
	AIC_table <- as.data.frame(matrix(sapply(modres,AIC), ncol=1, dimnames=list(model=names(modres), stat='AIC')))
	AIC_order <- AIC_table[order(AIC_table$AIC),,drop=F]

	## Basically troll through models until sanity is achieved
	sane = FALSE; eye=0
	while (!sane) {
		eye  = eye + 1
		mod  = substring(rownames(AIC_order)[eye], 1, 2)
		m    = modres[[as.numeric(mod)]]
		mnam = names(modres)[as.numeric(mod)]
		exam = sanity(m)
		sane = all(unlist(exam))
	}
	if (verbose) {
		m$sd_report
		summary(m)
		tidy(m,conf.int=TRUE)
		tidy(m,effects='fixed')
		tidy(m,effects='ran_pars',conf.int=TRUE)
	}
	#Look at residuals
	## Need to attach the survey series binary data object
	if (missing(datnam))
		datnam = paste0("data.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	attach(paste0(datnam,".rda"))
	surv_wUTM_TMB = surv_wUTM
	detach(2)  ## second place because it was just attached
#browser();return()
	
	set.seed(5)
	#If its a two-step model, 1=presence-absence, 2=positive catch rates
	surv_wUTM_TMB$resids <- residuals(m, type="mle-mvn", model=2)
	#surv_wUTM_TMB$resids <- residuals(m, type="mle-mvn")
	if (plot) {
		qqnorm(surv_wUTM_TMB$resids)
		abline(a=0,b=1)
	}
	## Create extrapolation grid for prediction to sum over to make eventual index 
	survey_grid <- gfplot::synoptic_grid |>
		filter(survey == survey_lab) |>
		select(X, Y, area = cell_area, depth=depth)

	nd <- sdmTMB::replicate_df(dat = survey_grid, time_name = "year", time_values = unique(surv_wUTM_TMB$year))
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
#browser();return()
	if (missing(outnam))
		outnam = paste0("index.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	save("index", "predictions", "nd", "survey_grid", "surv_wUTM_TMB", "m", "mnam", file=paste0(outnam,".rda"))
	return(index) ## data for geostatistical index
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~wickTMB

require(sdmTMB)
require(sp)
require(ggplot2)
require(dplyr)
require(akima)
source("wickSSD.r")  ## need data from scratch or loaded from binary
source("tmbModels.r")
wickTMB(survey_lab="SYN QCS", species="silvergray rockfish", tmbModels=tmbModels)
