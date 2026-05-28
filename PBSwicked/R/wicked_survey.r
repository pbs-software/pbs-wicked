##==============================================================================
## Module : wicked_survey
## -----------------------
## calcAE......Calculate ageing error using the function 'RunFn()' from 'AgeingError'
## getHBLL.....Get survey series indices for HBLL surveys
## predREBS....Predict probability of RER (BSR is the inverse)

##-----Supplementary hidden functions-----
##===============================================================================
## Note: not sure whether to rename functions so it's obvious that they are wicked
## See object `wickiverb' for possible guidance (RH 260506)

## calcAE-------------------------------2026-05-28
##  Calculate ageing error using the function 'run()'
##  from the R package 'AgeingError'.
## Note: function 'RunFn()' now deprecated (Ian Taylor, pers.commm.)
##  Adapted from Kendra Holt's code used in Petrale 2024
##  Attempt to de-hadley the code
##  ageDat  : use data from 'gfb_age_precision.sql'
##  atype   : AGE_READING_TYPE_CODE 2=Primary, 3=Secondary
##  species : common species name
## ---------------------------------------------RH
calcAE <- function (ageDat, atype=c(2,3), species="silvergray rockfish",
	hadley=FALSE, maxage=60, useTMB=TRUE, debug=FALSE,
	png=FALSE, pngres=400, PIN=c(8,8), lang="e")
{
	cwd = getwd()
	on.exit(setwd(cwd))
	spcode = switch(species,
		'silvergray rockfish'="405",
		'widow rockfish'="417",
		"999"
	) ## defunct
	#spcode = "AE" ## AE=AgeingError : just keep it the same for all species

	## From gfsynopsis res doc appendix A, it seems likely that age_reading_type_code = 2 is the initial primary age, while age_reading_type_code = 3 is the precision age.

	#require(ggplot2)
	#require(dplyr)
	#require(here)
	wd = getwd()

	## Rename RH columns to mesh with Kendra Holt
	colnames(ageDat) = c("specimen_id", "year", "species_code", "species_common_name", "age_reading_type_code", "ageing_method_desc", "specimen_age", "minimum_age", "maximum_age", "employee_id", "age_reading_id", "ageing_method")

	## Shortcut for all age reader types
	if (debug) {
		delim  <- "============================================================"
		master <- crossTab(ageDat, c("specimen_id","age_reading_type_code"), "specimen_age", mean)
		.flush.cat(delim, "\n", "DEBUG : see object called 'master' -- specimen ages by age reader type", "\n", delim, "\n", sep="")
		browser()
	}

	## Change to wide format, aligned by specimen ID =====================================
	if (hadley) {
		ageDat <- ageDat %>% filter(ageing_method_desc != "OTOLITH SURFACE ONLY") # remove surface ages (hadley)
		## --- extract precision ages
		dat.precAge <- ageDat %>% filter(age_reading_id == atype[2]) %>% select(-c("species_common_name", "species_code", "employee_id", "ageing_method", "age_reading_type_code"))
		dat.precAge <- dat.precAge %>% rename(Prec_Age = specimen_age, Prec_minAge = minimum_age, Prec_maxAge = maximum_age)# %>% mutate(age_reading_type = rep("precision", nrow(dat.precAge)))
		## --- extract primary ages
		dat.primAge <- ageDat %>% filter(age_reading_id == atype[1]) %>% select(-c("species_common_name", "species_code", "employee_id", "ageing_method", "age_reading_type_code"))
		dat.primAge <- dat.primAge %>% rename(Primary_Age = specimen_age, Primary_minAge = minimum_age, Primary_maxAge = maximum_age)# %>% mutate(age_reading_type = rep("primary", nrow(dat.primAge)))

	} else {
		ageDat = ageDat[!is.element(ageDat$ageing_method_desc, "OTOLITH SURFACE ONLY"),] # remove surface ages
		## --- extract precision ages (ARID=3)
		keep = match(c("species_common_name", "species_code", "ageing_method", "age_reading_type_code", "age_reading_id"), colnames(ageDat))
		dat.precAge = ageDat[is.element(ageDat$age_reading_type_code, atype[2]), -keep]
		colnames(dat.precAge)[ match(c("specimen_age", "minimum_age","maximum_age"), colnames(dat.precAge)) ] = c("Prec_Age", "Prec_minAge", "Prec_maxAge")
		## --- extract primary ages (ARID=2)
		keep = match(c("species_common_name", "species_code", "ageing_method", "age_reading_type_code", "age_reading_id"), colnames(ageDat))
		dat.primAge = ageDat[is.element(ageDat$age_reading_type_code, atype[1]), -keep]
		colnames(dat.primAge)[ match(c("specimen_age", "minimum_age","maximum_age"), colnames(dat.primAge)) ] = c("Primary_Age", "Primary_minAge", "Primary_maxAge")
	}
	## purging hadley to this point (won't bother further because function is now in PBSwicked)
	## --- combine
	dat.bySpecimen <- left_join(dat.primAge, dat.precAge, by = c("specimen_id", "year", "ageing_method_desc"))

	## Add columns showing calculated differences in read ages 
	dat.bySpecimen <- dat.bySpecimen %>% mutate(readerDiff = Primary_Age - Prec_Age)
	
	## Limit dataset to only double-reads
	doubleReadDat <- dat.bySpecimen %>% filter(is.na(readerDiff)==FALSE)
	
	## Calculate overall proportion of double-read otoliths with the same age assigned by two readers
	doubleReadDat <- doubleReadDat %>% mutate(indReaderSame=ifelse(doubleReadDat$readerDiff==0,1,0))
	propSame_btReader  <-  nrow(doubleReadDat[doubleReadDat$indReaderSame == 1,]) / nrow(doubleReadDat) 
	
	## Calculate proportion of double-read otoliths with assigned ages less than one year off
	doubleReadDat <- doubleReadDat %>% mutate(indReader1Year=ifelse(abs(doubleReadDat$readerDiff)<=1,1,0))  
	prop1Year_btReader <- nrow(doubleReadDat[doubleReadDat$indReader1Year == 1,]) / nrow(doubleReadDat) 
	
	## Create table of proportion of double-read otoliths with the same age assigned by two readers, as a function of primary age
	count_btReader_byAge <- doubleReadDat %>%  group_by(Primary_Age = as.factor(Primary_Age)) %>% summarise(totalCount=n(), sameCount=sum(indReaderSame)) %>% mutate(propSame=sameCount/totalCount)
	
	## Calculate table of proportion of double-read otoliths with the same age assigned by two readers, as a function of year
	count_btReader_byYear <- doubleReadDat %>%  group_by(year = as.factor(year)) %>% summarise(totalCount=n(), sameCount=sum(indReaderSame)) %>% mutate(propSame=sameCount/totalCount)
	
	## Look at Sample Sizes for Double-read Pairs
	if (hadley) {
		## Precision test sample sizes, by year
		print.data.frame(count_btReader_byYear %>% select(year, N = totalCount))
		# Precision test sample sizes, by final age
		print.data.frame(count_btReader_byAge %>% select(Primary_Age, N = totalCount))
	} else {
		## Precision test sample sizes, by year
		po1 = count_btReader_byYear[,c("year","totalCount")]; colnames(po1)[2]="N"; print.data.frame(po1)
		## Precision test sample sizes, by final age
		po2 = count_btReader_byAge[,c("Primary_Age","totalCount")]; colnames(po2)[2]="N"; print.data.frame(po2)
	}
#browser();return()

	## Exploratory Plots
	## --------------------------------------------
	## Distribution of differences between readers over time
	g_DistbyYear <- ggplot(doubleReadDat,aes(x=as.factor(year), y=readerDiff)) +geom_boxplot() + geom_hline(yintercept=0, color = "red") + coord_cartesian(ylim = c(-20, 20))  + ylab("Difference between readers") + xlab("Year")
	fout.e = paste0("calcAE(", spcode, ")-DistbyYear")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e'=paste0("./english/",fout.e), 'f'=paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		print(g_DistbyYear)
		if (png) dev.off()
	}; eop()
	
	## Reader Comparison
	g_Scatter <- ggplot(doubleReadDat,aes(x=Primary_Age, y=Prec_Age)) + geom_point() + ylab("Precision Age Estimate") + xlab("Primary Age Estimate") + geom_abline(intercept = 0, slope = 1)
	fout.e = paste0("calcAE(", spcode, ")-Scatter")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		print(g_Scatter)
		if (png) dev.off()
	}; eop()
#browser();return()

	## Fit Linear Regression
	mod.lm <- lm(Prec_Age ~ Primary_Age, doubleReadDat)
	mod_summary <- summary(mod.lm)

	## Try bootstrapping ages for precision readers from primary readers
	if (debug) {
		eval(parse(text="require(boot)"))
		## 1. Define a function that returns the predictions
		boot_predict <- function(data, indices) {
			# Resample the data
			d <- data[indices, ] 
			# Refit the model
			fit.lm  <- lm(Prec_Age ~ Primary_Age, data=d)
			# Return predictions for new data
			return(predict(fit.lm, newdata=data.frame('Primary_Age'=0:maxage)))
		}
		## 2. Run the bootstrap
		results <- boot(data=doubleReadDat, statistic=boot_predict, R=10)
		## 3. View the bootstrap distribution of predictions
		.flush.cat(delim, "\n", "DEBUG : see object called 'results' -- bootstrapped ages (results$t)", "\n", delim, "\n", sep="")
		browser()
	}

	## Plot regression residuals:
	fout.e = paste0("calcAE(", spcode, ")-lm(prec~prim)-res")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		plot(doubleReadDat$Primary_Age, mod_summary$residuals)
		abline(h=0, col="red")
		if (png) dev.off()
	}; eop()
#browser();return()
	createTdir()  ## create directories called 'tables' and 'data'

	residAge <- data.frame(Resids = mod_summary$residuals, Primary_Age = doubleReadDat$Primary_Age)
	residAge <- as_tibble(residAge)
	
	sigAge <- residAge %>% group_by(Primary_Age) %>% summarize(sigAtAge=sd(Resids), n=n())
	write.csv(sigAge, paste0("./tables/calcAE(", spcode, ")-sigAge0.csv") )
	sigAge <- sigAge%>%filter(n>5)
	createTdir()
	write.csv(sigAge, paste0("./tables/calcAE(", spcode, ")-sigAge.csv") )

	## Plot SD vs primary age
	fout.e = paste0("calcAE(", spcode, ")-sd-vs-primage")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		plot(sigAge$Primary_Age, sigAge$sigAtAge, typ="l", xlab="Primary Age", ylab="Residual std. dev.")
		if (png) dev.off()
	}; eop()
	
	## NWFSC Package AgeingError (Schnute/Punt AE algorithm)
	## -----------------------------------------------------
	#require(tidyverse)
	eval(parse(text="view = PBSmodelling::view"))  ## squash hadley wickham (but perhaps don't invoke his universe)

	## Get AgeingError if needed
	isAE = getVer(print=FALSE, all=TRUE, pat="AgeingError")
	if (grepl("AgeingError", isAE[,1])) {
		# Load package Punt's package
		eval(parse(text="require(AgeingError)"))
	} else {
		isDT = getVer(print=FALSE, all=TRUE, pat="devtools")
		if (grepl("devtools", isDT[,1])) {
			eval(parse(text="require(devtools)"))
		} else {
			eval(parse(text="install.packages(\"devtools\")"))
		}
		eval(parse(text="remotes::install_github(\"pfmc-assessments/AgeingError\")"))
		## using {pak}
		# pak::pak("pfmc-assessments/AgeingError")
	}
	## This is where all runs will be located
	AEdir    <- "./data"
	Nreaders <- 2 ## primary and precision
	dat      <- doubleReadDat
	dat      <- dat %>% select(Primary_Age, Prec_Age)
	AgeReads <- as.matrix(dat)
	
	## Format data
	Nreaders = ncol(AgeReads)
	AgeReads = ifelse(is.na(AgeReads),-999,AgeReads)  ## Change NA to -999 (which the Punt software considers missing data)
	## Potentially eliminate rows that are only read once
	## These rows have no information about reading error, but are potentially informative about latent age-structure
	## It is unknown whether eliminating these rows degrades estimation of error and bias, and is currently recommended to speed up computation
	#KeepRow = ifelse(rowSums(ifelse(AgeReads==-999,0,1),na.rm=TRUE)<=1,FALSE,TRUE)
	#AgeReads = AgeReads[KeepRow,]

	## Testing exclusion of ages > xx
	AgeReads  <- as.data.frame(AgeReads)  ## looks like this will be the input for the new 'run' function in 'AgeingError'
	AgeReads1 <- AgeReads[AgeReads$Primary_Age<=maxage, ]
	AgeReads1 <- as.matrix(AgeReads)

	## Combine duplicate rows (see 'RunFn' example : James Thorson)
	AgeReads2 = rMx(c(1, AgeReads1[1,])) # AgeReads2 is the correctly formatted data object
	for(RowI in 2:nrow(AgeReads1)){
		DupRow = NA
		for(PreviousRowJ in 1:nrow(AgeReads2)){
			if(all(AgeReads1[RowI,1:Nreaders]==AgeReads2[PreviousRowJ,1:Nreaders+1])) DupRow = PreviousRowJ
		}
		if(is.na(DupRow)) AgeReads2 = rbind(AgeReads2, c(1, AgeReads1[RowI,])) # Add new row to AgeReads2
		if(!is.na(DupRow)) AgeReads2[DupRow,1] = AgeReads2[DupRow,1] + 1 # Increment number of samples for the previous duplicate
	}
	## sigopt 3, curvilinear coefficient of variation, i.e., a 3-parameter Hollings-form relationship of coefficient of variation with true age.
	## note: sigopt 2 (curvilinear SD) estimates nonsense (at least for WWR)
	BiasOpt = c(0,0)
	SigOpt <- c(3,-1) 

	## Define minimum and maximum ages for integral across unobserved ages
	MinAge   = 1
	MaxAge   = max(AgeReads)
	KnotAges = list(NA, NA)

	write.csv(AgeReads,  paste0("./tables/calcAE(", spcode, ")-AgeReads.csv") )   ## (RH 251020)
	write.csv(AgeReads1, paste0("./tables/calcAE(", spcode, ")-AgeReads1.csv") )  ## (RH 260515)
	write.csv(AgeReads2, paste0("./tables/calcAE(", spcode, ")-AgeReads2.csv") )  ## (RH 251020)

	if (useTMB) {
		## Ian Taylor suggests we use 'run' instead of 'RunFn', which is deprecated
		write_files(dat=AgeReads, dir=AEdir, biasopt=BiasOpt, sigopt=SigOpt)
		modRes <- try(run(directory=AEdir))
		if(inherits(modRes, "try-error") ) {
		message("------------------------\n'run' failed for some reason") ; print(modRes); setwd(cwd); browser(); return() }
		
		# see estimated parameters
		modPar = modRes$model$par
		# see model selection results
		modSel = modRes$output$ModelSelection
		# see ageing error matrices
		modErr = modRes$output$ErrorAndBiasArray
		Age    = 0:maxage
		SD     = modErr["SD",paste("Age",Age),"Reader 1"]
#browser();return()
	} else {  ## old executable (deprecated, no longer tested)
		SourceFile <- paste(AEdir,"nwfscAgeingError_src",sep="/")
		modRes = try(RunFn(Data=AgeReads2, SigOpt=SigOpt, BiasOpt=BiasOpt, KnotAges=KnotAges,
			NDataSets=1, MinAge=MinAge, MaxAge=MaxAge, RefAge=floor(maxage/2),
			MinusAge=1, PlusAge=maxage, MaxSd=10, SaveFile=AEdir, AdmbFile=NULL, #SourceFile,
			EffSampleSize=0, Intern=TRUE, JustWrite=FALSE, verbose=FALSE))
		if(inherits(modRes, "try-error") ||  any(grepl("Error",modRes)) ) {
			message("------------------------\n'RunFn' failed for some reason") ; print(modRes); setwd(cwd); browser(); return() }
		## This stupid function goes into AEdir and stays there!
		setwd(cwd)
		# Plot output
		# Data = AgeReads2; MaxAge = MaxAge; SaveFile = DateFile; PlotType = "PDF"
		for (l in lang) {
			changeLangOpts(L=l)
			switch(l, 
				'e' = { JTdir="./english/"; ReaderNames=c("Primary reader", "Precision reader") },
				'f' = { JTdir="./french/"; ReaderNames=c("technicien principal", eval(parse(text=deparse("technicien de pr\u{00E9}cision"))) ) }
			)
			JTfig  = paste0(ReaderNames[1], " vs ", ReaderNames[2])
			JTfigs = paste0(c(JTfig, "Estimated vs Observed Age Structure", "True vs Reads (by reader)"), ".png") 
			## Need to hack 'PlotOutputFn' to get french fully implemented
			output <- PlotOutputFn(Data=AgeReads2, MaxAge=MaxAge, SaveFile=AEdir,  PlotType="PNG", ReaderNames=ReaderNames)
			AEdir = paste0(sub("/$","",AEdir),"/")
			JTfigs = list.files(path=AEdir, pattern="\\.png$")
			JTfigs.new = paste0("calcAE{", spcode, ")-", gsub("\\s+","_",JTfigs)) ## (RH 260210)
			clearFiles(paste0(JTdir, JTfigs.new))
			#file.copy(from=paste0(AEdir,JTfigs), to=JTdir, overwrite=T, copy.date=T)
			file.copy(from=paste0(AEdir,JTfigs), to=paste0(JTdir,JTfigs.new), overwrite=T, copy.date=T)
			JTtabs = list.files(path=AEdir, pattern="\\.csv$")
			JTtabs.new = paste0("calcAE{", spcode, ")-", gsub("\\s+","_",JTtabs)) ## (RH 260210)
			clearFiles(paste0("./tables/", JTtabs.new))
			file.copy(from=paste0(AEdir,JTtabs), to=paste0("./tables/",JTtabs.new), overwrite=T, copy.date=T)
			#file.remove(paste0(AEdir,JTfigs))
		}; eop() ## lang loop
		errorEsts <- output$ErrorAndBiasArray[,,1]
		print(errorEsts)
		write.csv(errorEsts, paste0("./tables/calcAE(", spcode, ")-ageErrorArray_model.csv") )
		
		Age  <- 1:(ncol(errorEsts)-1)
		SD   <- as.numeric(errorEsts["SD",-1])
		lfit <- lm(SD ~ Age)  ## not really useful for a curvilinear 
	
		## Fit the linear model to log-transformed data
		efit_lm <- lm(log(SD) ~ Age)
		b_lm    <- coef(efit_lm)["Age"]
		a_lm    <- exp(coef(efit_lm)["(Intercept)"])
		cv_lm   <- summary(efit_lm)$sigma / mean(log(SD)) * 100 
		## Fit the non-linear model
	#browser();return()
		SD_add   <- a_lm * exp(b_lm * Age) + rnorm(length(Age), sd = summary(efit_lm)$sigma)
		efit_nls <- nls(SD_add ~ a * exp(b * Age), start = list(a = 1, b = 0.1))
		a_nls    <- coef(efit_nls)["a"]
		b_nls    <- coef(efit_nls)["b"]
		cv_nls   <- summary(efit_nls)$sigma / mean(log(SD)) * 100 
	} ## end old executable

	fout.e = paste0("calcAE(", spcode, ")-model-sd-vs-age")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		expandGraph(mfrow=c(1,1), mar=c(3.5,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plot(Age, SD, pch=19, xlab=linguaFranca("Age",l), ylab=linguaFranca("Standard deviation",l), cex.axis=1.25, cex.lab=1.5, col="gainsboro")
		abline(a=0, b=0.1, lwd=2, lty=2, col="red")
		if (useTMB) {
			legtxt = legcol = leglty = legpch = legbg = NULL  ## in case we use stuff in debug
			if (debug) {
				.flush.cat(delim, "\n", "DEBUG : add in SD vs. Age for (min to max) data", "\n", delim, "\n", sep="")
				## Attempt to create distribution from min and max age for each primary age
				page = dat.primAge
				pbins = as.list(rep(NA,maxage+1)); names(pbins)=0:maxage
				for (i in 1:nrow(page)) {
					ii = page[i, "Primary_Age"]
					if (ii>maxage) next
					iii = as.character(min(ii,maxage))
					pbins[[iii]] = c(pbins[[iii]], page[i,"Primary_minAge"]:page[i,"Primary_maxAge"])
				}
				misguided = FALSE
				if (misguided) { ## misguided attempt to plot the range at each age (but plot space is [SD,Age])
					mino = crossTab(page, c("Primary_Age"), "Primary_minAge", function(x) {xx=x[x>0]; min(xx)})
					maxo = crossTab(page, c("Primary_Age"), "Primary_maxAge", function(x) {xx=x[x>0]; max(xx)})
					mnmx = cbind(mino,maxo)  ## assume the vectors are the same length for now
					xrng = as.numeric(rownames(mnmx))
					xrng = as.vector(sapply(xrng,function(x){c(x,x,NA)}))
					yrng = as.vector(t(cbind(mnmx, rep(NA,nrow(mnmx)))))
					lines(xrng, yrng, col="slategray")
					browser()
				}
				## Calculate SDs assuming uniform distribution between amin and amax
				sd.unif = sapply(pbins, sd, na.rm=T)
				## Assuming min and max represent roughly the 0.5% and 99.5% quantiles, 99% of data falls within range
				## qnorm(0.995) is approx 2.576
				sd_from_min <- sapply(1:maxage, function(a) { if (all(is.na(pbins[[a+1]]))) NA else (a - min(pbins[[a+1]],na.rm=T)) / qnorm(0.995) })
				sd_from_max <- sapply(1:maxage, function(a) { if (all(is.na(pbins[[a+1]]))) NA else (max(pbins[[a+1]],na.rm=T) - a) / qnorm(0.995) })
				## Average them if the mean is perfectly centered
				sd.norm <- c(NA, apply(cbind(sd_from_min, sd_from_max), 1, mean, na.rm=T) )
				points(0:maxage, sd.norm, pch=17, col="purple")
				## really need to do a cumulative legend but cannot be bothered for something not useful
				legtxt = c("Approx. SD from (min, max)")
				legcol="purple"; leglty=NA; legpch=17; legbg=NA
			} ## end debug
			## Use observed precision ages by primary age
			age.obs = split(doubleReadDat$Prec_Age, doubleReadDat$Primary_Age)
			sd.obs  = sapply(age.obs,sd)
			use.obs = is.element(as.numeric(names(sd.obs)), 0:maxage)
			sd.use  = sd.obs[use.obs]
			age.use = as.numeric(names(sd.use))
			points(age.use, sd.use, pch=16, col="deepskyblue")
			## Fit the linear model to log-transformed data
			fit.lm  = lm(log(sd.use) ~ age.use)
			lines(age.use, exp(predict(fit.lm, newdata=data.frame('age.use'=age.use))), col="blue", lty=5, lwd=2) # Transform predictions back to original scale
			## Fit a non-linear model
			fit.nls = nls(sd.use ~ a * exp(b * age.use), start=list(a=1, b=0.1))
			#lines(age.use, predict(fit.nls, newdata=data.frame('age.use'=age.use)), col="blue", lty=3, lwd=2) ## not really much different than linerar fit

			legtxt = c(legtxt, "SD of precision ages at primary age", "linear model fit to log SD", "AgeingError SD fitting to CV (not points)", "CASAL CV = 0.10")
			legcol = c(legcol, "deepskyblue","blue","black","red")
			leglty = c(leglty, NA,5,1,2)
			legpch = c(legpch, 16,NA,21,NA)
			legbg  = c(legbg, NA,NA,"green",NA)
		} else {
			lines(Age, predict(efit_nls), col="chocolate1", lty=1, lwd=2)
			lines(Age, exp(predict(efit_lm)), col = "blue", lty=5, lwd = 2) # Transform predictions back to original scale
			legtxt = c(paste0("LM fitted CV ~ ",  round(cv_lm), "%"), paste0("NLS fitted CV ~ ",  round(cv_nls), "%"), "CASAL CV = 10%")
			legcol = c("chocolate1","blue","red")
			leglty = c(1,5,3)
			legpch = c(NA,NA,NA)
			legbg  = c(NA,NA,NA)
		}
		## Add the main AgeingError fit
		lines(Age, SD, col="black", lwd=2)
		points(Age, SD, pch=21, col="black", bg="green", cex=1.2)
		addLegend(0.025, 0.975, legend=legtxt, pch=legpch, pt.bg=legbg, lty=leglty, col=legcol, pt.cex=1.2, pt.lwd=1, lwd=2, seg.len=3, cex=1.2, bty="n", xjust=0)
		if (png) dev.off()
	}; eop() ## end lang loop
	if (debug) { browser(); return() }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcAE


## getHBLL------------------------------2026-05-04
## Outside HBLL Survey - Originally used for Canary
## Code supplied by Dana Haggarty (2022-04-11)
## Includes hook competition adjustment.
## Moved from PBStools to PBSwicked because code uses Hadley nonsense. (RH 260504)
## ------------------------------------------DH|RH
getHBLL <- function(strSpp="417", timestamp="260502")
{
	cook <- function(parsnip="sumting", spice="elapsed", meat="beef") {
		.flush.cat("Grabbing ", meat, " ...", sep="", "\n")
		start.time = proc.time()
		stew <- eval(parse(text=parsnip))
		end.time = proc.time()
		exlax = difftime(time1=end.time[spice], time2=start.time[spice])
#browser();return()
		msg   = paste0("..... time wasted : ", paste(round(exlax,2), units(exlax)))
		.flush.cat(msg, "\n")
		return(invisible(stew))
	}
	##-----------------------------end subfunctions

	#rm(list=setdiff(ls(all=T), c(".First",".SavedPlots","clr","clr.rgb","qu","so" )))  #Erase all previously saved objects
	data("species", package="PBSdata")
	spn = toUpper(species[strSpp,"name"])
	spc = species[strSpp,"code3"]

	#library(dplyr)
	#library(gfdata)
	createTdir()

	## Use 'save' to retain tibble nonsense
	d1.csv = paste0("./tables/", spc,"_LLhookdata_", timestamp, ".csv")
	d1.rda = paste0("./data/", spc, "_LLhookdata_", timestamp, ".rda")
	if (file.exists(d1.rda)) {
		load(d1.rda)
	} else {
		recipe = paste0("gfdata::get_ll_hook_data(species = \"", strSpp, "\" ,  ssid = c(22, 36))")
		d1 = cook(recipe, meat="LL hook and line data")
		write.csv(d1, d1.csv)
		save("d1", file=d1.rda)
	}
	## Use 'save' to retain tibble nonsense
	d2.csv = paste0("./tables/", spc, "_outsideHBLL_", timestamp, ".csv")
	d2.rda = paste0("./data/", spc, "_outsideHBLL_", timestamp, ".rda")
	if (file.exists(d2.rda)) {
		load(d2.rda)
	} else {
		## This query takes a long time to execute -- perhaps recode in future?
		recipe = paste0("gfdata::get_survey_sets(\"", spn, "\", ssid = c(22, 36))")
		d2 = cook(recipe, meat="survey sets")
		write.csv(d2, d2.csv )
		save("d2", file=d2.rda)
	}
#browser();return()

	## The following code should deal with the issue of 0 baited hooks being observed.
	adjust <- d1 %>%
		group_by(year, fishing_event_id) %>%
		mutate(total_hooks = count_target_species + count_non_target_species +
			count_bait_only + count_empty_hooks - count_bent_broken) %>%
		mutate(count_bait_only = replace(count_bait_only, which(count_bait_only == 0), 1)) %>%
		mutate(prop_bait_hooks = count_bait_only / total_hooks) %>%
		mutate(hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks)) %>%
		mutate(expected_catch = round(count_target_species * hook_adjust_factor))

	## find out which events are in hook data (d1) but not set data: (d2)
	d1.not.d2 = setdiff(d1$fishing_event_id, d2$fishing_event_id)

	hook_adjusted_data <- left_join(d2, adjust, by = c("fishing_event_id", "year"))
	## Final output
	outnam = paste0("./data/", spc, "_outsideHBLL_hook_adjusted_", timestamp)
	##saveRDS(hook_adjusted_data,   file = paste0(outnam, ".rds"))  ## largely useless format
	save("hook_adjusted_data",    file = paste0(outnam, ".rda"))
	write.csv(hook_adjusted_data, file = paste0(sub("data","tables",outnam), ".csv"))
	hbll = as.data.frame(hook_adjusted_data); save("hbll",file=paste0("./data/hbll", strSpp,".rda"))  ## for use in calcHBLL
	return(hook_adjusted_data)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getHBLL

## predREBS-----------------------------2026-05-07
## Predict probability of RER (BSR is the inverse)
## ------------------------------------------NF|RH
predREBS <- function(dat, indvar=c("X","Y"), cutoff=5, rmSM=TRUE, 
   col=c("navyblue","blue","lightseagreen","gold","orangered"),
   ncells=250, addpts=FALSE, outnam="predREBS",
   png=FALSE, pngres=400, PIN=c(9,9.5), lang=c("f","e"))
{
	loco = lenv()
	## Remove seamounts
	if (rmSM) {
		zMML = is.element(colnames(dat), c("MAJOR_STAT_AREA_CODE", "MINOR_STAT_AREA_CODE", "LOCALITY_CODE"))
		if (sum(zMML)!=3) {
			message("data file does not contain spatial coordinates (major, minor, locality) to remove seamounts)")
			browser(); return()
		}
		dat$major = dat[,"MAJOR_STAT_AREA_CODE"]
		dat$minor = dat[,"MINOR_STAT_AREA_CODE"]
		dat$locality = dat[,"LOCALITY_CODE"]
		dat$locality[is.na(dat$locality)] = 0  ## zapSeamounts seems to abort if any values are NA
		dat = zapSeamounts(dat)
	}
	## sdmTMB model (2-D space) : make a mesh
	mesh <- make_mesh(dat, xy_cols=indvar, cutoff=cutoff)

	## Fit an sdmTMB model (NF 260507)
	fit_space <- sdmTMB(species ~ 1,
		data = dat,
		mesh = mesh,
		family = binomial(link = "logit"),
		spatial = "on"
	)
	
	## Hypothetical new data, NEEDS TO HAVE COLUMNS X AND Y IN UTM COORDINATES
	#new_dat <- data.frame(X=500, Y=5600)
	#new_dat$P_RER <- predict(fit_space, newdata=new_dat, type="response")$est

	## mesh$loc_xy same coordinates as in fit_space$data
	#new_dat = as.data.frame(mesh$loc_xy[sample(1:nrow(mesh$loc_xy),100),]) ## for debugging
	new_dat = as.data.frame(mesh$loc_xy) ## whole shebang
	poop <- predict(fit_space, newdata=new_dat, type="response")

	## Stick to using akima package
	data("nepacLL", package="PBSmapping", envir=loco)
	attr(nepacLL, "zone") = 9
	nepacUTM = convUL(nepacLL)
	unpackList(poop)
	xlim = range(X); ylim=range(Y)

	interp_data <- interp(x=X, y=Y, z=est, 
		xo = seq(xlim[1], xlim[2], length=ncells),
		yo = seq(ylim[1], ylim[2], length=ncells),
		linear=TRUE, duplicate="mean", extrap=FALSE)
	brush = colorRampPalette(col)

	## Start plotting this mess around
	if (png) createFdir(lang, dir=".")
	fout.e = outnam
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(file=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
		}
		expandGraph(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plotMap(nepacUTM, xlim=xlim, ylim=ylim, col="transparent", border="gainsboro", plt=NULL, cex.axis=1.2, cex.lab=1.5, lwd=0.2)

		image(interp_data, 
			main = "", add=TRUE, xlab = "", ylab = "", 
			col=brush(ncells)) # Add color palette
		if (addpts) {
			points(X, Y, pch=21, cex=0.6, col="black", bg="gainsboro", lwd=0.5)
		}
		addPolys(nepacUTM, col="honeydew", border="grey", lwd=0.2)
		addStrip(0.05, 0.35, col=col, lab=show0(round(seq(min(est),max(est),length=length(col)),2),2), xwidth=0.01, yheight=0.3, adj=-0.2)
		addLabel(0.95, 0.95, txt=linguaFranca("Proportion RER (BSR = 1 -RER)",l), cex=1, adj=c(1,0))
		box()
		if (png) dev.off()
	}; eop()  ## end lang loop
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~predREBS


