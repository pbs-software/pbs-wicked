##==============================================================================
## Module : wicked_survey
## -----------------------
## calcAE...........Calculate ageing error using the function 'RunFn()' from 'AgeingError'
## calcTweedie......Calculate Tweedie CPUE (adapted from Gulf Hake, augmented by Gemini)
## getHBLL..........Get survey series indices for HBLL surveys
## predREBS.........Predict probability of RER (BSR is the inverse)

##-----Supplementary hidden functions-----
##===============================================================================
## Note: not sure whether to rename functions so it's obvious that they are wicked
## See object `wickiverb' for possible guidance (RH 260506)

## calcAE-------------------------------2026-09-22
##  Calculate ageing error using the function 'run()'
##  from the R package 'AgeingError'.
## Note: function 'RunFn()' now deprecated (Ian Taylor, pers.commm.)
##  Adapted from Kendra Holt's code used in Petrale 2024
##  Attempt to de-hadley the code
##  ageDat  : use data from 'gfb_age_precision.sql' (via 'PBStools::plotAgeError.r')
##  atype   : AGE_READING_TYPE_CODE (1=Final, 2=Primary, 3=Precision, 4=Secondary)
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
	)
	#spcode = "AE" ## AE=AgeingError : just keep it the same for all species

	fr = function(x) { eval(parse(text=deparse(x))) }  ## convert text with unicode to accented nonsense

	## KH: From gfsynopsis res doc appendix A, it seems likely that 
	## age_reading_type_code = 2 is the initial primary age, 
	## while age_reading_type_code = 3 is the precision age.

	## require (PBSwicked)
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

	## Exploratory Plots
	## --------------------------------------------
	## Stupid themes for stupid ggplots
	themeAE = theme(
		panel.background = element_rect(fill="white"), 
		panel.border = element_rect(colour="black", fill=NA, size=1),
		axis.title = element_text(size=18), ## Changes both X and Y titles
		axis.text  = element_text(size=12)  ## Changes both X and Y tick labels
	)

	## Distribution of differences between readers over time
	g_DistbyYear <- ggplot(doubleReadDat,aes(x=as.factor(year), y=readerDiff)) + geom_boxplot(fill="chartreuse") + geom_hline(yintercept=0, color = "red") + coord_cartesian(ylim = c(-15, 15))  + ylab("Difference between readers") + xlab("Year")

	## Plot Dist by Year
	fout.e = paste0("calcAE(", spcode, ")-DistbyYear")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e'=paste0("./english/",fout.e), 'f'=paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		g_DistbyYear_out = g_DistbyYear + themeAE
		if (l=="f") g_DistbyYear_out = g_DistbyYear_out + xlab(fr("ann\u{00E9}e")) + ylab(fr("diff\u{00E9}rence entre les techniciens"))
		print(g_DistbyYear_out)
		if (png) dev.off()
	}; eop()

	## Plot Reader Comparison
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
		g_Scatter_out = g_Scatter + themeAE
		if (l=="f") g_Scatter_out = g_Scatter_out + xlab(fr("estimation de l'\u{00E2}	ge primaire")) + ylab(fr("estimation d'\u{00E2}ge pr\u{00E9}cise"))
		print(g_Scatter_out)
		if (png) dev.off()
	}; eop()

	## Fit Linear Regression : precision age as a function of primary age 
	mod.lm <- lm(Prec_Age ~ Primary_Age, doubleReadDat)
	mod_summary <- summary(mod.lm)

	## Try bootstrapping ages for precision readers from primary readers (experimental)
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
	} ## end debug

	## Plot regression residuals
	fout.e = paste0("calcAE(", spcode, ")-lm(prec~prim)-res")
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			createFdir(lang=l)
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		plot(doubleReadDat$Primary_Age, mod_summary$residuals, type="n", cex.axis=1.2, cex.lab=1.5, xlab=linguaFranca("Age by primary reader",l), ylab=linguaFranca("Residuals",l))
		abline(h=0, col="blue")
		zpos = mod_summary$residuals > 0
		zneg = mod_summary$residuals <= 0
		points(doubleReadDat$Primary_Age[zpos], mod_summary$residuals[zpos], pch=24, col="green4", bg="green")
		points(doubleReadDat$Primary_Age[zneg], mod_summary$residuals[zneg], pch=25, col="red", bg="pink")
		if (png) dev.off()
	}; eop()

	createTdir()  ## create directories called 'tables' and 'data'

	residAge <- data.frame(Resids = mod_summary$residuals, Primary_Age = doubleReadDat$Primary_Age)
	residAge <- as_tibble(residAge)

	## Write SD at age to CSV files
	sigAge <- residAge %>% group_by(Primary_Age) %>% summarize(sigAtAge=sd(Resids), n=n())
	write.csv(sigAge, paste0("./tables/calcAE(", spcode, ")-sigAge0.csv") )
	sigAge <- sigAge%>%filter(n>5)
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
		plot(sigAge$Primary_Age, sigAge$sigAtAge, type="l", xlab=switch(l, 'e'="Primary Age", 'f'=fr("\u{00E2}ge primaire")), ylab=switch(l, 'e'="Residual SD at age", 'f'=fr("\u{00C9}T r\u{00E9}siduel \u{00E0} l'\u{00E2}ge")) )
		points(sigAge$Primary_Age, sigAge$sigAtAge, pch=21, cex=1.2, col="black", bg="gold")
		if (png) dev.off()
	}; eop()

	## NWFSC Package AgeingError (Schnute/Punt AE algorithm)
	## -----------------------------------------------------
	#require(tidyverse)  ## ain't no one got time for this
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
	}
	## This is where all runs will be located
	AEdir    <- "./data"
	Nreaders <- 2 ## primary and precision
	dat      <- doubleReadDat
	dat      <- dat %>% select(Primary_Age, Prec_Age)
	AgeReads <- as.matrix(dat)

	## PJS wants to see 
	dat4pjs  <- doubleReadDat %>% select(specimen_id, year, Primary_Age, Prec_Age)
	dat4pjs$precDiff  <- dat4pjs$Prec_Age - dat4pjs$Primary_Age
	minPrecAnom <- 5
	dat4slc <- dat4pjs[abs(dat4pjs$precDiff) >= minPrecAnom,]
#browser();return()

	## Format data
	Nreaders = ncol(AgeReads)
	AgeReads = ifelse(is.na(AgeReads), -999, AgeReads)  ## Change NA to -999 (which the Punt software considers missing data)
	## Potentially eliminate rows that are only read once
	## These rows have no information about reading error, but are potentially informative about latent age-structure
	## It is unknown whether eliminating these rows degrades estimation of error and bias, and is currently recommended to speed up computation
	#KeepRow = ifelse(rowSums(ifelse(AgeReads==-999,0,1),na.rm=TRUE)<=1,FALSE,TRUE)
	#AgeReads = AgeReads[KeepRow,]

	## Testing exclusion of ages > xx (not sure why we would want to do this)
	AgeReads  <- as.data.frame(AgeReads)  ## looks like this will be the input for the new 'run' function in 'AgeingError'
	AgeReads1 <- AgeReads[AgeReads$Primary_Age <= maxage*2, ]  ## arbitrarily assign a high cutoff
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
	SigOpt <- c(ifelse(strSpp%in%c("405"),2,3), -1)

	## Define minimum and maximum ages for integral across unobserved ages
	MinAge   = 1
	MaxAge   = max(AgeReads)
	KnotAges = list(NA, NA)
	ReaderNames = c("Primary reader", "Precision reader")

	write.csv(AgeReads,  paste0("./tables/calcAE(", spcode, ")-AgeReads.csv") )   ## (RH 251020)
	write.csv(AgeReads1, paste0("./tables/calcAE(", spcode, ")-AgeReads1.csv") )  ## (RH 260515)
	write.csv(AgeReads2, paste0("./tables/calcAE(", spcode, ")-AgeReads2.csv") )  ## (RH 251020)

	if (useTMB) {
		## Ian Taylor suggests we use 'run' instead of 'RunFn', which is deprecated
		write_files(dat=AgeReads, dir=AEdir, biasopt=BiasOpt, sigopt=SigOpt)
		modRes <- try(run(directory=AEdir))
		if(inherits(modRes, "try-error") ) {
		message("------------------------\n'run' failed for some reason") ; print(modRes); setwd(cwd); browser(); return() }
		## run() generates png files  automatically

		## Hacked into AgeingErorr to figure out how to use the function 'plot_output'
		## Note that both run() and ProcessResults() will produce these plots.
		load(file.path(AEdir, "AgeingError.lda"))  ## contains SaveAll
		for (l in lang) {
			## french is only partly there (e.g., ReaderNames and delimiters) but will have to do for now
			changeLangOpts(L=l)
			switch(l, 
				'e' = { SaveDir="./english/"; ReaderNames=c("Primary reader", "Precision reader") },
				'f' = { SaveDir="./french/";  ReaderNames=c("technicien principal", eval(parse(text=deparse("technicien de pr\u{00E9}cision"))) ) }
			)
			plot_output(Data=AgeReads2, IDataSet=1, MaxAge=SaveAll$data$MaxAge, Report=SaveAll$report, 
				Nparameters=length(SaveAll$parameters$SDPar), LogLike=SaveAll$report$Obj_fun, 
				subplot=1:3, ReaderNames=ReaderNames, Species="AgeingError", SaveDir=SaveDir)
		} ; eop()
		## see estimated parameters
		modPar = modRes$model$par
		## see model selection results
		modSel = modRes$output$ModelSelection
		## see ageing error matrices
		modErr = modRes$output$ErrorAndBiasArray
		Age    = 0:maxage
		SD     = modErr["SD",paste("Age",Age),"Reader 1"]
	} else {  ## old executable (deprecated, no longer tested, last used for SGR 2025)
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
		errorEsts <- output$ErrorAndBiasArray[,,1]  ## Reader 1
		print(errorEsts)
		write.csv(errorEsts, paste0("./tables/calcAE(", spcode, ")-ageErrorArray_model.csv") )

		Age  <- 0:maxage
		SD   <- errorEsts["SD", paste("Age",Age)]

		## The code from here down to 'end old executable' is deprecated
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#lfit <- lm(SD ~ Age)  ## not really useful for a curvilinear 
		## Fit the linear model to log-transformed data
		#efit_lm <- lm(log(SD) ~ Age)
		#b_lm    <- coef(efit_lm)["Age"]
		#a_lm    <- exp(coef(efit_lm)["(Intercept)"])
		#cv_lm   <- summary(efit_lm)$sigma / mean(log(SD)) * 100 
		## Fit the non-linear model
		#SD_add   <- a_lm * exp(b_lm * Age) + rnorm(length(Age), sd = summary(efit_lm)$sigma)
		#efit_nls <- nls(SD_add ~ a * exp(b * Age), start = list(a = 1, b = 0.1))
		#a_nls    <- coef(efit_nls)["a"]
		#b_nls    <- coef(efit_nls)["b"]
		#cv_nls   <- summary(efit_nls)$sigma / mean(log(SD)) * 100 
	} ## end old executable

	## Use observed precision ages by primary age (for old and new methods)
	age.obs = split(doubleReadDat$Prec_Age, doubleReadDat$Primary_Age)
	sd.obs  = sapply(age.obs,sd)
	use.obs = is.element(as.numeric(names(sd.obs)), 0:maxage)  ## restricts data to less than that used by AgeigError
	use.obs = rep(TRUE, length(age.obs))  ## use all data in linear model
	sd.use  = sd.obs[use.obs & !is.na(sd.obs) & sd.obs > 0]  ## (RH 260922) expanded the criteria for useable data in lm fit
	#sd.use[is.na(sd.use)] = GT0(0)  ## complete rubbish
	age.use = as.numeric(names(sd.use))
	## Fit the linear model to log-transformed data
#browser();return()
	fit.lm  = lm(log(sd.use) ~ age.use)


	## Cleanup the mess left behind from AgeingError
	## Remove png files:
	pngs = normalizePath(list.files(path=AEdir, pattern="\\.png$", full.names=TRUE), winslash="/")
	rubbish = file.remove(pngs)
	## Move csv files to './tables' (we'll have to assume that AEdir is relative to the cwd)
	for (idir in c("./french", "./english", AEdir)) {
		csvs = list.files(path=idir, pattern="\\.csv$")
		if (length(csvs)>0) {
			CSVdir  = sub("\\./.*", "./tables", idir)
			rubbish = file.copy(from=file.path(idir,csvs), to=file.path(CSVdir,csvs), overwrite=TRUE, copy.date=TRUE)
			rubbish = file.remove(file.path(idir,csvs))
		}
	}

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
		ylim = c(0, max(c(SD, sd.use), na.rm=T))
		plot(Age, SD, pch=19, ylim=ylim, xlab=linguaFranca("Age",l), ylab=linguaFranca("Standard deviation",l), cex.axis=1.25, cex.lab=1.5, col="gainsboro")
		abline(a=0, b=0.1, lwd=2, lty=2, col="red")  ## CASAL CV=10%
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
				legtxt = switch(l, 'e'="Approx. SD from (min, max)", 'f'=fr("\u{00C9}T approximatif (min, max)"))
				legcol="purple"; leglty=NA; legpch=17; legbg=NA
		} ## end debug
		## add SD of observed precision ages by primary age
		points(age.use, sd.use, pch=16, col="deepskyblue")
		## add lm fit to previous points
		lines(age.use, exp(predict(fit.lm, newdata=data.frame('age.use'=age.use))), col="blue", lty=5, lwd=2) # transform predictions back to original scale
		legtxt = switch(l,
			'e'={c(legtxt, "SD of precision ages at primary age",
				"linear model fit to log SD", 
				paste0("AgeingError SD using functional form ", SigOpt[1]),
				"Constant CV = 0.10")},
			'f'={c(legtxt, fr("\u{00C9}T des \u{00E2}ges de pr\u{00E9}cision \u{00E0} l'\u{00E2}ge primaire"), 
				fr("Mod\u{00E8}le lin\u{00E9}aire ajust\u{00E9} au logarithme de l'\u{00C9}T"), 
				fr(paste0("AgeingError \u{00C9}T utilisant la forme fonctionnelle ", SigOpt[1])),
				"CV constant = 0,10")}
		)
		legcol = c(legcol, "deepskyblue","blue","black","red")
		leglty = c(leglty, NA,5,1,2)
		legpch = c(legpch, 16,NA,21,NA)
		legbg  = c(legbg, NA,NA,"green",NA)

		## Add the main AgeingError fit
		lines(Age, SD, col="black", lwd=2)
		points(Age, SD, pch=21, col="black", bg="green", cex=1.2)
		addLegend(0.025, 0.975, legend=legtxt, pch=legpch, pt.bg=legbg, lty=leglty, col=legcol, pt.cex=1.2, pt.lwd=1, lwd=2, seg.len=3, cex=1.2, bty="n", xjust=0)
#browser();return()
		if (png) dev.off()
	}; eop() ## end lang loop
	if (debug) { browser(); return() }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcAE


## calcTweedie -------------------------2026-10-21
##  Calculate Tweedie CPUE (adapted from Gulf Hake)
##  Substantial help from Google Gemini (GG)
## -----------------------------------PJS|SA|GG|RH
calcTweedie <- function(strSpp="417", gear="bottom trawl", 
   facs=c("location", "depth", "major", "vessel", "latitude"),
   covs=NULL, rand=NULL, ylim=NULL, refit=TRUE,
   yrs=1996:2025, major=3:9, vess.yrs=8, vess.trips=7, 
   drange=c(75,400), dbin=25, effort=6, 
   manualExtract=FALSE, rawIndex=FALSE, ex0=FALSE,
   year.type="CY", gulf=FALSE, pjs.index=NULL, showCatch=FALSE,
   png=FALSE, pngres=400, PIN=c(9,6), outnam="tweedie.cpue", lang=c("f","e"))
{
	ipak = strSpp=="225" && gulf  ## logical : TRUE if inshore (gulf) hake
	## Gemini's reworking of the original calcFY()
	calcFY <- function(dt, 
		startM = if(ipak) 11 else 4, 
		startD = if(ipak) 15 else 1, 
		endM   = if(ipak) 4 else 3, 
		endD   = if(ipak) 15 else 31, 
		nextYr = ipak)
	{
		## Format vector cleanly using base R sprintf instead of a custom pad0
		yrmoda <- substring(dt, 1, 10)
		yr     <- as.numeric(substring(yrmoda, 1, 4))
		moda   <- substring(yrmoda, 6, 10)

		## Construct standardized MM-DD comparison strings
		startMD <- sprintf("%02d-%02d", as.numeric(startM), as.numeric(startD))
		endMD   <- sprintf("%02d-%02d", as.numeric(endM), as.numeric(endD))

		## Initialize the output vector directly from the calendar years
		fyr <- yr

		## 1. Handle missing data safely up front
		is_bad <- is.na(moda) | nchar(yrmoda) < 10

		## 2. Vectorized Boolean subsetting arrays
		## Check if the fiscal year crosses a calendar boundary (e.g., Nov to Apr)
		crosses_year <- startMD > endMD
		if (crosses_year) {
			## Case A: Fiscal year spans across Dec 31 / Jan 1 boundary (e.g., Nov 15 - Apr 15)
			## A date belongs to the 'Late Segment' (Nov 15 - Dec 31) or 'Early Segment' (Jan 1 - Apr 15)
			is_late_segment  <- (moda >= startMD) & !is_bad
			is_early_segment <- (moda <= endMD) & !is_bad
			if (nextYr) {
				# If nextYr = TRUE, the fiscal year is named after the calendar year it ends in
				fyr[is_late_segment]  <- yr[is_late_segment] + 1
				fyr[is_early_segment] <- yr[is_early_segment]
			} else {
				# If nextYr = FALSE, the fiscal year is named after the calendar year it starts in
				fyr[is_late_segment]  <- yr[is_late_segment]
				fyr[is_early_segment] <- yr[is_early_segment] - 1
			}
		} else {
			## Case B: Fiscal year stays entirely inside one calendar year (e.g., Apr 1 - Mar 31)
			is_inside_fy <- (moda >= startMD & moda <= endMD) & !is_bad
			if (nextYr) {
				fyr[is_inside_fy]  <- yr[is_inside_fy] + 1
				fyr[!is_inside_fy] <- yr[!is_inside_fy]
			} else {
				fyr[is_inside_fy]  <- yr[is_inside_fy]
				fyr[!is_inside_fy] <- yr[!is_inside_fy] - 1
			}
		}
		## Force malformed inputs to NA instead of 0
		fyr[is_bad] <- NA
		names(fyr) <- yrmoda
		return(fyr)
	}
	## Main function starts here-------------------
	loco = lenv()
	data("species", package="PBSdata", envir=loco)
	spp.code  = tolower(species[strSpp,"code3"])
	gear.code = paste0(sapply(strsplit(gear, split=" +"), function(x){paste0(substring(x,1,1),collapse="")}),collapse="+")
	fnam = paste0("./data/",spp.code, ".", gear.code, ".rda")
	dnam = paste(spp.code, gear.code, sep=".")
	outnam = paste0(outnam,".",spp.code)

	if (!file.exists(fnam)){
		#wwr.bt = get_cpue_index(gear="bottom trawl", min_cpue_year = 1996, major=pad0(3:9,2)) ## SA gfdata function doesn' allow multiple majors
		getData("get-cpue-index.sql",dbName="GFFOS",strSpp="417",major=3:9, gear=paste0("'", toupper(gear), "'"))
		colnames(PBSdat) = tolower(colnames(PBSdat))
		## use convert 20 or 120 for ODBC canonical yyyy-mm-dd hh:mi:ss (24h)
		## then in R : as.POSIXct(strptime(junk$FE_START_DATE, "%Y-%m-%d %H:%M:%S"))
		PBSdat$best_date     = as.POSIXct(strptime(PBSdat$best_date, "%Y-%m-%d %H:%M:%S"))
		PBSdat$fe_start_date = as.POSIXct(strptime(PBSdat$fe_start_date, "%Y-%m-%d %H:%M:%S"))
		PBSdat$fe_end_date   = as.POSIXct(strptime(PBSdat$fe_end_date, "%Y-%m-%d %H:%M:%S"))
		assign(dnam,  as_tibble(PBSdat), envir=loco )  ## only to keep format that SA uses
		save(list=dnam, file=fnam)
	} else {
		load(fnam)
	}
	eval(parse(text=paste0("spp.dat <- ", dnam)))
	class(spp.dat) = "data.frame"
	spp.cat = spp.dat  ## spp.cat for use later
#browser();return()

	if (!refit && "saved_fits" %in% tlisp()) {
		unpackList(ttcall(saved_fits))
	} else {
		## Calendar/Fishing year (main response, needs to be called 'year_factor for 'glmmTMB::glmmTMB')
		if (year.type=="FY") {
			spp.dat$fyear = calcFY(spp.dat$best_date)
			spp.dat = spp.dat[is.element(spp.dat$fyear, yrs),]
			#spp.dat = spp.dat[spp.dat$fyear>0 & !is.na(spp.dat$fyear),]
			.flush.cat("filter by 'fyear' = ", nrow(spp.dat), " records\n")
			if (strSpp=="225" && gulf) { ## gulf hake
				## fishing years from 15 November 2005 to 15 April 2013 dropped because hake catches from the defined area and period were minimal (PJS, BM)
				spp.dat = spp.dat[!is.element(spp.dat$fyear, c(2007:2013, 2023)),]  ## years for Jan-Apr period (inlc. Nov-Dec)
				#spp.dat = spp.dat[!is.element(spp.dat$fyear, c(2023)),]  ## years for Jan-Apr period (inlc. Nov-Dec)
				.flush.cat("filter by dropping no-hake 'fyear' = ", nrow(spp.dat), " records\n")
			}
			spp.dat$year_factor = as.factor(spp.dat$fyear)
		} else {
			spp.dat = spp.dat[is.element(spp.dat$year, yrs),]
			.flush.cat("filter by 'cyear' = ", nrow(spp.dat), " records\n")
			spp.dat$year_factor = as.factor(spp.dat$year)
		}
		## Start qualifying data
		## ---------------------
		spp.dat = spp.dat[is.element(spp.dat$species_code, strSpp),]  ## redundant but WTH
		.flush.cat("filter by 'species_code' = ", nrow(spp.dat), " records\n")
	
		## PJS filter 1 : years (redundant but just in case)
		spp.dat = spp.dat[is.element(spp.dat$year_factor, yrs),]
		.flush.cat("filter by 'year_factor' = ", nrow(spp.dat), " records\n")
	
		## PJS filter 2 : major areas
		if (gulf) {
			major = 1
			## Area 12 and Area 20 have been excluded from this list because of the potential for hake catches in these areas to include 'outside' hake.
			minor.areas = c(13, 14, 15, 16, 17, 18, 19, 29)
			names(minor.areas) = c("Johnstone Strait", "Qualicum", "Powell River", "Sechelt", "Porlier Pass", "Satellite Channel", "Victoria", "Fraser River")
			spp.dat = spp.dat[is.element(spp.dat$minor_stat_area_code, minor.areas),]
			.flush.cat("filter by inside 'minor area' = ", nrow(spp.dat), " records\n")
		} else {
			spp.dat = spp.dat[is.element(spp.dat$major_stat_area_code, major),]
			.flush.cat("filter by major PMFC area = ", nrow(spp.dat), " records\n")
		}
	
		## PJS filter 3 : vessels
		## Vessel participation (using Gemini) use vectorized version:
		spp.dat$vessel_id = spp.dat$vessel_registration_number ## (VRN)
		## 0. Filter for records where species is Widow Rockfish and catch is greater than 0
		target_catch <- subset(spp.dat, species_code == strSpp & catch_kg > 0)
		## 1. Create a matrix of trip counts by vessel (rows) and year (columns)
		vessel_year_matrix <- table(target_catch$vessel_id, target_catch$year_factor)
		## 2. Convert counts to logical matrix: TRUE if >= 7 trips, otherwise FALSE
		is_active_matrix <- vessel_year_matrix >= vess.yrs
		## 3. Sum active years per row and subset vessels with >= 7 active years
		active_year_counts <- rowSums(is_active_matrix)
		target_vessels_alt <- names(active_year_counts[active_year_counts >= vess.trips])
		spp.dat = spp.dat[is.element(spp.dat$vessel_id, target_vessels_alt),]
		.flush.cat("filter by 'VRN' = ", nrow(spp.dat), " records\n")
	#browser();return()
	
		## PJS filter 4 : depths
		spp.dat = spp.dat[!is.na(spp.dat$best_depth),]
		spp.dat = spp.dat[spp.dat$best_depth > drange[1] & spp.dat$best_depth <= drange[2],]
		.flush.cat("filter by 'best_depth' = ", nrow(spp.dat), " records\n")
	
		## PJS filter 5 : gear
		spp.dat = spp.dat[is.element(spp.dat$gear, toupper(gear)),]
		.flush.cat("filter by 'gear' = ", nrow(spp.dat), " records\n")
	
		## PJS filter 6 : latitude (and longitude)
		spp.dat = spp.dat[spp.dat$latitude>=48.3  & spp.dat$latitude<=54.6  & !is.na(spp.dat$latitude),]
		spp.dat = spp.dat[spp.dat$longitude>=-134 & spp.dat$longitude<=-122 & !is.na(spp.dat$longitude),]
		.flush.cat("filter by 'latitude-longitude' = ", nrow(spp.dat), " records\n")
		spp.dat$X = spp.dat$longitude
		spp.dat$Y = spp.dat$latitude
	
		## PJS filter 7 : effort
		## Valid estimate of time towed that was > 0 hours and <= x hours.
		spp.dat = spp.dat[!is.na(spp.dat$fe_start_date) & !is.na(spp.dat$fe_end_date),]
		.flush.cat("filter by 'start-end dates' = ", nrow(spp.dat), " records\n")
		spp.dat$effort = round(as.numeric(difftime(spp.dat$fe_end_date, spp.dat$fe_start_date, units = "hours")),3)
		spp.dat = spp.dat[spp.dat$effort>0 & spp.dat$effort<=effort & !is.na(spp.dat$effort),]
		.flush.cat("filter by 'effort' = ", nrow(spp.dat), " records\n")
	
		## PJS filter 8 : catch, but only if excluding zero-catch
		## Catch (t) and CPUE
		spp.dat$catch_t = spp.dat$catch_kg / 1000.
		spp.dat$cpue    = spp.dat$catch_kg / spp.dat$effort
		if (ex0) {
			## exclude zero catches
			outnam = paste0(outnam,".ex0")
			spp.dat = spp.dat[spp.dat$catch_kg>0 & !is.na(spp.dat$catch_kg),]
			.flush.cat("filter by 'excluding zero-catch' = ", nrow(spp.dat), " records\n")
		}
		if (gulf) {
			## Annual catch in Gulf region
			spp.cat = spp.cat[is.element(spp.cat$minor_stat_area_code, minor.areas),]
			spp.cat$year = as.numeric(substring(spp.cat$best_date,1,4))
			ycat = sapply(split(spp.cat$catch_kg, spp.cat$year),sum)/1000. ## tonnes
		}
		save ("spp.dat", file="./data/spp.dat.rda")
	
		## Collect user-specified factors (apply the factor after qualification)
		## ------------------------------
		if ("location" %in% facs) {
			spp.dat$location = paste(spp.dat$major_stat_area_code, pad0(spp.dat$minor_stat_area_code,2), pad0(spp.dat$locality_code,2), sep=".")
			loc.cat = split(spp.dat$catch_kg, spp.dat$location)
			loc.good = sapply(loc.cat,countVec)>=100
			loc.bad  = !loc.good
			spp.dat$location[is.element(spp.dat$location, names(loc.bad)[loc.bad])] = 'plus'
			spp.dat$location_fac = as.factor(spp.dat$location)
		}
		if ("depth" %in% facs) {
			spp.dat$depth_fac = as.factor(ceiling(spp.dat$best_depth / dbin) * dbin)
		}
		if ("major" %in% facs) {
			spp.dat$major_fac = as.factor(spp.dat$major_stat_area_code)
		}
		if ("vessel" %in% facs) {
			spp.dat$vessel_fac = as.factor(spp.dat$vessel_id)
		}
		if ("latitude" %in% facs) {
			spp.dat$latitude = ceiling(spp.dat$Y / 0.1) * 0.1
			lat.cat = split(spp.dat$catch_kg, spp.dat$latitude)
			lat.good = sapply(lat.cat,countVec)>=100
			lat.bad  = !lat.good
			spp.dat$latitude[is.element(spp.dat$latitude, names(lat.bad)[lat.bad])] = 99
			spp.dat$latitude_fac = as.factor(spp.dat$latitude)
		}
		if ("month" %in% facs) {
			spp.dat$month = as.numeric(substring(spp.dat$best_date,6,7))
			spp.dat$month[is.element(spp.dat$month,12)] = 11
			spp.dat$month_fac = as.factor(spp.dat$month)
		}
		## Collect user-specified continuous covariates
		## --------------------------------------------
		if ("effort" %in% covs) {  ## not a factor but a continuous covariate
			spp.dat$hours_fished = spp.dat$effort
			covs = sub("effort", "hours_fished", covs)
		}
		if ("depth" %in% covs) {   ## not a factor but a continuous covariate
			spp.dat$depth_numeric = spp.dat$best_depth
			## 1. Pre-calculate the orthogonal polynomial components in your dataset
			depth_poly <- poly(scale(spp.dat$depth_numeric), 2)
			spp.dat$depth_p1 <- depth_poly[, 1]
			spp.dat$depth_p2 <- depth_poly[, 2]
			covs = c(covs, "depth_p1", "depth_p2")
			covs = setdiff(covs, "depth")
		}
		## Collect user-specified random effects
		## -------------------------------------
		if ("vessel" %in% rand) {  ## not a factor but a random effect
			spp.dat$vessel_raneff = spp.dat$vessel_id
			rand = sub("vessel", "vessel_raneff", rand)
		}

		## Superfluous table but keep for general information
		junk = split(spp.dat,spp.dat$year_factor)
		spp.tab = sapply(junk,function(x){
			xout  = c(
				nvess = length(unique(x$vessel_id)),
				ntid  = length(unique(x$trip_id)),
				#nfeid = length(unique(x$fishing_event_id)),  ## early years, FEIDs were not unique (e.g., many FEID=1)
				nex0  = countVec(x$catch_kg),
				nall  = length(x$catch_kg),
				pct0  = 100*(length(x$catch_kg) - countVec(x$catch_kg)) / length(x$catch_kg),
				tcat  = sum(x$catch_t),  ## express in tonnes
				teff  = sum(x$effort[x$catch_kg>0]),
				cpue  = sum(x$catch_kg[x$catch_kg>0]) / sum(x$effort[x$catch_kg>0])
			)
			return(xout)
		})
		spp.tab = as.data.frame(t(spp.tab))
		if (gulf) {
			spp.tab$gcat = ycat[intersect(rownames(spp.tab), names(ycat))] ## add in annual Gulf catch
		}
		write.csv(spp.tab,file="./tables/spp.tab.csv",row.names=TRUE)

		predictors = c("year_factor", paste0(facs, "_fac"), covs, rand)
		spp.tdf = as_tibble(spp.dat[,c("catch_kg", predictors)])
	
		## Gemini's R code for raw manual extraction
		## -----------------------------------------
		## 1. Fit the model without an intercept (0 +) to get absolute year effects
		mess = paste0("spp.fit <- glmmTMB::glmmTMB(catch_kg ~ 0 + ", paste0(predictors, collapse=" + "), ", data=spp.tdf, family=glmmTMB::tweedie(link=\"log\"), 	map=list(psi=factor(NA)), control=glmmTMB::glmmTMBControl(profile=TRUE, rank_check=\"adjust\", optCtrl=list(iter.max=1000, eval.max=1000)) )")
		mess = sub("hours_fished", "offset(log(hours_fished))", mess)
		mess = sub("depth_numeric", "poly(scale(depth_numeric), 2)", mess)
		mess = sub("vessel_raneff", "(1 | vessel_raneff)", mess)      ## <-- Vessel treated as a random intercept
		model.call = mess
		eval(parse(text=model.call))

		## 2. Pull the fixed effects conditional matrix cleanly
		coef_matrix <- summary(spp.fit)$coefficients$cond

		## 3. Locate the rows corresponding to your years
		year_rows <- grep("year_factor", rownames(coef_matrix))

		## 4. Force numeric vectors explicitly to safeguard mathematical operators
		est_log   <- as.numeric(coef_matrix[year_rows, "Estimate"])
		se_log    <- as.numeric(coef_matrix[year_rows, "Std. Error"])
		years_vec <- as.numeric(gsub("year_factor", "", rownames(coef_matrix)[year_rows]))

		## 5. Safely apply exp() transformations to generate absolute and centered bounds
		## Center the log estimates by subtracting their mean
		est_log_centered <- est_log - mean(est_log)

		## Collect results from manual extraction
		spp.res <- data.frame(
			year          = years_vec,
			est_link      = est_log,
			se_link       = se_log,

			## Original absolute scales
			est           = exp(est_log),
			lwr           = exp(est_log - 1.96 * se_log),
			upr           = exp(est_log + 1.96 * se_log),

			# Centered relative scales (geometric mean = 1)
			est_centered  = exp(est_log_centered),
			lwr_centered  = exp(est_log_centered - 1.96 * se_log),
			upr_centered  = exp(est_log_centered + 1.96 * se_log)
		)

		# Sort chronologically 
		spp.res <- spp.res[order(spp.res$year), ]
		write.csv(spp.res, file="./tables/spp.res.csv",row.names=TRUE)

		## --------------------------------------------------------------------------
		## 'emmeans' derivation is significantly better and more rigorous for generating
		##  an abundance index than the raw manual extraction method (above).
		## --------------------------------------------------------------------------

		## emmeans index calculation (Gemini 261016)
		tweedie_index <- emmeans::emmeans(
			spp.fit, 
			specs = ~ year_factor, 
			nuisance = paste0(facs, "_fac"), # Averages across space and time efficiently
			type = "response"
		)
		# Convert to a dataframe to look at your index
		tweedie.index <- as.data.frame(tweedie_index)

		# 1. Calculate the overall mean of the standardized responses
		grand_mean <- mean(tweedie.index$response)

		# 2. Divide the estimates and intervals by that same grand mean
		tweedie.index$index_centered <- tweedie.index$response  / grand_mean
		tweedie.index$LCL_centered   <- tweedie.index$asymp.LCL / grand_mean
		tweedie.index$UCL_centered   <- tweedie.index$asymp.UCL / grand_mean

		# 3. Calculate the log-scale standard error for SS3
		# (Note: This is scale-invariant, so using response or index_centered yields identical results)
		tweedie.index$ss3_log_se <- sqrt(log(1 + (tweedie.index$SE / tweedie.index$response)^2))

		write.csv(tweedie.index, file="./tables/tweedie.index.csv", row.names=TRUE)

		saved_fits = list(spp.tdf=spp.tdf, spp.tab=spp.tab, spp.fit=spp.fit, spp.res=spp.res, tweedie.index=tweedie.index, model.call=model.call)  ## in case user wants to bypass refitting
		ttput(saved_fits)
	} ## end refit

	## ==============================================================
	## 6. RECREATE THE PLOT IN BASE R (with shaded ribbon boundaries)
	## ==============================================================

	if (manualExtract) {
		model_used = spp.res
	} else {
		model_used = tweedie.index
		cnams = colnames(model_used)
		cnams = 
			sub("UCL_centered",   "upr_centered",
			sub("LCL_centered",   "lwr_centered",
			sub("index_centered", "est_centered",
			sub("asymp.UCL",     "upr",
			sub("asymp.LCL",     "lwr",
			sub("response",       "est", 
			sub("year_factor",    "year", cnams)))))))
		colnames(model_used) = cnams
	}
	## Establish dynamic canvas boundaries
	if (rawIndex) {
		model_est = model_used$est
		model_lwr = model_used$lwr
		model_upr = model_used$upr
	} else {
		model_est = model_used$est_centered
		model_lwr = model_used$lwr_centered
		model_upr = model_used$upr_centered
		outnam = paste0(outnam,".centered")
	}
	y_max <- max(model_upr, na.rm=TRUE) * 1.05
	x_lim <- range(spp.res$year)

	xall = x_lim[1]:x_lim[2]
	xy   = array(NA, dim=c(length(xall),5), dimnames=list(year=xall, index=c("est","lwr","upr","cat","eff")))
	for (j in 1:3) {
		jj   = c("est","lwr","upr")[j]
		jval = switch(jj, 'est'=model_est, 'lwr'=model_lwr, 'upr'=model_upr)
		names(jval) = as.numeric(as.character(model_used$year))
		ii   = as.character(intersect(names(jval),xall))
		xy[ii,jj] = jval[ii]
	}
	yrs.common = intersect(rownames(xy), rownames(spp.tab))
	xy[yrs.common, "cat"] = spp.tab[yrs.common,"tcat"]  ## total annual catch (tonnes)
	xy[yrs.common, "eff"] = spp.tab[yrs.common,"teff"]  ## total annual effort (hours)

	x = as.numeric(rownames(xy)); xlim = range(x)
	if (is.null(ylim))
		ylim = c(0, y_max)

	createFdir(lang) ## create 'english' and 'french' subdirectories for png figures
	fout.e = outnam
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(filename=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
		}
		## Initiate empty frame canvas
		expandGraph(mar=c(3.2,3.2,1,1), mgp=c(2,0.5,0))
		plot(spp.res$year, model_est, type="n", ylim=ylim, xlim=x_lim,
			xlab=linguaFranca(paste0( ifelse(year.type=="CY", "Calendar", "Fishing")," Year"),l), ylab=linguaFranca("Abundance Index",l),
			panel.first = grid(nx=NULL, ny=NULL, col = "lightgray", lty = "dotted"), cex.axis=1.2, cex.lab=1.5, las=1)
		if (!rawIndex) {
			abline(h=1, col="grey")
		}
		## Add PJS delta-lognormal CPUE (if available)
		if (!is.null(pjs.index)) {
			pjs = pjs.index
			x.pjs = as.vector(rbind(pjs$year, pjs$year, rep(NA, length(pjs$year))))
			y.pjs = as.vector(rbind(pjs$lwr, pjs$upr, rep(NA, length(pjs$year))))
			lines(x.pjs, y.pjs, col="black", lwd=2)
			points(pjs$year, pjs$est, pch=22, bg="red", col="black", cex=1.2, lwd=0.75)
		}
		## Render the semi-transparent confidence interval ribbon
		## Note: rev() traces back along the top edge to close the shape path loop cleanly
		if (length(xall) != nrow(model_used)) { ## then years are not continous (may need debugging)
			## Need to deconstruct broken time series because 'polygon' cannot interpret it correctly (RH 221212)
			xpol = c(x,rev(x)); ypol = c(xy[,"lwr"],rev(xy[,"upr"]))
			xyz  = as.data.frame(cbind(xpol,ypol,zpol=!is.na(ypol)))
			## https://stackoverflow.com/questions/55606323/create-counter-for-runs-of-true-among-false-and-na-by-group
			xyz$group = as.numeric(as.factor(cumsum(is.na(xyz$zpol^NA)) + xyz$zpol^NA))
			ngroup = max(length(table(xyz$group))-1, 1)
			for (i in ngroup:1) {
				if (i==0) next  ## shouldn't be needed but just in case
				## start with central group
				ii = is.element(xyz$group, i)
				ixpol = xyz$xpol[ii]; iypol = xyz$ypol[ii]
				## gather disjointed groups
				if (i < ngroup) {
					iii = is.element(xyz$group, i + ngroup)
					ixpol = c(ixpol,xyz$xpol[iii]); iypol = c(iypol,xyz$ypol[iii])
				}
				if (any(iypol > par()$usr[4])){
					polys = data.frame(PID=rep(1,length(iypol)), POS=1:length(iypol), X=ixpol, Y=iypol)
					cpoly = clipPolys(polys, xlim=par()$usr[1:2], ylim=par()$usr[3:4])
					ixpol = cpoly$X
					iypol = cpoly$Y
				}
				polygon(ixpol, iypol, col=lucent("blue",0.10), border=NA)
			}
			lines(x, xy[,"lwr"], lwd=2, lty=3, col="blue")
			lines(x, xy[,"upr"], lwd=2, lty=3, col="blue")
			lines(x, xy[,"est"], col="navy", lwd=2)
			points(x, xy[,"est"], col="navy", bg="yellow", pch=24, cex=1)
		} else {
			## Use Gemini's code for continuous years
			polygon(
				x = c(spp.res$year, rev(spp.res$year)), 
				y = c(model_lwr, rev(model_upr)), 
				col = rgb(0.1, 0.2, 0.8, 0.15),  # Slate blue with 15% opacity
				border = NA
			)
			## Draw upper and lower standard error bounds (matching your plot's blue dotted lines)
			lines(spp.res$year, model_lwr, col = "blue", lty = 3, lwd = 1)
			lines(spp.res$year, model_upr, col = "blue", lty = 3, lwd = 1)
			## Superimpose the primary trend line
			lines(spp.res$year, model_est, col = "darkblue", lwd = 3)
			## Overlay the exact marker style from your plot (yellow-filled triangles)
			points(spp.res$year, model_est, pch=24, bg="yellow", col="darkblue", cex = 1.2)
		}
		addLegend(0.05, 0.95, legend=linguaFranca(c("Tweedie CPUE index", "Delta-lognormal CPUE index"),l), pch=c(24,22), pt.cex=1.5, pt.bg=c("yellow","red"), col=c("darkblue","black"), pt.lwd=1.2, lty=c(1,1), lwd=1.5, seg.len=4, bty="n")
		box()
		if (png) dev.off()
	}; eop()

	if (showCatch) {  ## leftover from Gulf Hake (has not neen tested here)
		## Compare standardised index with catch
		z = !is.na(xy[,"est"])
		idx = xy[z,"est"] / calcGM(xy[z,"est"])
		xy[names(idx),"idx"] = idx
	
		ycat = ycat / calcGM(ycat)
		#ycat = spp.tab$tcat / calcGM(spp.tab$tcat)
		#names(ycat) = rownames(spp.tab)
	#browser();return()
		zz = intersect(names(ycat),rownames(xy))
		xy[zz,"cat"] = ycat[zz]
	
		ylim= c(0, max(xy[,c("idx","cat")], na.rm=T))
		if (png)
			png(filename=paste0("idx.vs.cat.png"), width=PIN[1], height=PIN[2]*2/3, units="in", res=pngres)
		expandGraph(mar=c(3,3,1,1), mgp=c(1.75,0.5,0))
		plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="Fishing Year End (Nov15-Apr15)", ylab="Relative Index", cex.axis=1, cex.lab=1.2, las=1)
		abline(h=1, col="gainsboro")
		lines(x, xy[,"cat"], col="red", lwd=2, lty=5)
		lines(x, xy[,"idx"], col="black", lwd=2)
		addLegend(0.6, 0.95, col=c("red","black"), lty=c(5,1), lwd=2, legend=c("Std. Annual Catch","Std. Fyear Index"), xjust=1, bty="n", seg.len=3, cex=1.2)
		if (png) dev.off()
	}
#browser();return()
	return(list(ntows=table(spp.dat$year), manual=spp.res, tweedie=tweedie.index, model=model.call))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcTweedie


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


