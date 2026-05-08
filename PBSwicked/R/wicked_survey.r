##==============================================================================
## Module : wicked_survey
## -----------------------
## getHBLL..........Get survey series indices for HBLL surveys
## predREBS.........Predict probability of RER (BSR is the inverse)

##-----Supplementary hidden functions-----
##===============================================================================
## Note: not sure whether to rename functions so it's obvious that they are wicked
## See object `wickiverb' for possible guidance (RH 260506)

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


