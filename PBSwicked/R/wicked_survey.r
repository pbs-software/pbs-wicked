##==============================================================================
## Module : wicked_survey
## -----------------------
## getHBLL..........Get survey series indices for HBLL surveys

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
