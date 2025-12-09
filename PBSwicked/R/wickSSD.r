## wickSSD -----------------------------2025-12-05
## Get survey series data (SSD) for a species
##  using wicked Hadley Wickham code
## Code from NF to geo-standardize synoptic indices
## ------------------------------------------NF|RH
wickSSD <- function(survey_lab="SYN QCS", species="silvergray rockfish", outnam, plot=TRUE)
{
	if (missing(outnam))
		outnam = paste0("data.",gsub("\\s+",".",survey_lab),".",strsplit(species,split="\\s+")[[1]][1])
	if (!file.exists(paste0(outnam,".rda"))) {
		ssids = gfdata::get_ssids()
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

		rawindex <- gfdata::get_survey_index(species, ssid = ssid)  ## do this now rather than later in 'wickIDX()'

		## Get data, note that WCHG removed 2014 data because of < half of tows were completed. 
		surv <- gfdata::get_survey_sets(species, ssid = ssid)  ## is this the query that takes forever?
		if (survey_lab=="SYN WCHG"){
			surv<-surv[-which(surv$year==2014),]
		}
		surv$present<-NA
		surv$present[surv$density_kgpm2>0]<-1
		surv$present[surv$density_kgpm2==0]<-0
		surv$area_swept = surv$tow_length_m * surv$doorspread_m 
		surv$area_swept1 = surv$doorspread_m  * (surv$speed_mpm  * surv$duration_min)
		surv$area_swept<-ifelse(is.na(surv$area_swept),surv$area_swept1,surv$area_swept)
		surv$log_area_swept<-log(surv$area_swept/1e5)

		surv<-as.data.frame(cbind(surv$year,surv$catch_weight,surv$log_area_swept,surv$present,surv$depth_m,surv$longitude,surv$latitude))
		colnames(surv)<-c('year','catch_weight','log_area_swept','present','depth','X','Y')

		## Code to fill in missing depths
		surv_wutm <- gfplot:::ll2utm(surv, utm_zone = 9) 
		get_depth <- gfplot:::interp_survey_bathymetry(surv_wutm)
		#Line that actually fills it in
		surv$depth[is.na(surv$depth)] <- get_depth$data$akima_depth[is.na(surv$depth)]

		#logging and standardizing depth
		surv$log_depth_c<-(log(surv$depth)-mean(log(surv$depth)))/sd(log(surv$depth))
		surv$log_depth_c2<-surv$log_depth_c^2
		
		#Getting UTMs back in
		surv_wUTM<-gfplot:::ll2utm(surv, utm_zone = 9) 

		if (survey_lab=="SYN WCHG"){
			mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 8)  #10km mesh except for WCHG, which had mesh cutoff of 8km, these cutoffs are based on default for gfsynopsis report
		} else { 
			mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 10)  #10km mesh except for WCHG, which had mesh cutoff of 8km
		}
		## save this mess
		save("ssids", "rawindex", "surv", "surv_wutm", "surv_wUTM", "mesh", file=paste0(outnam,".rda"))
	} else {
		load(paste0(outnam, ".rda"))
	}
	## Mesh (RH 251205) cannot seem to change appearance of mesh with par() arguments
	if (plot)
		plot(mesh)
#browser();return()
}
require(sdmTMB)
require(sp)
require(ggplot2)
require(dplyr)
require(akima)

wickSSD(survey_lab="SYN QCS", species="silvergray rockfish")

