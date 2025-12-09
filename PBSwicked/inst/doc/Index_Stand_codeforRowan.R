
#Code for Rowan to geo-standardize synoptic indices

library(sdmTMB)
library(sp)
library(ggplot2)
library(dplyr)
library(akima)

gfdata::get_ssids()

#Select species and survey
survey_lab<-'SYN QCS'

if(survey_lab=="SYN QCS"){
  ssid<-1} else if (survey_lab=="SYN HS"){
    ssid<-3} else if (survey_lab=="SYN WCVI"){
      ssid<-4} else if (survey_lab=="SYN WCHG"){
        ssid<-16}
#16=West Coast Haida Gwaii Synoptic Bottom Trawl SYN WCHG (8km mesh)
#1=Queen Charlotte Sound Synoptic Bottom Trawl
#3=Hecate Strait Synoptic Bottom Trawl
#4=West Coast Vancouver Island Synoptic Bottom Trawl

#Here is where one would change species
species<-'silvergray rockfish'

#Get data, note that WCHG removed 2014 data because of < half of tows were completed. 
surv <- gfdata::get_survey_sets(species, ssid = ssid)
if (survey_lab=="SYN WCHG"){surv<-surv[-which(surv$year==2014),]}

surv$present<-NA
surv$present[surv$density_kgpm2>0]<-1
surv$present[surv$density_kgpm2==0]<-0
surv$area_swept = surv$tow_length_m * surv$doorspread_m 
surv$area_swept1 = surv$doorspread_m  * (surv$speed_mpm  * surv$duration_min)
surv$area_swept<-ifelse(is.na(surv$area_swept),surv$area_swept1,surv$area_swept)
surv$log_area_swept<-log(surv$area_swept/1e5)

surv<-as.data.frame(cbind(surv$year,surv$catch_weight,surv$log_area_swept,surv$present,surv$depth_m,surv$longitude,surv$latitude))
colnames(surv)<-c('year','catch_weight','log_area_swept','present','depth','X','Y')

#Code to fill in missing depths
surv_wutm <- gfplot:::ll2utm(surv, utm_zone = 9) 
get_depth <- gfplot:::interp_survey_bathymetry(surv_wutm)
#Line that actually fills it in
surv$depth[is.na(surv$depth)] <- get_depth$data$akima_depth[is.na(surv$depth)]

#logging and standardizing depth
surv$log_depth_c<-(log(surv$depth)-mean(log(surv$depth)))/sd(log(surv$depth))
surv$log_depth_c2<-surv$log_depth_c^2

#Getting UTMs back in
surv_wUTM<-gfplot:::ll2utm(surv, utm_zone = 9) 

#Mesh
if (survey_lab=="SYN WCHG"){
  mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 8)  #10km mesh except for WCHG, which had mesh cutoff of 8km, these cutoffs are based on default for gfsynopsis report
} else{ mesh <- make_mesh(surv_wUTM, xy_cols = c("X", "Y"), cutoff = 10)  #10km mesh except for WCHG, which had mesh cutoff of 8km
}
plot(mesh)

#Ok so now fitting models. The way they choose which likelihood to move forward with in gfsynopsis report is one that passes all sanity checks and has lowest AIC
#Delta Lognormal
m_delta_lognormal <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_lognormal(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta Lognormal with depth
m_delta_lognormal_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_lognormal(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Tweedie
m_tweedie <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = tweedie(link = "log"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Tweedie with depth
m_tweedie_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = tweedie(link = "log"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta gamma
m_delta_gamma <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_gamma(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta gamma with depth
m_delta_gamma_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_gamma(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson link gamma
m_delta_poislinkgamma <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson link gamma with depth
m_delta_poislinkgamma_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson link lognormal
m_delta_poislinklognorm <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson link lognormal with depth
m_delta_poislinklognorm_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta generalized gamma distribution
m_delta_gengamma <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_gengamma(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta generalized gamma distribution, with depth
m_delta_gengamma_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_gengamma(),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson-link generalized gamma distribution
m_delta_poislinkgengamma <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year),
  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Delta poisson-link generalized gamma distribution, with depth
m_delta_poislinkgengamma_depth <- sdmTMB(
  data = surv_wUTM,
  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
  offset="log_area_swept",
  spatial='on',
  spatiotemporal='iid',
  anisotropy=TRUE,
  silent=FALSE,
  control = sdmTMBcontrol(newton_loops = 1))

#Use AIC to compare and pick best model
AIC_table<-AIC(m_delta_lognormal,m_delta_lognormal_depth,m_tweedie,m_tweedie_depth,m_delta_gamma,m_delta_gamma_depth,
               m_delta_poislinkgamma,m_delta_poislinkgamma_depth,m_delta_poislinklognorm,m_delta_poislinklognorm_depth,
               m_delta_gengamma,m_delta_gengamma_depth,m_delta_poislinkgengamma,m_delta_poislinkgengamma_depth)
AIC_table[order(AIC_table$AIC),]

#Top AIC model that passes sanity checks
m<-m_delta_gengamma_depth

#However if you have time, should do cross validation instead

#Make sure it passes sanity checks, if not look at next best model
sanity(m)

m$sd_report
summary(m)
tidy(m,conf.int=TRUE)
tidy(m,effects='fixed')
tidy(m,effects='ran_pars',conf.int=TRUE)

#Look at residuals
set.seed(5)
#If its a two-step model, 1=presence-absence, 2=positive catch rates
surv_wUTM$resids <- residuals(m, type="mle-mvn", model=2)
#surv_wUTM$resids <- residuals(m, type="mle-mvn")
qqnorm(surv_wUTM$resids)
abline(a=0,b=1)

#Create extrapolation grid for prediction to sum over to make eventual index 
survey_grid <- gfplot::synoptic_grid |>
  filter(survey == survey_lab) |>
  select(X, Y, area = cell_area, depth=depth)

nd <- sdmTMB::replicate_df(dat = survey_grid, time_name = "year", time_values = unique(surv_wUTM$year))
#Making the standardized variable for the prediction. Remember this is using the mean log depth of the original dataset
nd$log_depth_c <- (log(nd$depth)-mean(log(surv_wUTM$depth)))/sd(log(surv_wUTM$depth))
nd$log_depth_c2 <- nd$log_depth_c^2

#Predicting on new grid
predictions <- predict(m, newdata = nd, return_tmb_object = TRUE)
#Getting the index
index <- get_index(predictions, area = nd$area, bias_correct = TRUE)

#Here is the data you would want for your assessment
index

#plotting
theme_set(ggsidekick::theme_sleek())

plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(ind_pred$data, exp(est2)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(ind_pred$data, exp(est_non_rf2)) +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

plot_map(ind_pred$data, omega_s2) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

plot_map(ind_pred$data, epsilon_st2) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

#Raw Index
rawindex <- gfdata::get_survey_index(species, ssid = ssid)
raw<-as.data.frame(cbind(rawindex$year,rawindex$biomass,rawindex$lowerci,rawindex$upperci))
colnames(raw)<-c('year','est','lwr','upr')
raw$type<-'Design'
index$type<-'Model'
index<-full_join(raw,index)

#plotting by Geomean
plot(index$year[index$type=="Design"]-0.1, index$est[index$type=="Design"]/exp(mean(log(index$est[index$type=="Design"]))), type="p", ylim=c(0,6), las=1, ylab="Index", xlab="Year", xlim=c(2002,2025), col=4, pch=16)
points(index$year[index$type=="Model"], index$est[index$type=="Model"]/exp(mean(log(index$est[index$type=="Model"]))),type="b", pch=16)
arrows(index$year[index$type=="Design"]-0.1, index$lwr[index$type=="Design"]/exp(mean(log(index$est[index$type=="Design"]))), index$year[index$type=="Design"]-0.1, index$upr[index$type=="Design"]/exp(mean(log(index$est[index$type=="Design"]))), angle=0, lwd=2, col=4)
lines(index$year[index$type=="Model"], index$lwr[index$type=="Model"]/exp(mean(log(index$est[index$type=="Model"]))), lty=2)
lines(index$year[index$type=="Model"], index$upr[index$type=="Model"]/exp(mean(log(index$est[index$type=="Model"]))), lty=2)
legend("top", c("Geostatistical","Design"), col=c(1,4), lwd=2)




