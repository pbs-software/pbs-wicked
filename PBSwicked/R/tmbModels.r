## NB : Ok so now fitting 14 models.
## The way they choose which likelihood to move forward with in 
##  gfsynopsis report is one that passes all sanity checks and has lowest AIC
## (RH 251205) Crikey! this will take forever!
tmbModels <- list()
tmbModels[['01_Delta_Lognormal']] <- 
	m_delta_lognormal <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_lognormal(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['02_Delta_Lognormal_with_with_depth']] <- 
	m_delta_lognormal_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_lognormal(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['03_Tweedie']] <- 
	m_tweedie <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = tweedie(link = "log"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['04_Tweedie_with_depth']] <- 
	m_tweedie_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = tweedie(link = "log"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['05_Delta_gamma']] <- 
	m_delta_gamma <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['06_Delta_gamma_with_depth']] <- 
	m_delta_gamma_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['07_Delta_poisson_link_gamma']] <- 
	m_delta_poislinkgamma <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['08_Delta_poisson_link_gamma_with_depth']] <- 
	m_delta_poislinkgamma_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['09_Delta_poisson_link_lognormal']] <- 
	m_delta_poislinklognorm <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['10_Delta_poisson_link_lognormal_with_depth']] <- 
	m_delta_poislinklognorm_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_lognormal(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['11_Delta_generalized_gamma_distribution']] <- 
	m_delta_gengamma <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gengamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['12_Delta_generalized_gamma_distribution_with_depth']] <- 
	m_delta_gengamma_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gengamma(),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['13_Delta_poisson_link_generalized_gamma_distribution']] <- 
	m_delta_poislinkgengamma <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year),
	  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}
tmbModels[['14_Delta_poisson_link_generalized_gamma_distribution_with_depth']] <- 
	m_delta_poislinkgengamma_depth <- function(surv_wUTM) {
	sdmTMB(
	  data = surv_wUTM,
	  formula = catch_weight ~ 0 + as.factor(year) + log_depth_c + log_depth_c2,
	  time = "year", mesh = mesh, family = delta_gengamma(type="poisson-link"),
	  offset="log_area_swept",
	  spatial='on',
	  spatiotemporal='iid',
	  anisotropy=TRUE,
	  silent=FALSE,
	  control = sdmTMBcontrol(newton_loops = 1)
	)}

