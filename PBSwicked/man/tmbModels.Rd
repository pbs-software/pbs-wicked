\name{tmbModels}
\alias{tmbModels}
\alias{m_delta_gamma}
\alias{m_delta_gamma_depth}
\alias{m_delta_gengamma}
\alias{m_delta_gengamma_depth}
\alias{m_delta_lognormal}
\alias{m_delta_lognormal_depth}
\alias{m_delta_poislinkgamma}
\alias{m_delta_poislinkgamma_depth}
\alias{m_delta_poislinkgengamma}
\alias{m_delta_poislinkgengamma_depth}
\alias{m_delta_poislinklognorm}
\alias{m_delta_poislinklognorm_depth}
\alias{m_tweedie}
\alias{m_tweedie_depth}
\docType{data}
\title{
  TMB Models
}
\description{
  A list of sdmTMB model functions used to fit geospatial data.
}
\usage{tmbModels}
\format{
  The format is:\cr
  List of 14 functions, all with argument \code{surv_wUTM}\cr
  01_Delta_Lognormal\cr
  02_Delta_Lognormal_with_depth\cr
  03_Tweedie\cr
  04_Tweedie_with_depth\cr
  05_Delta_gamma\cr
  06_Delta_gamma_with_depth\cr
  07_Delta_poisson_link_gamma\cr
  08_Delta_poisson_link_gamma_with_depth\cr
  09_Delta_poisson_link_lognormal\cr
  10_Delta_poisson_link_lognormal_with_depth\cr
  11_Delta_generalized_gamma_distribution\cr
  12_Delta_generalized_gamma_distribution_with_depth\cr
  13_Delta_poisson_link_generalized_gamma_distribution\cr
  14_Delta_poisson_link_generalized_gamma_distribution_with_depth\cr
}
\details{
  Fit a spatial or spatiotemporal generalized linear mixed effects model (GLMM)
  with the TMB (Template Model Builder) R package and the SPDE (stochastic 
  partial differential equation) approximation to Gaussian random fields.
}
\source{
  CRAN: \href{https://cran.r-project.org/package=sdmTMB}{sdmTMB: Spatial and Spatiotemporal SPDE-Based GLMMs with `TMB'}
}
\references{
	Anderson, SC, Ward, EJ, English, PA, Barnett, LAK, and Thorson, JT. (2024)
	\href{https://doi.org/10.1101/2022.03.24.485545}{%
	sdmTMB: An R package for fast, flexible, and user-friendly Generalized Linear Mixed Effects Models with spatial and spatiotemporal random fields.}
	\emph{bioRxiv}, elocation-id 2022.03.24.485545
}
\author{
  \href{mailto:nick.fisch@dfo-mpo.gc.ca}{Nick Fisch}, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-12-12}
}
\keyword{datasets}
\concept{wicked_geostat.r}

