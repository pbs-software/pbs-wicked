\name{pickIDX}
\alias{pickIDX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create Geospatial Index Figure
}
\description{
  Pick (create pictures of) prediction maps and geospatial index series.
}
\usage{
pickIDX(survey_lab="SYN QCS", species="silvergray rockfish", 
   datnam, idxnam, pixnam, plot.maps=FALSE, 
   png=FALSE, pngres=400, PIN=c(8, 6), lang=c("e"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survey_lab}{\code{character} -- label for survey series}
  \item{species}{\code{character} -- genus species for specific rockfish}
  \item{datnam}{\code{character} -- input name of binary RDA file 
    which contains data from \code{wickSSD}}
  \item{idxnam}{\code{character} -- input name of binary RDA file 
    which contains chosen index results from \code{wickTMB}}
  \item{pixnam}{\code{character} -- output name of binary RDA file to store figures}
  \item{plot.maps}{\code{logical} -- if \code{TRUE}, plot the predictiion maps}
  \item{png}{\code{logical} -- if \code{TRUE}, save the figures 
    (prediction maps, geospatial index series) to a PNG file}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of PNG image}
  \item{PIN}{\code{numeric} -- width and height (inches) of PNG image file}
  \item{lang}{\code{character} -- language of text in PNG image: `e'=english, `f'=french}
}
\details{
  Geospatial index series plot compares swept-area, design-based 
  indices (PJS plots) with index series fitted by a geostatistical
  model (annual estimates and CI bounds).
  The geostatistical model selected is the one (out of 14) with the
  lowest Akaike Information Criterion (AIC), a metric that balances
  goodness of fit with model complexity.
  Additionally, the selected model must be deemed `sane' (via sanity
  checks), using supplemental diagnostics to ensure the model's
  outputs, performance, and data inputs are reasonable and reliable.
}
\value{
  Nothing | rien | nada

  Saves RDA containing images of prediction maps
  created by Wickham's galactic algorithm \code{ggplot}.\cr
  RDA example: \code{plot.SYN.QCS.silvergray.rda}
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
\note{
  Original code from Nick Fisch to determine geo-standardized synoptic survey indices.\cr
  Uses Hadley Wickham pseudo-R packages.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  In package \pkg{PBSwicked}:\cr
  \code{\link[PBSwicked]{wickSSD}},
  \code{\link[PBSwicked]{wickTMB}},
  \code{\link[PBSwicked]{tmbModels}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_geostat.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
