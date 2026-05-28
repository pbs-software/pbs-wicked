\name{calcAE}
\alias{calcAE}
\title{
  Calculate Ageing Error
}
\description{
  Calculate ageing error (AE) using the Punt (from Schnute) algorithm in the package \code{'AgeingError'}.
}
\usage{
calcAE(ageDat, atype=c(2,3), species="silvergray rockfish",
   hadley=FALSE, maxage=60, useTMB=TRUE, debug=FALSE,
   png=FALSE, pngres=400, PIN=c(8,8), lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ageDat}{\code{data.frame} -- output from a query to \code{'gfb_age_precision.sql'};\cr
    for example: \code{'aerr417.rda'} for Widow Rockfish (WWR) }
  \item{atype}{\code{numeric} -- age reader type, where
    \tabular{rcl}{
    \code{0} \tab = \tab Unknown\cr
    \code{1} \tab = \tab Final\cr
    \code{2} \tab = \tab Primary\cr
    \code{3} \tab = \tab Precision test\cr
    \code{4} \tab = \tab Secondary
    }}
  \item{species}{\code{character} -- species proper name (lower case);
    follows \code{'gf'} packages convention}
  \item{hadley}{\code{logical} -- if \code{TRUE}, use hadley code}
  \item{maxage}{\code{numeric} -- maximum age used in SS3 (usually the plus group)}
  \item{useTMB}{\code{logical} -- if \code{TRUE}, use a TMB model to fit the age
    data from primary and precision readers via the function \code{run()} rather
    than using the deprecated function \code{RunFn()}}
  \item{debug}{\code{logical} -- if \code{TRUE}, use embedded code
    (not normally executed) to try various routines; code pauses using the
    \code{browser()} function (press `c' to continue when ready)}
  \item{png}{\code{logical} -- if \code{TRUE}, save the figures to a PNG file}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of PNG image}
  \item{PIN}{\code{numeric} -- width and height (inches) of PNG image file}
  \item{lang}{\code{character} -- language of text in PNG image: `e'=english, `f'=french}
}
\details{
  Query to get AE data suitable for this function:
\preformatted{
strSpp = "417" ## Widow Rockfish
getData(fqtName="gfb_age_precision.sql", dbName="GFBioSQL", strSpp=strSpp)
mess = paste0("aerr", strSpp, " ageDat = PBSdat; ", "save(\"aerr", strSpp, 
   "\", file=\"./data/aerr", strSpp, ".rda\")")
eval(parse(text=mess))
}

  Running the function (using Widow as an example):
\preformatted{
species = "widow rockfish"
maxage  = 45
zNA     = is.na(ageDat$age)
agedat  = ageDat[!zNA,]
calcAE(ageDat=agedat, atype=c(2,3), species=species, maxage=maxage)
}
}
\value{
  Look for figures and tables in subdirectories :
  \code{'english'} for figures, \code{'tables'} for ascii tables, 
  \code{'data'} for binary output files.
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K. and Robertson, S. (2008)
\href{https://doi.org/10.1139/F08-111}{%
Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australia's southern and eastern scalefish and shark fishery}.
Can. J. Fish. Aquat. Sci. 65. 1991--2005.

Richards, L.J., Schnute, J.T., Kronlund, A.R. and Beamish, R.J. (1992)
\href{https://doi.org/10.1139/f92-200}{%
Statistical models for the analysis of ageing error}.
Can. J. Fish. Aquat. Sci. 49(9). 1801--1815.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-05-28}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{compAE}}\cr
In package \pkg{PBStools}:\cr
  \code{\link[PBStools]{plotAgeErr}}\cr
In package \pkg{stats}:\cr
  \code{\link[stats]{filter}},
  \code{\link[stats]{lm}},
  \code{\link[stats]{nls}}\cr
In package \pkg{dplyr}:\cr
  \code{\link[dplyr]{select}},
  \code{\link[dplyr]{rename}},
  \code{\link[dplyr]{left_join}},
  \code{\link[dplyr]{mutate}},
  \code{\link[dplyr]{group_by}},
  \code{\link[dplyr]{summarise}}\cr
In package \pkg{ggplot2}:\cr
  \code{\link[ggplot2]{ggplot}}\cr
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_surveys}
% \concept{ ~cpt2 }
% Use only one concept per line.
