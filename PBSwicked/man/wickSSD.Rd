\name{wickSSD}
\alias{wickSSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wick Survey Series Data (SSD)
}
\description{
  Wick (extract) survey series data (SSD) for a species from the
  GFBioSQL database (housed at DFO).
}
\usage{
wickSSD(survey_lab="SYN QCS", species="silvergray rockfish", datnam, plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survey_lab}{\code{character} -- label for survey series}
  \item{species}{\code{character} -- genus species for specific rockfish}
  \item{datnam}{\code{character} -- optional output name of binary RDA file to store the data}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot the results}
}
\details{
  Wicks (extracts) survey data to derive geostatistical (geospatial)
  synoptic survey index series.
}
\value{
  Returns nothing, but saves data objects to a binary RDA file, for example:\cr
  \code{data.SYN.QCS.silvergray.rda}

  List of objects returned:
  \itemize{
    \item \code{ssids} -- survey series ID numbers in GFBioSQL
    \item \code{rawindex} -- swept-area, design-based index series
    \item \code{surv} -- survey sets from GFBioSQL
    \item \code{surv_wutm} -- survey sets but with missing depths interpolated
    \item \code{surv_wUTM} -- survey sets with UTM coordinates finalised
    \item \code{mesh} -- polygonal mesh for interpolation; 
      looks like a Voronoi (Dirichlet) tesselation
  }
}
\author{
  \href{mailto:nick.fisch@dfo-mpo.gc.ca}{Nick Fisch}, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-13}
}
\note{
  Original code from Nick Fisch to determine geo-standardized synoptic survey indices.\cr
  Uses Hadley Wickham pseudo-R packages and pixie dust.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Function Calls}{
	To fit geospatial model to survey series, call the following functions:\cr
	\enumerate{
		\item \code{wickSSD(survey_lab="SYN QCS", species="widow rockfish")}
			\itemize{
				\item extracts survey series data (SSD) -- \emph{takes forever}
				\item saves SSD to a binary file to avoid running the query over and over
			}
		\item \code{wickTMB(survey_lab="SYN QCS", species="widow rockfish", tmbModels=tmbModels)}
			\itemize{
				\item fits 14 TMB models to the geospatial data -- \emph{takes forever}
				\item chooses the best model using AIC and sanity checks
			}
		\item \code{pickIDX(survey_lab="SYN QCS", species="widow rockfish", plot.maps=T, png=T)}
			\itemize{
				\item plots abundance indices, both geospatial and design-based
				\item writes the indices (standardised by geometric mean) to a CSV
			}
		}
}
\seealso{
  In package \pkg{PBSwicked}:\cr
  \code{\link[PBSwicked]{wickTMB}},
  \code{\link[PBSwicked]{tmbModels}},
  \code{\link[PBSwicked]{pickIDX}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_geostat}
% \concept{ ~cpt2 }
% Use only one concept per line.
