\name{wickTMB}
\alias{wickTMB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wick TMB Models
}
\description{
  Run Template Model Builder (TMB) geostatistical models and choose the best one.
}
\usage{
wickTMB(survey_lab="SYN QCS", species="silvergray rockfish", 
   tmbModels, datnam, modnam, idxnam, verbose=TRUE, plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survey_lab}{\code{character} -- label for survey series}
  \item{species}{\code{character} -- genus species for specific rockfish}
  \item{tmbModels}{\code{list} -- list of TMB model functions}
  \item{datnam}{\code{character} -- input name of binary RDA file which contains data from \code{wickSSD}}
  \item{modnam}{\code{character} -- output name of binary RDA file to store all model results}
  \item{idxnam}{\code{character} -- output name of binary RDA file to store chosen index results}
  \item{verbose}{\code{logical} -- if \code{TRUE}, spew summary results onto the R console}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot the results}
}
\details{
  Wicks (fits) 14 TMB models to fit survey data using geospatial information.
}
\value{
  An index data frame with annual records and colums of index-related measures:
  \itemize{
    \item \code{year} -- calendar year (numeric)
    \item \code{est} -- estimated abundance index (numeric)
    \item \code{lwr} -- lower CI of estimate (numeric)
    \item \code{upr} -- upper CI of estimate (numeric)
    \item \code{log_est} -- natural log of estimate (numeric)
    \item \code{se} -- standard error of LN estimate (numeric)
    \item \code{type} -- type of estimate (character), e.g., `index'
  }
  
  Saves binary RDA containing named list of model results for all 14 models.\cr
  RDA example: \code{model.SYN.QCS.silvergray.rda}

  Saves binary RDA containing model predictions (geostatistical indices) for the chosen model.\cr
  RDA example: \code{index.SYN.QCS.silvergray.rda}

  List of objects saved in the index RDA:
  \itemize{
    \item \code{index} -- predicted (fitted) survey series index values (data.frame)
    \item \code{predictions} -- prediction list object with names:
      `data', `report', `obj', `fit_obj', `pred_tmb_data'
    \item \code{nd} -- new data for predictions (data.frame)
    \item \code{survey_grid} -- four columns: `X', `Y', `area', `depth' (data.frame)
    \item \code{surv_wUTM_TMB} -- survey sets with UTM coordinates finalised (data.frame)
    \item \code{m} -- chosen model (class sdmTMB)
    \item \code{mnam} -- name of chosen model
  }
}
\author{
  \href{mailto:nick.fisch@dfo-mpo.gc.ca}{Nick Fisch}, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-12-11}
}
\note{
  Original code from Nick Fisch to determine geo-standardized synoptic survey indices.\cr
  Uses Hadley Wickham pseudo-R packages.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  In package \pkg{PBSwicked}:\cr
  \code{\link[PBSwicked]{wickSSD}},
  \code{\link[PBSwicked]{tmbModels}},
  \code{\link[PBSwicked]{pickIDX}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_geostat.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
