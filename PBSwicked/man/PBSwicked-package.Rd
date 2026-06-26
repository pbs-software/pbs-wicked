\name{PBSwicked-package}
\alias{PBSwicked-package}
\alias{PBSwicked}
\docType{package}
\title{
  Wicked Wonders from Another Universe
}
\description{
  By the pricking of my thumbs, something wicked this way comes.
}
\details{
  Provides pseudo-R code (largely based on the Hadley Wickham R-verse that parallels our universe)
  to bedazzle fisheries stock assessment scientists. 
  It brings to mind the Thomas Dolby 1982 song \emph{She Blinded Me with Science}. 
  The functions herein are kept separate from other PBS packages, 
  largely because there are so many R package dependencies.
}
\section{Hadley Wickham}{
  Although the maintainer castigates the use of R-overlay language, 
  largely authored by \href{https://hadley.nz/}{Hadley Wickham}, 
  there is tremendous respect for the man. 
  We share a few commonalities (country of birth, friends of Dorothy, 
  data visualization), which transcend our differences. 
  While Hadley prefers the \href{https://en.wikipedia.org/wiki/Tidyverse}{tidyverse},
  the maintainer prefers an `untidyverse'. 
  The reference to `wicked' herein, a not-so-subtle play on `Wickham', 
  is a tribute to the man's wicked R code in all senses of the word -- 
  \emph{awesome}, \emph{excellent}, \emph{skillful}, \emph{complex}, and
  \emph{slightly bad in a good way}. 
  Most users prefer the tidyverse, and so this package houses morsels from that realm.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-06-26}
}
\keyword{package}
%% Uncomment below to imitate parts of library(help = PBSwicked)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{PBSwicked}}
%\section{Documentation Index}{\packageIndices{PBSwicked}}
\seealso{
In package \pkg{PBSwicked}:\cr
  \code{\link[PBSwicked]{calcAE}} -- calculate ageing error using the NOAA package \code{'AgeingError'}\cr
  \code{\link[PBSwicked]{getHBLL}} -- get the hard-bottom longline survey data for a selected species\cr
  \code{\link[PBSwicked]{pickIDX}} -- pick (create pictures of) prediction maps and geospatial index series\cr
  \code{\link[PBSwicked]{predREBS}} -- predict the proportion of Rougheye Rockfish (RER) from REBS data\cr
  \code{\link[PBSwicked]{tmbModels}} -- list of sdmTMB model functions used to fit geospatial data\cr
  \code{\link[PBSwicked]{wickiverb}} --  a named character vector containing 26 verb choices\cr
  \code{\link[PBSwicked]{wickSSD}} -- wick (extract) survey series data (SSD) for a species from the GFBioSQL database\cr
  \code{\link[PBSwicked]{wickTMB}} --  wick (run) Template Model Builder geostatistical models and choose the best one
}
