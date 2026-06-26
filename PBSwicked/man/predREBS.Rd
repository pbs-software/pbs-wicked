\name{predREBS}
\alias{predREBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict REBS from Location
}
\description{
  Predict the proportion of Rougheye Rockfish (RER) from REBS data
  using two-dimensional space (longitude, latitude) along the BC coast.
}
\usage{
predREBS(dat, indvar=c("X","Y"), cutoff=5, rmSM=TRUE, 
   col=c("navyblue", "blue", "lightseagreen", "gold", "orangered"), 
   ncells=250, addpts=FALSE, outnam="predREBS", 
   png=FALSE, pngres=400, PIN=c(9, 9.5), lang=c("f","e"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{\code{data.frame} -- data from an SQL query with 
    subsequent massage (contact NF); data needs column labelled
    \code{'species'}, which consists of 0=not RER and 1=RER}
  \item{indvar}{\code{character} -- colnames that indicate longitude
    and latitude, (in this case, \code{'X'} indicates UTM Easting and
    \code{'Y'} indicates UTM Northing, both in km)}
  \item{cutoff}{\code{numeric} -- the minimum distance between mesh
    vertices in units of the x and y coordinates}
  \item{rmSM}{\code{logical} -- if \code{TRUE}, remove seamount data}
  \item{col}{\code{character} -- vector of color names used in 
    \code{'colorRampPalette'}}
  \item{ncells}{\code{numeric} -- number of grid cells to cover the
    extent of the data}
  \item{addpts}{\code{logical} -- if \code{TRUE}, add the data point locations}
  \item{outnam}{\code{character} -- name of output png file}
  \item{png}{\code{logical} -- if \code{TRUE}, send the figure to a \code{png} file}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of PNG image}
  \item{PIN}{\code{numeric} -- width and height (inches) of PNG image file}
  \item{lang}{\code{character} -- language of text in PNG image: `e'=english, `f'=french}
}
\details{
  Fits a spatial model to presence-absence data, where species=1
  indicates RER detemined by genetics.
  
}
\value{
  Nada pour le moment.
}
\author{
  \href{mailto:nick.fisch@dfo-mpo.gc.ca}{Nick Fisch}, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-05-07}
}
\note{
  Very basic model for illustrative purposes. 
  Proportion of Blackspotted Rockfish (BSR) is implied to 
  equal 1-RER, but \code{species}=0 may include hybrids
  in addition to BSR.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  In package \pkg{PBSwicked}:\cr
  \code{\link[PBSwicked]{getHBLL}},
  \code{\link[PBSwicked]{pickIDX}}\cr
  In package \pkg{sdmTMB}:\cr
  \code{\link[sdmTMB]{make_mesh}},
  \code{\link[sdmTMB]{sdmTMB}}\cr
  In package \pkg{akima}:\cr
  \code{\link[akima]{interp}}\cr
  In package \pkg{graphics}:\cr
  \code{\link[graphics]{image}}\cr
  In package \pkg{PBStools}:\cr
  \code{\link[PBStools]{addStrip}},
  \code{\link[PBStools]{changeLangOpts}},
  \code{\link[PBStools]{clearFiles}},
  \code{\link[PBStools]{createFdir}},
  \code{\link[PBStools]{zapSeamounts}}\cr
  In package \pkg{PBSmapping}:\cr
  \code{\link[PBSmapping]{addPolys}},
  \code{\link[PBSmapping]{convUL}},
  \code{\link[PBSmapping]{plotMap}}\cr
  In package \pkg{PBSmodelling}:\cr
  \code{\link[PBSmodelling]{addLabel}},
  \code{\link[PBSmodelling]{expandGraph}},
  \code{\link[PBSmodelling]{unpackList}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
\concept{wicked_survey}
% Use only one concept per line.
