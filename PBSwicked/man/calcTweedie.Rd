\name{calcTweedie}
\alias{calcTweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Tweedie CPUE Index Series
}
\description{
  Calculate a Tweedie CPUE abundance index series by fitting
  catch data (including zero-catches) using \code{glmmTMB::glmmTMB}.
}
\usage{
calcTweedie(strSpp="417", gear="bottom trawl", 
   facs=c("location", "depth", "major", "vessel", "latitude"), 
   covs=NULL, rand=NULL, ylim=NULL, refit=TRUE, 
   yrs=1996:2025, major=3:9, vess.yrs=8, vess.trips=7, 
   drange=c(75, 400), dbin=25, effort=6, 
   manualExtract=FALSE, rawIndex=FALSE, ex0=FALSE, 
   year.type="CY", gulf=FALSE, pjs.index=NULL, showCatch=FALSE, 
   png=FALSE, pngres=400, PIN=c(9,6), outnam="tweedie.cpue", lang=c("f","e"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strSpp}{\code{character} -- page number in Hart (1973) that features the selected fish}
  \item{gear}{\code{character} -- description of gear that appears in the
    \code{GF_MERGED_CATCH} table in database \code{GFBioSQL}}
  \item{facs}{\code{character} -- vector of factors to include in the
     generalized linear mixed model (GLMM)}
  \item{covs}{\code{character} -- vector of continuous covariates to
    include in the GLMM}
  \item{rand}{\code{character} -- vector of random effects to include
    in the GLMM}
  \item{ylim}{\code{numeric} -- explicit y-axis limits for CPUE index plot}
  \item{refit}{\code{logical} -- if \code{TRUE}, refit the GLMM;
    otherwise, use fitted model stored in \code{.PBStoolEnv} from previous fit
    (only available during the user's R session)}
  \item{yrs}{\code{numeric} -- vector of years to use in the GLMM}
  \item{major}{\code{numeric} -- vector of major area codes (to subset data)}
  \item{vess.yrs}{\code{numeric} -- minimum number of years that vessels
    must be active in the \code{strSpp} fishery (to subset data)}
  \item{vess.trips}{\code{numeric} -- minimum number of trips per year
    that vessels must catch \code{strSpp} (to subset data)}
  \item{drange}{\code{numeric} -- depth range of \code{strSpp} (to subset data)}
  \item{dbin}{\code{numeric} -- size of depth bin (meters) to bin continuous depths}
  \item{effort}{\code{numeric} -- maximum number of hours allowed per tow (to subset data)}
  \item{manualExtract}{\code{logical} -- if \code{TRUE}, use a manual extract
    method to get the CPUE index series; otherwise, use the \code{emmeans}
    (estimated marginal means) package in R}
  \item{rawIndex}{\code{logical} -- if \code{TRUE}, extract the raw CPUE index series
    from the fitted mode; otherwise, extract the centered (standardized) index series}
  \item{ex0}{\code{logical} -- if \code{TRUE}, exclude zero-catch values, which defeats the 
    purpose of using a Tweedie but sometimes this may be warranted (e.g., Gulf Hake)}
  \item{year.type}{\code{character} -- either \code{'CY'} for calendar year or
    \code{'FY'} for fishing year}
  \item{gulf}{\code{logical} -- if \code{TRUE}, use specific tweaks to fit the sparse
    data available for Pacific Hake in the Stait of Georgia (aka `Gulf')}
  \item{pjs.index}{\code{numeric} -- data frame containing fitted indices from a 
    delta-lognormal (aka hurdle) model; column names : \code{year est lwr upr}}
  \item{showCatch}{\code{logical} -- if \code{TRUE}, plot comparison of reported catch to 
    fitted catch (not tested recently); used in past for Gulf Hake}
  \item{png}{\code{logical} -- if \code{TRUE}, send the figures to \code{PNG} files}
  \item{pngres}{\code{numeric} -- resolution of \code{PNG} figure (pixels/inch)}
  \item{PIN}{\code{numeric} -- vector (width, height) of \code{PNG} figure (inches)}
  \item{outnam}{\code{character} -- base name for \code{PNG} file}
  \item{lang}{\code{character} -- vector of languages (\code{'e'}=english, \code{'f'}=french)
    for text and numeric delimiters in \code{PNG} file}
}
\details{
  Mimics code used by Sean Anderson in the packages \code{'gfdata'} and \code{'gfplots'}.
  The SQL query \code{\Sexpr[stage=render]{system.file( "sql/get-cpue-index.sql", package = "PBSwicked" )}}
  replaces \code{'gfdata::get_cpue_index'} because the \code{gfdata} does not appear to 
  support multiple \code{'major'} codes. Regardless, the supplied SQL script herein
  provides more flexibility and control.
}
\value{
  List of fitted objects :\cr
  \code{ntows} -- number of tows per year\cr
  \code{manual} -- data frame of annual model-derived values using a manual extraction\cr
  \code{tweedie} -- data frame of annual model-derived values using estimated marginal means\cr
  \code{model} -- the model call used to derive the fits
}
\references{
Bolker, B.M., Brooks, M.E., Clark, C.J., Geange, S.W., Poulsen, J.R., Stevens, M.H.H., and White, J.S. (2009).
Generalized linear mixed models: a practical guide for ecology and evolution. 
\emph{Trends in Ecology & Evolution} \bold{24}, 127--135.

Lenth, R. and Piaskowski, J. (2025).
\href{https://CRAN.R-project.org/package=emmeans}{%
emmeans: estimated marginal means, aka least-squares means}.
R package version 2.0.1

McGillycuddy M, Warton DI, Popovic G, and Bolker BM (2025).
\href{https://doi.org/10.18637/jss.v112.i01}{%
Parsimoniously fitting large multivariate random effects in glmmTMB}.
\emph{Journal of Statistical Software}, \bold{112}(1), 1--19.

Brooks, M.E., Kristensen, K., van Benthem, K.J., Magnusson, A., Berg, C.W., Nielsen, A., Skaug, H.J., Maechler, M., and Bolker, B.M. (2017).
\href{https://doi.org/10.32614/RJ-2017-066}{%
glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling}. 
\emph{The R Journal}, \bold{9}(2), 378--400.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish (until 2026)\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-10-26}
}
\note{
  Gemini helped code some of the modelling and graphical sections. 
  He also re-coded the subfunction \code{'calcFY'} (calculate fishing year).
}
\seealso{
In package \pkg{glmmTMB}:\cr
\code{\link[glmmTMB]{glmmTMB}}\cr
In package \pkg{emmeans}:\cr
\code{\link[emmeans]{emmeans}}\cr
In package \pkg{PBSmapping}:\cr
\code{\link[PBSmapping]{clipPolys}}\cr
In package \pkg{PBSmodelling}:\cr
\code{\link[PBSmodelling]{addLegend}}\cr
In package \pkg{PBStools}:\cr
\code{\link[PBStools]{countVec}},
\code{\link[PBStools]{getData}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
\keyword{optimize}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_survey}
% \concept{ ~cpt2 }
% Use only one concept per line.
