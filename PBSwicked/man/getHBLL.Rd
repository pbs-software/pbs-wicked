\name{getHBLL}
\alias{getHBLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get HBLL Survey Data
}
\description{
  Get the Hard-bottom Longline survey data for a selected species.
}
\usage{
getHBLL(strSpp="417", timestamp="260502")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strSpp}{\code{character} -- Hart code (page number in Hart) for
  groundfish or alpha-numeric code for others like invertebrates}
  \item{timestamp}{\code{character} -- time stamp, typicaly using the format
    'YYMMDD', but can be whatever string the user wishes}
}
\details{
  The algorithm used by Dana Haggarty in 2022.
  The code is largely reliant on Wickham's composting methods.
  Dana's routine accounts for hook competition along a longline set.
}
\value{
  The final H&L data set where each record is a longline set for the
  species selected, specifying hooks (empty, baited, total) and the
  adjustment factor to account for hook competition.
}
\author{
  \href{mailto:dana.haggarty@dfo-mpo.gc.ca}{Dana Haggarty}, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-05-06}
}
\note{
  Takes an ungodly amount of time to query the database.
  Function saves the query results so one does not have to 
  repeat the two queries if debugging.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSdata}:\cr
  \code{\link[PBSdata]{species}}\cr
In package \pkg{gfdata}:\cr
  \code{\link[gfdata]{get_ll_hook_data}}, 
  \code{\link[gfdata]{get_survey_sets}}\cr
In package \pkg{dplyr}:\cr
  \code{\link[dplyr]{left_join}}, 
  \code{\link[dplyr]{mutate}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{IO}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{wicked_survey}
% \concept{ ~cpt2 }
% Use only one concept per line.
