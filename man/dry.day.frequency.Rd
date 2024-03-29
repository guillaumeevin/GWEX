% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWexPrec_lib.r
\name{dry.day.frequency}
\alias{dry.day.frequency}
\title{dry.day.frequency}
\usage{
dry.day.frequency(mat.prec, th)
}
\arguments{
\item{mat.prec}{matrix of precipitation (possibly for one month/period)}

\item{th}{threshold above which we consider that a day is wet (e.g. 0.2 mm)}
}
\value{
\item{vector of numeric}{dry day frequencies}
}
\description{
Estimate the dry day frequency (proportion of dry days) for all stations
}
\author{
Guillaume Evin
}
