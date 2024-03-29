% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWexPrec_lib.r
\name{find.autocor}
\alias{find.autocor}
\title{find.autocor}
\usage{
find.autocor(autocor.emp, nChainFit, Xt, parMargin, typeMargin)
}
\arguments{
\item{autocor.emp}{target correlation between intensities}

\item{nChainFit}{number of simulations}

\item{Xt}{simulated occurrences, nChainFit x 2 matrix}

\item{parMargin}{parameters of the margins 2 x 3}

\item{typeMargin}{type of marginal distribution: 'EGPD' or 'mixExp'}
}
\value{
\item{scalar}{needed correlation}
}
\description{
finds the autocorrelation leading to observed autocorrelation
}
\author{
Guillaume Evin
}
