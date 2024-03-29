% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWexPrec_lib.r
\name{autocor.emp.int}
\alias{autocor.emp.int}
\title{autocor.emp.int}
\usage{
autocor.emp.int(rho, nChainFit, Xt, parMargin, typeMargin)
}
\arguments{
\item{rho}{autocorrelation of the AR(1) process}

\item{nChainFit}{number of simulated variates}

\item{Xt}{simulated occurrences, nChainFit x 2 matrix}

\item{parMargin}{parameters of the margins 2 x 3}

\item{typeMargin}{type of marginal distribution: 'EGPD' or 'mixExp'}
}
\value{
\item{scalar}{correlation between simulated intensities}
}
\description{
Finds empirical autocorrelations (lag-1) between intensities corresponding
to a degree of autocorrelation of an AR(1) process
}
\author{
Guillaume Evin
}
