% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWex_lib.r
\name{fitGwexModel}
\alias{fitGwexModel}
\title{fitGwexModel: fit a GWex model to observations.}
\usage{
fitGwexModel(objGwexObs, parMargin = NULL, listOption = NULL)
}
\arguments{
\item{objGwexObs}{an object of class \code{\linkS4class{GwexObs}}}

\item{parMargin}{(not required for temperature) list parMargin where
and each element corresponds to a month (1...12) and contains a 
matrix nStation x 3 of pre-estimated parameters of the marginal distributions
(EGPD or Mixture of Exponentials)}

\item{listOption}{for precipitation, a list with the following fields:
\itemize{
  \item \strong{th}: threshold value in mm above which precipitation observations are considered to be non-zero (=0.2 by default)
  \item \strong{nLag}: order of the Markov chain for the transitions between dry and wet states (=2 by default)
  \item \strong{typeMargin}: 'EGPD' (Extended GPD) or 'mixExp' (Mixture of Exponentials). 'mixExp' by default
  \item \strong{copulaInt}: 'Gaussian' or 'Student': type of dependence for amounts (='Gaussian' by default)
  \item \strong{isMAR}: logical value, do we apply a Autoregressive Multivariate Autoregressive model (order 1) = FALSE by default
  \item \strong{is3Damount}: logical value, do we apply the model on 3D-amount. =FALSE by default
  \item \strong{nChainFit}: integer, length of the runs which are generated during the fitting procedure. =100000 by default
  \item \strong{nCluster}: integer, number of clusters which can be used for the parallel computation
}

and for temperature, a list with the following fields:
\itemize{
  \item \strong{hasTrend}: logical value, do we fit a linear trend for the long-term change, =FALSE by default
  \item \strong{objGwexPrec}: object of class \code{\linkS4class{GwexObs}} containing precipitation observations. If provided, we assume 
  that temperature must be modelled and simulated according to the precipitation states 'dry' and 'wet'. For each state, a seasonal cycle is fitted (mean and sd).
  \item \strong{typeMargin}: 'SGED' (default) or 'Gaussian': type of marginal distribution.
  \item \strong{depStation}: 'MAR1' (default) or 'Gaussian': MAR1 (Multivariate Autoregressive model order 1) for the spatial and temporal dependence or 'Gaussian' 
  for the spatial dependence only.
}}
}
\value{
Return an object of class \code{\linkS4class{GwexFit}} with:
\itemize{
  \item \strong{p}: The number of station,
  \item \strong{version}: package version,
  \item \strong{variable}: the type of variable,
  \item \strong{fit}: a list containing the list of options \code{listOption}
   and the list of estimated parameters \code{listPar}.
}
}
\description{
fitGwexModel: fit a GWex model to observations.
}
\examples{
# Format dates corresponding to daily observations of precipitation and temperature
vecDates = seq(from=as.Date("01/01/2005",format="\%d/\%m/\%Y"),
to=as.Date("31/12/2014",format="\%d/\%m/\%Y"),by='day')

###############################################################
#               FIT THE PRECIPITATION MODEL
###############################################################

# Format observations: create a Gwex object for one station only to show a quick
# example. The syntax is similar for multi-site applications.
myObsPrec = GwexObs(variable='Prec',date=vecDates,obs=dailyPrecipGWEX[,1,drop=FALSE])

# Fit precipitation model with a threshold of 0.5 mm to distinguish wet and dry 
# states (th) and keep default options otherwise, e.g. a Gaussian
# copula for the spatial dependence (copulaInt) and a mixExp distribution for 
# marginal intensities ('typeMargin')
myParPrec = fitGwexModel(myObsPrec,listOption=list(th=0.5))
myParPrec # print object

###############################################################
#     FIT THE TEMPERATURE MODEL, COND. TO PRECIPITATION
###############################################################
# Format observations: create a G-Wex object
myObsTemp = GwexObs(variable='Temp',date=vecDates,obs=dailyTemperGWEX)

# Fit temperature model with a long-term linear trend ('hasTrend'), Gaussian margins 
# ('typeMargin') and Gaussian spatial dependence ('depStation')
myParTemp = fitGwexModel(myObsTemp,listOption=list(hasTrend=TRUE,typeMargin='Gaussian',
depStation='Gaussian'))
myParTemp # print object
}
\author{
Guillaume Evin
}
