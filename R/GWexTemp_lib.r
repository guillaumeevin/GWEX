###===============================###===============================###
### Guillaume Evin
### 19/05/2017, Grenoble
###  IRSTEA
### guillaume.evin@irstea.fr
###
###  Provide utilities for the estimation and simulation of the GWex
### multi-site temperature model.
###===============================###===============================###


#==============================================================================
# mySgedFit
#
# fit nu and xi parameters of the SGED distribution
#
# @param x vector of numerical values
# @export
# @return \item{estimated parameters}{Vector of two parameters: \eqn{nu} and \eqn{xi}}
#
# @author Guillaume Evin
mySgedFit = function(x){
  start = c(nu = 2, xi = 1)
  loglik = function(param, x) {
    f = -sum(log(fGarch::dsged(x, 0, 1, param[1], param[2])))
  }
  fit = nlminb(start = start, objective = loglik,
               lower = c(0, 0), upper = c(Inf, Inf), x = x)
  return(fit$par)
}


#==============================================================================
# get.nonleap.dm
#
# Day/Month for a non-leap year
#
# @return \item{Vector of 365 dates for a non-leap years}{String with the format day/month}
#
# @author Guillaume Evin
get.nonleap.dm = function(){
  # vec Dates for a non-leap year
  vec.Date.seas = seq(from=as.Date("2001/01/01"),to=as.Date("2001/12/31"),by=1)
  vec.Date.seas.dm = format(vec.Date.seas,'%d%m')

  # return results
  return(vec.Date.seas.dm)
}


#==============================================================================
# get.seasonal
#
# Estimate the annual cycle of temperature. For each day, we apply the function \code{fun} to all temperature data for this day (ex. the mean
# of all January, 1st). This cycle is smoothed with the \code{\link{loess}} function.
#
# @param x temperature data
# @param vec.Dates vector of Dates associated to x
# @param myfun function apply for each day
#
# @export
#
# @return A vector of length 365, the annual cycle of temperature
#
# @author Guillaume Evin
get.seasonal = function(x,vec.Dates,myfun='mean'){
  # choose span values, smoother cycle are preferred for the standard deviation
  if(myfun=="mean"){
    myspan = 0.2
  }else if(myfun=='sd'){
    myspan = 0.5
  }
  
  # day/month for the observed dates
  vec.Dates.dm = format(vec.Dates,'%d%m')

  # day/month for a non-leap year
  nonleap.dm = get.nonleap.dm()

  # for 'normal' days
  seas = vector(length=365)
  for(i.d in 1:365){
    is.d = vec.Dates.dm==nonleap.dm[i.d]
    seas[i.d] = do.call(myfun,list(x=x[is.d],na.rm=T))
  }
  
  # if fun="sd", when there is not many data for one day, the standard
  # deviation can be equal to zero (or when bootstrap methods are applied).
  # In that case, we remove zero values and test if we still have enough
  # nonzero standard deviation
  if(myfun=='sd'){
    seas[seas==0] = NA
    if(sum(is.na(seas))>200){
      stop("not enough data to fit a seasonal cycle for this state")
    }
  }

  # smooth trend
  df = data.frame(seas=seas,d=1:365)
  dfz = stats::na.omit(df)
  loess.out = loess(seas~d,data=dfz,span = myspan)
  seas.smooth = predict(loess.out, df)
  return(seas.smooth)
}


#==============================================================================
# predictTempCycle
#
# Return the seasonal cycle for an arbitrary period. The seasonal trend \code{season.fit} is repeated to match the sequence given in \code{vec.Dates}. For leap years, the trend returns equivalent values for the 28/02 and 29/02 days.
#
# @param season.fit seasonal trend: vector 365
# @param vec.Dates vector of Dates to predict
#
# @export
#
# @return A vector giving the seasonal cycles for the time period in vec.Dates
#
# @author Guillaume Evin
predictTempCycle = function(season.fit,vec.Dates){
  # length of the simulations
  n = length(vec.Dates)

  # day/month for the dates to simulate
  vec.Dates.dm = format(vec.Dates,'%d%m')

  # day/month for a non-leap year
  nonleap.dm = get.nonleap.dm()

  # for 'normal' days
  seas.sim = vector(length=n)
  for(i.d in 1:365){
    is.d = vec.Dates.dm==nonleap.dm[i.d]
    seas.sim[is.d] = season.fit[i.d]
  }

  # for leap years
  is.28feb = which(nonleap.dm=="2802")
  is.leap = vec.Dates.dm=="2902"
  seas.sim[is.leap] = season.fit[is.28feb]

  return(seas.sim)
}


#==============================================================================
# get.period.fitting.temp
#
# Get 3-month moving window for one month
#
# @param m.char 3-letter name of a month (e.g. 'JAN')
#
# @return return 3-month period corresponding to this month (ex c(1,2,3) for February)
#
# @author Guillaume Evin
get.period.fitting.temp =  function(m.char){
  # m.char:

  # list of months
  list.m = get.list.month()

  # index of this month
  i.m = which(list.m==m.char)

  # return 3-month period corresponding to this month
  vec.m.ext = c(11,12,1:12,1,2)
  return(vec.m.ext[i.m:(i.m+4)])
}


#==============================================================================
# predicttrend
#
# Return the trend for an arbitrary period
#
# @param vec.slope slopes for the 12 months
# @param vec.Dates vector of Dates to predict
#
# @export
#
# @return return a vector giving the long-term trend for the time period in vec.Dates
#
# @author Guillaume Evin
predicttrend = function(vec.slope,vec.Dates){
  # length of the simulations
  n = length(vec.Dates)

  # month
  vec.m = as.numeric(format(vec.Dates,'%m'))

  # years
  vec.y.raw = as.numeric(format(vec.Dates,'%Y'))
  vec.y = vec.y.raw - min(vec.y.raw)
  vec.y.unique = unique(vec.y)

  # compute trends
  trend.sim = vector(length=n)
  for(y in vec.y.unique){
    for(m in 1:12){
      is.per = vec.y==y & vec.m==m
      trend.sim[is.per] = y*vec.slope[m]
    }
  }

  return(trend.sim)
}


#==============================================================================
# fit.GWex.temp
#
# estimate all the parameters for the G-Wex model of temperature
# @references   Evin, G., A.-C. Favre, and B. Hingray. 2018. “Stochastic Generators of Multi-Site Daily Temperature: Comparison of Performances in Various Applications.” 
# Theoretical and Applied Climatology, 1–14. doi.org/10.1007/s00704-018-2404-x.
#
# @param objGwexObs object of class \code{\linkS4class{GwexObs}}
# @param listOption  list with the following fields:
# \itemize{
#   \item \strong{hasTrend}: logical value, do we fit a linear trend for the long-term change, =FALSE by default
#   \item \strong{objGwexPrec}: object of class \code{\linkS4class{GwexObs}} containing precipitation observations. If provided, we assume that temperature must be modelled and simulated according to the precipitation states 'dry' and 'wet'. For each state, a seasonal cycle is fitted (mean and sd).
#   \item \strong{th}: if objGwexPrec is in listOption, th is the threshold above which we consider that a day is wet (e.g. 0.2 mm)
#   \item \strong{typeMargin}: 'SGED' (default) or 'Gaussian': type of marginal distribution.
#   \item \strong{depStation}: logical value, do we apply a Autoregressive Multivariate Autoregressive model (order 1) =TRUE by default
# }
#
# @export
#
# @return a list containing the list of options \code{listOption} and the list of estimated parameters \code{listPar}:
#
# \itemize{
#   \item \strong{list.trend}: Annual trend for the average temperature, for each month: Raw (in \code{lm}) and smoothed (in \code{smooth}).
#   \item \strong{list.par.margin}: List with one element per station. Each element contains the seasonal cycle for the average temperature \code{season.mean}, the corresponding cycle for the standard deviation \code{season.std} and the parameters of the marginal distributions for each month \code{SkewNormal.par}.
#   \item \strong{list.par.dep}: Lst with one element per month. Each element contains the matrix of correlations \eqn{M_0} for the spatial dependence.
# }
# @author Guillaume Evin
fit.GWex.temp = function(objGwexObs,listOption=NULL){
  ######### Check inputs and assign default values ##########
  if(is.null(listOption)){
    listOption = list()
  }else{
    if(!is.list(listOption)) stop('listOption must be a list')
  }
  
  # hasTrend
  if('hasTrend' %in% names(listOption)){
    hasTrend = listOption[['hasTrend']]
    if(!is.logical(hasTrend)) stop('hasTrend must be logical')
  }else{
    hasTrend = T
    listOption[['hasTrend']] = hasTrend
  }
  
  # th: thresold for precipitation (used only if objGwexPrec is present)
  if('th' %in% names(listOption)){
    th = listOption[['th']]
    if(!(is.numeric(th)&(th>=0))) stop('wrong value for th')
  }else{
    th = 0.2
    listOption[['th']] = th
  }
  
  # objGwexPrec: if objGwexPrec is present, we condition the temperature model to precipitation
  # observations/simulations (see Wilks, 2009)
  if('objGwexPrec' %in% names(listOption)){
    condPrec = T
    objGwexPrec = listOption[['objGwexPrec']]
  }else{
    condPrec = F
  }
  listOption[['condPrec']] = condPrec
  
  # isParallel: not for temperature
  listOption[['isParallel']] = FALSE
  
  # typeMargin: 'SGED' (default) or 'Gaussian'
  if('typeMargin' %in% names(listOption)){
    typeMargin = listOption[['typeMargin']]
    if(!typeMargin%in%c('SGED','Gaussian')) stop('typeMargin must be equal to SGED or Gaussian')
  }else{
    typeMargin = 'SGED'
    listOption[['typeMargin']] = typeMargin
  }
  
  # depStation: 
  # - 'MAR1': applies a Multivariate Autoregressive process to include temporal and spatial dependences.
  # - 'Gaussian': Just include a spatial dependence.
  if('depStation' %in% names(listOption)){
    depStation = listOption[['depStation']]
    if(!depStation%in%c('MAR1','Gaussian')) stop('depStation must be equal to MAR1 or Gaussian')
  }else{
    depStation = 'MAR1'
    listOption[['depStation']] = depStation
  }
  
  ######### Initialize some objects ##########
  # Precipitation matrix
  mat.T = objGwexObs@obs
  
  # number of stations
  p = ncol(mat.T)
  
  # Dates
  vec.dates = objGwexObs@date
  n.day = length(vec.dates)
  
  # Years
  vec.y = strftime(vec.dates, "%Y")
  n.y = length(unique(vec.y))
  
  # Months
  vec.month = as.numeric(strftime(vec.dates, "%m"))
  
  # liste des mois
  vec.month.char = get.list.month()
  n.m = length(vec.month.char)
  
  
  # initialise some objects
  list.par.margin = list()
  u = matrix(nrow = n.day, ncol=p)
  
  
  ######## NON-STATIONARITY TREND ###########
  if(hasTrend){
    # estimate the trend by season
    lm.slope = vector(length=n.m)
    for(m in 1:n.m){
      # period for this month
      per.m = get.period.fitting.month(vec.month.char[m])
      is.per = vec.month%in%per.m
      # regional mean
      mat.T.per = apply(mat.T[is.per,,drop=F],1,mean,na.rm=T)
      # annual averages
      t.mean = aggregate(mat.T.per,by=list(y=vec.y[is.per]),FUN=mean)$x
      # apply a linear regression
      lm.model = lm(y~x,data = list(x=1:n.y,y=t.mean))
      # retrieve slope
      lm.slope[m] = lm.model$coefficients[2]
    }
    
    # smooth these slopes
    smooth.slope = lowess(lm.slope)$y
    
    # predicted trend
    t.trend = predicttrend(smooth.slope,vec.dates)
    
    # return list
    list.trend = list(lm = lm.slope, smooth = smooth.slope)
  }else{
    t.trend = NULL
    list.trend = list()
  }
  
  
  #====================== MARGINS ========================
  # for the progress bar
  pb <- txtProgressBar()
  
  Tdetrend = matrix(nrow=n.day,ncol=p)
  for(i.st in 1:p){
    # vector of temperature for this station
    t.st = mat.T[,i.st]
    
    
    ######## NON-STATIONARITY TREND ###########
    if(hasTrend){
      # remove the trend
      t.detrend = t.st - t.trend
    }else{
      # otherwise we do not apply a trend
      t.detrend = t.st
    }
    
    Tdetrend[,i.st] = t.detrend
    
    ######## SEASONALITY ###########
    # empirical estimate of the seasonal cycle for the mean and sd of the temperature
    # we first remove this seasonality for each station
    
    # if we condition on precipitation values, we fit a seasonal cycles for two precipitatio states:
    # (wet / dry)
    if(condPrec){
      t.mean.season = t.std.season = list()
      t.std = vector(length=n.day)
      for(state in c('dry','wet')){
        # filter temperature obs corresponding to this precipitation states
        t.sel = t.detrend
        if(state == 'dry'){
          is.state = objGwexPrec@obs[,i.st]<=th
        }else if(state == 'wet'){
          is.state = objGwexPrec@obs[,i.st]>th
        }
        is.state[is.na(is.state)] = F
        t.sel[!is.state] = NA
        
        # seasonal cycle of the daily mean: average for each day of the year
        t.mean.season[[state]] = get.seasonal(t.sel,vec.dates,"mean")
        t.mean.pred = predictTempCycle(t.mean.season[[state]],vec.dates)
        t.sh = t.sel - t.mean.pred # remove mean
        
        # seasonal cycle of the daily sd: smooth estimate of sd and average for each day of the year
        t.std.season[[state]] = get.seasonal(t.sh,vec.dates,"sd")
        t.std.pred = predictTempCycle(t.std.season[[state]],vec.dates)
        t.std[is.state] = t.sh[is.state]/t.std.pred[is.state] # standardise
      }
    }else{
      # seasonal cycle of the daily mean: average for each day of the year
      t.mean.season = get.seasonal(t.detrend,vec.dates,"mean")
      t.mean.pred = predictTempCycle(t.mean.season,vec.dates)
      t.sh = t.detrend - t.mean.pred # remove mean
      
      # seasonal cycle of the daily sd: smooth estimate of sd and average for each day of the year
      t.std.season = get.seasonal(t.sh,vec.dates,"sd")
      t.std.pred = predictTempCycle(t.std.season,vec.dates)
      t.std = t.sh/t.std.pred # standardise
    }
    
    
    ######## MARGINAL DISTRIBUTION ###########
    # A Skew normal distribution is fitted to the standardized temperature
    
    # initialize list
    list.SkewNormal.par = list()
    
    if(typeMargin=='SGED'){
      # for each month
      for(m in vec.month.char){
        # three month period
        per.m = get.period.fitting.month(m)
        
        # donnees filtrees
        t.std.per = t.std[vec.month%in%per.m]
        t.filt = t.std.per[!is.na(t.std.per)]
        
        # fit SkewNormal distribution
        list.SkewNormal.par[[m]] = mySgedFit(t.filt)
      }
    }
    
    
    # PIT transform (u are the inputs for copula functions)
    if(typeMargin=='Gaussian'){
      u[,i.st] = pnorm(t.std,mean=0, sd=1)
    }else if(typeMargin=='SGED'){
      for(m in 1:12){
        is.m = (vec.month==m)
        par.m = list.SkewNormal.par[[vec.month.char[m]]]
        u[is.m,i.st] = fGarch::psged(t.std[is.m],mean=0, sd=1, nu=par.m[1], xi=par.m[2])
      }
    }
    
    
    ######## ALL PARAMETERS FOR THE MARGINS ###########
    list.par.margin[[i.st]] = list(season.mean = t.mean.season,
                                   season.std = t.std.season,
                                   SkewNormal.par = list.SkewNormal.par)
    
    # progress bar
    setTxtProgressBar(pb, i.st/(p+12))
  }
  
  # detrend temperature data
  list.trend[['Tdetrend']] = Tdetrend
  
  
  #========== TEMPORAL AND SPATIAL DEPENDENCE ============
  list.par.dep = list()
  
  # Gaussian quantiles
  q.gau = qnorm(u)
  
  # for each month
  for(m in vec.month.char){
    # three month period
    per.m = get.period.fitting.month(m)
    
    # donnees filtrees
    q.gau.per = q.gau[vec.month%in%per.m,,drop=F]
    n.day = nrow(q.gau.per)
    
    if(depStation=='MAR1'){
      # inter-site correlations between all pairs of
      # stations, at lag-0 and lag-1
      q.lag = cbind(q.gau.per[2:n.day,],q.gau.per[1:(n.day-1),])
      
      # pairwise estimation of a correlation matrix with the Kendall tau
      corALL = cor(q.lag, method="pearson", use="pairwise.complete.obs")
      
      # inter-site Pearson correlations + its inverse
      M0 = corALL[1:p,1:p]
      M0inv = MASS::ginv(M0)
      
      # lag-1 correlations between pairs of stations
      M1 = corALL[1:p,(p+1):(2*p)]
      
      # covariance matrices of the MAR(1) process (Matalas, 1967)
      A = M1%*%M0inv
      covZ = M0 - M1%*%M0inv%*%t(M1)
      
      
      # ALL PARAMETERS FOR THE MAR(1)
      list.par.dep[[m]] = list(M0=M0,M1=M1,A=A,covZ=covZ)
    }else if(depStation=='Gaussian'){
      # pairwise estimation of a correlation matrix with the Kendall tau
      M0 = cor(q.gau.per, method="pearson", use="pairwise.complete.obs")
      list.par.dep[[m]] = list(M0=M0)
    }
    
    # progress bar
    i.st = i.st+1
    setTxtProgressBar(pb, i.st/(p+12))
  }
  
  # close progress bar
  close(pb)
  
  
  # return options and estimated parameters
  listPar=list(Xt = q.gau,
               list.trend = list.trend,
               list.par.margin=list.par.margin,
               list.par.dep=list.par.dep,
               p=p)
  return(list(listOption=listOption,listPar=listPar))
}


#==============================================================================
# sim.GWex.temp.1it
#
# Simulate one scenario of temperatures from the GWex model
#
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of dates
# @param myseed seed of the random generation, to be fixed if the results need to be replicated
# @param matSimPrec optional: matrix of precipitation simulation if temperature are generated conditionally to precipitation states "dry" and "wet"
#
# @export
# @return \item{matrix}{Temperature simulated for the dates contained in vec.Dates at the different stations}
#
# @author Guillaume Evin
sim.GWex.temp.1it = function(objGwexFit,vec.Dates,myseed,matSimPrec){
  # set seed of random generation
  set.seed(myseed)

  # number of stations
  p = getNbStations(objGwexFit)

  # retrieve option (model conditional to precipitation?)
  condPrec = objGwexFit@fit$listOption$condPrec
  
  # NOTE 04/10/2024: if temperatures are simulated conditionnally to the precipitation, we 
  # define a threshold to identify wet and dry days. In the simulations, we
  # assume a threshold of 0 mm to preserve the characteristics of the GWEX
  # precipitation model (e.g. probability of having a dry day)
  if(condPrec){
    th = 0
  }

  # caracteristics of the time series generated
  n = length(vec.Dates)
  vec.month = as.numeric(strftime(vec.Dates, "%m"))

  # liste des mois
  vec.month.char = get.list.month()
  n.m = length(vec.month.char)

  # initialise quantities
  Yt.Gau = Yt.Pr = Yt.std = Yt.detrend = Yt = matrix(nrow=n,ncol=p)



  ###____ Spatial and temporal dependence between the stations _____###
  # type of dependence
  depStation =  objGwexFit@fit$listOption$depStation

  # If we have an autoregressive process, we simulate one variate for each time step
  if(depStation=='MAR1'){
    # difference of periods
    vec.per = get.list.month()[vec.month]
    change.per = c(TRUE,vec.per[2:n]!=vec.per[1:(n-1)])
    
    # Parameters for the first iteration
    PAR.DEP = objGwexFit@fit$listPar$list.par.dep[[vec.per[1]]]
    
    # iteration 1 for the spatial dependence
    Yt.Gau[1,] = MASS::mvrnorm(n=1, mu=rep(0,p), Sigma=PAR.DEP[['M0']])
    
    for(t in 2:n){
      if(change.per[t]){
        PAR.DEP = objGwexFit@fit$listPar$list.par.dep[[vec.per[t]]]
      }
      
      # generate from the corresponding multivariate distribution
      # t-1
      Yt.Gau.prev = t(PAR.DEP$A%*%Yt.Gau[t-1,])
      
      # generate from a multivariate Gaussian
      inno = MASS::mvrnorm(n=1, mu=rep(0,p), Sigma=PAR.DEP[['covZ']])
      
      # MAR(1)
      Yt.Gau[t,] = Yt.Gau.prev + inno
    }
  }else{
    # Spatial dependence only, we simulate month by month
    for(m in 1:n.m){
      # retrieve parameters
      PAR.DEP = objGwexFit@fit$listPar$list.par.dep[[m]]
      # for all days concerned by this month
      is.m = (vec.month==m)
      n.sim = sum(is.m)
      # simulate the spatial dependence
      Yt.Gau[is.m,] = MASS::mvrnorm(n=n.sim, mu=rep(0,p), Sigma=PAR.DEP[['M0']])
    }
  }
  
  # transformation in probability for Gaussian variates
  Yt.Pr = pnorm(Yt.Gau)
  
  ###____________ Inverse-CDF ____________###
  # parameters
  typeMargin = objGwexFit@fit$listOption$typeMargin
  PAR.margin = objGwexFit@fit$listPar$list.par.margin

  for(i.st in 1:p){
    if(typeMargin=='Gaussian'){
      Yt.std[,i.st] = qnorm(p=Yt.Pr[,i.st], mean=0, sd=1)
    }else if(typeMargin=='SGED'){
      for(m in 1:n.m){
        par.m = PAR.margin[[i.st]]$SkewNormal.par[[m]]
        is.m = (vec.month==m)
        Yt.std[is.m,i.st] = fGarch::qsged(p=Yt.Pr[is.m,i.st], mean=0, sd=1, nu=par.m[1], xi=par.m[2])
      }
    }
  }

  ###___ Add cycle and un-standardize _____###
  for(i.st in 1:p){
    # retrieve seasonal cycles
    season.mean = PAR.margin[[i.st]]$season.mean
    season.std = PAR.margin[[i.st]]$season.std

    # predict seasonal cycles (conditionally to precipitation sim or not)
    if(condPrec){
      x.mean.pred = x.std.pred = vector(length=n)
      for(state in c("dry","wet")){
        if(state == 'dry'){
          is.state = matSimPrec[,i.st,]<=th
        }else if(state == 'wet'){
          is.state = matSimPrec[,i.st,]>th
        }
        x.mean.pred[is.state] = predictTempCycle(season.mean[[state]],vec.Dates[is.state])
        x.std.pred[is.state] = predictTempCycle(season.std[[state]],vec.Dates[is.state])
      }
    }else{
      x.mean.pred = predictTempCycle(season.mean,vec.Dates)
      x.std.pred = predictTempCycle(season.std,vec.Dates)
    }

    # unstandardize
    Yt.detrend[,i.st] = Yt.std[,i.st]*x.std.pred + x.mean.pred
  }

  ###___ Add long-term trends _____###
  hasTrend = objGwexFit@fit$listOption$hasTrend
  if(hasTrend){
    # trend for the simulaed period
    t.trend = predicttrend(objGwexFit@fit$listPar$list.trend$smooth,vec.Dates)
    for(i.st in 1:p){
      # add trend
      Yt[,i.st] = Yt.detrend[,i.st] + t.trend
    }
  }else{
    # otherwise, we simply return Yt.detrend
    Yt = Yt.detrend
  }

  # return results
  return(list(Yt=Yt,Xt=Yt.Gau,Zt=Yt.std,Tdetrend=Yt.detrend))
}