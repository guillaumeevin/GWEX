###===============================###===============================###
### Guillaume Evin
### 02/12/2022, Grenoble
###  INRAE-UR ETNA
### guillaume.evin@inrae.fr
###
###  Provide utilities for the estimation and simulation of the GWex
### multi-site precipitation model.
###
### Wilks, D.S. (1998) "Multisite generalization of a daily stochastic
### precipitation generation model", J Hydrol, 210: 178-191
###
### Evin, G., A.-C. Favre, and B. Hingray. 2018. “Stochastic Generation of 
### Multi-Site Daily Precipitation Focusing on Extreme Events.” Hydrol. Earth 
### Syst. Sci. 22 (1): 655–72. https://doi.org/10.5194/hess-22-655-2018.
###===============================###===============================###

#==============================================================================
#' dry.day.frequency
#'
#' Estimate the dry day frequency (proportion of dry days) for all stations
#'
#' @param mat.prec matrix of precipitation (possibly for one month/period)
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#'
#' @return \item{vector of numeric}{dry day frequencies}
#'
#' @author Guillaume Evin
dry.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # dry day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x<=th,na.rm=T),th))
}


#==============================================================================
#' wet.day.frequency
#'
#' Estimate the wet day frequency (proportion of wet days) for all stations
#'
#' @param mat.prec matrix of precipitation (possibly for one month/period)
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#'
#' @return \item{vector of numeric}{wet day frequencies}
#'
#' @author Guillaume Evin
wet.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # wet day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x>th,na.rm=T),th))
}

#==============================================================================
#' lag.trans.proba.vector.byClass
#'
#' Estimate the transition probabilities between wet and dry states, for nlag previous days,
#' for one station
#'
#' @param vec.prec vector nx1 of precipitation for one station
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param nlag number of lag days
#' @param dayScale time resolution
#'
#' @return \item{matrix}{matrix nLag^2 x (nLag+1) of transition probability between dry/wet state.
#'  The first nLag columns indicate the wet/dry states for the previous nLag days}
#'
#' @author Guillaume Evin
lag.trans.proba.vector.byClass = function(vec.prec, # vector nx1 of precipitation for one station
                                          isPeriod, # vector nx1 of logical
                                          th, # threshold
                                          nlag, # number of lag days
                                          dayScale # time resolution
){
  # number lags + 1 
  ndays = nlag+1
  
  # length of input vectors
  nVec = length(vec.prec)
  
  # filter dates where we have the previous nlag days & concerned by this class
  # NOTE: we consider that only the current day should be concerned by the class,
  # not the previous days (too restrictive)
  filter.first.days = 1:nVec>nlag
  isClass = isPeriod
  ind.isClass = which(isClass&filter.first.days)
  
  # is.wet.lag gives the boolean values indicating if the time step 0, 1, ... ,0+nlag
  # are wet or not
  is.wet.lag = matrix(nrow=length(ind.isClass),ncol=ndays)
  for(i.lag in 0:nlag) is.wet.lag[,(i.lag+1)] = vec.prec[ind.isClass+i.lag]>th
  
  # no nan
  nz = apply(is.wet.lag,1,function(x) !any(is.na(x)))
  x.nz = is.wet.lag[nz,]
  
  # matrix of joints probabilities
  comb.pr = expand.grid(lapply(numeric(ndays), function(x) c(F,T)))
  names(comb.pr) = paste0('t',(-nlag:0))
  comb.grid = as.matrix(comb.pr)
  comb.pr$P = NA
  
  for(i.comb in 1:nrow(comb.grid)){
    comb = comb.grid[i.comb,]
    comb.pr$P[i.comb] = mean(apply(x.nz,1, function(x) all(x==comb)))
  }
  
  # matrix of conditional probabilities Pr(i0=wet|i-1=.,i-2=.,...)
  nPr = nrow(comb.pr)
  prF = comb.pr$P[1:(nPr/2)]
  prT = comb.pr$P[(nPr/2+1):nPr]
  PrCond = prT/(prF+prT)
  comb.prT = comb.pr[(nPr/2+1):nPr,1:nlag,drop=F]
  comb.prT$P = PrCond
  
  # return  conditional probabilities
  return(comb.prT)
}


#==============================================================================
#' lag.trans.proba.mat.byClass
#'
#' Estimate the transition probabilites between wet and dry states, for nlag previous days,
#' for all stations
#'
#' @param mat.prec matrix of precipitation
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param nlag number of lag days
#' @param dayScale time resolution
#'
#' @return \item{list}{list with one item per station, where each item is a
#' matrix nLag^2 x (nLag+1) of transition probability between dry/wet state.
#'  The first nLag columns indicate the wet/dry states for the previous nLag days}
#'
#' @author Guillaume Evin
lag.trans.proba.mat.byClass = function(mat.prec, # matrix of precipitation
                                       isPeriod, # vector of logical
                                       th, # threshold above which we consider that a day is wet
                                       nlag, # number of lag days
                                       dayScale # time resolution
){
  # lag.trans.proba.proba (computeStat_lib)
  return(apply(mat.prec,2,lag.trans.proba.vector.byClass,isPeriod,th,nlag,dayScale))
}


#==============================================================================
#' modify.cor.matrix
#'
#' Modify a non-positive definite correlation matrix in order to have a positive 
#' definite matrix
#'
#' @param cor.matrix possibly non-positive definite correlation matrix
#'
#' @return positive definite correlation matrix
#'
#' @references Rousseeuw, P. J. and G. Molenberghs. 1993. Transformation of non positive semidefinite
#' correlation matrices. Communications in Statistics: Theory and Methods 22(4):965-984.
#' @references Rebonato, R., & Jackel, P. (2000). The most general methodology to create a valid 
#' correlation matrix for risk management and option pricing purposes. J. Risk, 2(2), 17-26.
#'
#' @author Guillaume Evin
#========================================================================
modify.cor.matrix = function(cor.matrix){
  # eigen values and vectors
  eigen.out = eigen(cor.matrix,symmetric=T)
  eig.val.diag = eigen.out$values
  eig.vec = eigen.out$vectors
  
  # is there negative eigen values?
  is.neg.eigen.val = eig.val.diag<0
  
  # if yes, replace this negative eigen values by small positive values
  if(any(is.neg.eigen.val)){
    eig.val.diag[is.neg.eigen.val] = 10^-10
    eig.val = diag(eig.val.diag)
    # recontruct correlation matrix
    cor.matrix.r = eig.vec %*% eig.val %*% solve(eig.vec)
    # return normalized correlation matrix
    cor.matrix.out = cor.matrix.r / (diag(cor.matrix.r) %*% t(diag(cor.matrix.r)))
  }else{
    # if there are no negative eigen values, return the original correlation matrix
    cor.matrix.out = cor.matrix
  }
  
  # diagonal values to be equal to 1, it appears that there is no tolerance for some functions
  # (generation of Student variates) and these values differ from 1 by an order of 10^-16
  diag(cor.matrix.out)=1
  return(cor.matrix.out)
}


#==============================================================================
#' get.df.Student
#'
#' Estimates the nu parameter (degrees of freedom) of the multivariate Student
#' distribution when the correlation matrix Sig is given 
#'
#' @param P matrix of non-zero precipitation (zero precipitation are set to NA)
#' @param Sig correlation matrix
#' @param max.df maximum degrees of freedom tested (default=20)
#'
#' @return nu estimate
#'
#' @references McNeil et al. (2005) "Quantitative Risk Management"
#'
#' @author Guillaume Evin
get.df.Student = function(P,Sig,max.df=20){
  # transformation PIT for the margins
  U = get.emp.cdf.matrix(P)
  
  # Compute likelihood for every df value, from 1 to 20, for more than 20, we can approximate by a Gaussian dist.
  vec.lk = rep(NA,max.df)
  for(df in 1:max.df){
    t.data = apply(U, 2, qt, df = df)
    lk = mvtnorm::dmvt(x=t.data,sigma=Sig,df=df,log=T) - apply(dt(x=t.data,df=df,log=T),1,sum)
    nz = !is.na(lk)
    vec.lk[df] = sum(lk[nz])
  }
  
  # We take the df value which maximizes the likelihood
  df.hat = which.max(vec.lk)
  
  # if the likelihood return finite values
  if(!is.infinite(vec.lk[df.hat])){
    return(df.hat)
  }else{
    return(max.df+10)
  }
}


#==============================================================================
#' get.df.Student
#'
#' get the cdf values (empirical distribution) of positive precipitation  
#'
#' @param X matrix of positive precipitation  
#' 
#' @return matrix with cdf values (NA if zero precipitation)
#'
#' @author Guillaume Evin
get.emp.cdf.matrix = function(X){
  # number of columns
  p = ncol(X)
  
  # prepare output
  Y = matrix(NA,nrow=nrow(X),ncol=p)
  
  # loop over the stations
  for(i in 1:p){
    # all data
    X.i = X[,i]
    # filter
    nz = !is.na(X.i)
    # cdf values (Gringorten prob (optimized for gumbel dist))
    X.cdf = lmomco::pp(X.i[nz],a=0.44,sort=F)
    # gaussian variates
    Y[which(nz),i] = X.cdf
  }
  
  return(Y)
}


#==============================================================================
#' get.list.month
#' 
#' return a vector of 3-char tags of the 12 months
get.list.month = function(){
  return(c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AOU','SEP','OCT','NOV','DEC'))
}


#==============================================================================
#' get.period.fitting.month
#' 
#' @param m.char 3-letter name of a month (e.g. 'JAN')
#' 
#' return the 3 indices corresponding to the 3-month period of a month ('JAN')
get.period.fitting.month = function(m.char){
  # m.char: 3-letter name of a month (e.g. 'JAN')
  
  # list of months
  list.m = get.list.month()
  
  # index of this month
  i.m = which(list.m==m.char)
  
  # return 3-month period corresponding to this month
  vec.m.ext = c(12,1:12,1)
  return(vec.m.ext[i.m:(i.m+2)])
}


#==============================================================================
#' get.list.season
#'
#' get the vector of the four seasons c('DJF','MAM','JJA','SON')
#'
#' @author Guillaume Evin
get.list.season = function(){
  return(c('DJF','MAM','JJA','SON'))
}


#==============================================================================
#' month2season
#'
#' transform vector of months to seasons
#'
#' @param vecMonth a vector of months given as integers 1:12
#'
#' @author Guillaume Evin
month2season = function(vecMonth){
  iSeason = c(1,1,2,2,2,3,3,3,4,4,4,1)
  return(iSeason[vecMonth])
}


#==============================================================================
#' EGPD.pGI, EGPD.dGI, EGPD.qGI
#'
#' First parametric family for G(v) = v^kappa: distribution, density and quantile function
#'
#' @param v probability
#' @param kappa transformation parameter greater than 0
#' @param p probability
#'
#' @return distribution, density and quantile of EGPD
#'
#' @author Guillaume Evin
#' @name functions.EGPD.GI
NULL

#' @rdname functions.EGPD.GI
EGPD.pGI = function(v,kappa) return(v^kappa)
#' @rdname functions.EGPD.GI
EGPD.dGI = function(v,kappa) return(kappa*v^(kappa-1))
#' @rdname functions.EGPD.GI
EGPD.qGI = function(p,kappa) return(p^(1/kappa))


#==============================================================================
#' dEGPD.GI, pEGPD.GI, qEGPD.GI, rEGPD.GI
#'
#' Density function, distribution function, quantile function, random generation
#' for the unified EGPD distribution
#'
#' @param x Vector of quantiles
#' @param p Vector of probabilities
#' @param n Number of observations
#' @param kappa transformation parameter greater than 0
#' @param sig Scale parameter
#' @param xi Shape parameter
#'
#' @return dEGPD.GI gives the density function, pEGPD.GI gives the distribution function, 
#' qEGPD.GI gives the quantile function, and rEGPD.GI generates random deviates.
#'
#' @author Guillaume Evin
#' @name dist.functions.EGPD.GI
NULL

#' @rdname dist.functions.EGPD.GI
dEGPD.GI = function(x,kappa,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  dH = Renext::dGPD(x=x,scale=sig,shape=xi)
  dEGPD = dH*EGPD.dGI(pH,kappa)
  return(dEGPD)
}
#' @rdname dist.functions.EGPD.GI
pEGPD.GI = function(x,kappa,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  return(EGPD.pGI(pH,kappa))
}
#' @rdname dist.functions.EGPD.GI
qEGPD.GI = function(p,kappa,sig,xi){
  qG = EGPD.qGI(p,kappa)
  qH = Renext::qGPD(p=qG,scale=sig,shape=xi)
  return(qH)
}
#' @rdname dist.functions.EGPD.GI
rEGPD.GI = function(n,kappa,sig,xi){
  u = runif(n=n)
  rEGPD = qEGPD.GI(u,kappa,sig,xi)
  return(rEGPD)
}


#==============================================================================
#' EGPD.GI.mu0, EGPD.GI.mu1, EGPD.GI.mu2
#'
#' Probability Weighted Moments of order 0, 1 and 2 of the unified EGPD distribution
#'
#' @param kappa transformation parameter greater than 0
#' @param sig Scale parameter
#' @param xi Shape parameter
#'
#' @return Probability Weighted Moments
#'
#' @author Guillaume Evin
#' @name PWM.EGPD.GI
NULL

#' @rdname PWM.EGPD.GI
EGPD.GI.mu0 = function(kappa,sig,xi){
  mu0 = (sig/xi)*(kappa*beta(kappa,1-xi)-1)
  return(mu0)
}

#' @rdname PWM.EGPD.GI
EGPD.GI.mu1 = function(kappa,sig,xi){
  mu1 = (sig/xi)*(kappa*(beta(kappa,1-xi)-beta(2*kappa,1-xi))-1/2)
  return(mu1)
}

#' @rdname PWM.EGPD.GI
EGPD.GI.mu2 = function(kappa,sig,xi){
  mu2 = (sig/xi)*(kappa*(beta(kappa,1-xi)-2*beta(2*kappa,1-xi)+beta(3*kappa,1-xi))-1/3)
  return(mu2)
}


#==============================================================================
#' EGPD.GI.fPWM
#'
#' Parameter estimation of the unified EGPD distribution with the PWM method.
#' Set of equations which have to be equal to zero
#'
#' @param par vector of parameters kappa,sig,xi
#' @param pwm set of probability weighted moments of order 0, 1 and 2
#' @param xi shape parameter
#'
#' @return differences between expected and target weighted moments
#'
#' @author Guillaume Evin
# set of equations which have to be equal to zero
EGPD.GI.fPWM =  function(par,pwm,xi){
  kappa = par[1]
  sig = par[2]
  
  y = numeric(2)
  y[1] = EGPD.GI.mu0(kappa,sig,xi) - pwm[1]
  y[2] = EGPD.GI.mu1(kappa,sig,xi) - pwm[2]
  
  return(y)
}

#==============================================================================
#' EGPD.GI.fit.PWM
#'
#' Parameter estimation of the unified EGPD distribution with the PWM method.
#' Numerical solver of the system of nonlinear equations
#'
#' @param x vector of parameters kappa,sig
#' @param xi shape parameter
#'
#' @return estimated parameters kappa, sig, xi
#'
#' @author Guillaume Evin
EGPD.GI.fit.PWM = function(x,xi=0.05){
  sam.pwm = c(EnvStats::pwMoment(x,k=0),
              EnvStats::pwMoment(x,k=1))
  fit.out = nleqslv::nleqslv(c(2,sd(x)),EGPD.GI.fPWM,jac=NULL,pwm=sam.pwm,xi=xi)
  fit.out$x = c(fit.out$x,xi)
  return(fit.out)
}

#==============================================================================
#' agg.matrix
#'
#' Simple accumulation of a matrix of precipitation
#'
#' @param mat matrix nDates x nStations to be aggregated
#' @param k number of days for the accumulation
#' @param average logical: should we average over the different periods (default=F)
#'
#' @return aggregated matrix
#'
#' @author Guillaume Evin
agg.matrix = function(mat,k,average=F){
  # number of steps
  n = floor(nrow(mat)/k)
  mat.agg = matrix(0,nrow=n,ncol=ncol(mat))
  
  # aggregate
  i0 = (seq(from=1,to=n)-1)*k
  for(i in 1:k){
    mat.agg = mat.agg + mat[i0+i,]
  }
  
  # if we wnt the average on the period
  if(average) mat.agg = mat.agg/k
  
  return(mat.agg)
}

#==============================================================================
#' get.listOption
#'
#' get default options and check values proposed by the user 
#'
#' @param listOption list containing fields corr. to the different options. Can be NULL if no options are set
#'
#' @return \item{listOption}{list of options}
#'
#' @author Guillaume Evin
get.listOption = function(listOption){
  
  ######### Check inputs and assign default values ##########
  if(is.null(listOption)){
    listOption = list()
  }
  
  #_____ options for occurrence
  
  # th
  if('th' %in% names(listOption)){
    th = listOption[['th']]
    if(!(is.numeric(th)&(th>=0))) stop('wrong value for th')
  }else{
    th = 0.2
    listOption[['th']] = th
  }
  
  # nLag
  if('nLag' %in% names(listOption)){
    nLag = listOption[['nLag']]
    if(!any(nLag==c(1,2,3,4,5))) stop('nLag must be between 1 and 5')
  }else{
    nLag = 2
    listOption[['nLag']] = nLag
  }
  
  
  #_____ options for amount
  
  # typeMargin
  if('typeMargin' %in% names(listOption)){
    typeMargin = listOption[['typeMargin']]
    if(!(typeMargin %in% c('mixExp','EGPD'))) stop('typeMargin ust be equal to mixExp or EGPD')
  }else{
    typeMargin = 'EGPD'
    listOption[['typeMargin']] = typeMargin
  }
  
  # copulaInt
  if('copulaInt' %in% names(listOption)){
    copulaInt = listOption[['copulaInt']]
    if(!copulaInt%in%c('Gaussian','Student')) stop('copulaInt must be equal to Gaussian or Student')
  }else{
    copulaInt = 'Student'
    listOption[['copulaInt']] = copulaInt
  }
  
  # isMAR
  if('isMAR' %in% names(listOption)){
    isMAR = listOption[['isMAR']]
    if(!is.logical(isMAR)) stop('isMAR must be logical')
  }else{
    isMAR = T
    listOption[['isMAR']] = isMAR
  }
  
  # is3Damount
  if('is3Damount' %in% names(listOption)){
    is3Damount = listOption[['is3Damount']]
    if(!is.logical(is3Damount)) stop('is3Damount must be logical')
  }else{
    is3Damount = F
    listOption[['is3Damount']] = is3Damount
  }
  
  # nChainFit
  if('nChainFit' %in% names(listOption)){
    nChainFit = listOption[['nChainFit']]
    if(!is.numeric(nChainFit)) stop('nChainFit must be an integer')
  }else{
    nChainFit = 100000
    listOption[['nChainFit']] = nChainFit
  }
  
  # nCluster
  if('nCluster' %in% names(listOption)){
    nCluster = listOption[['nCluster']]
    if(!is.numeric(nCluster)) stop('nCluster must be an integer')
  }else{
    nCluster = 1
    listOption[['nCluster']] = nCluster
  }
  
  # isParallel: parallelize results if nCluster is greater than 1
  isParallel = nCluster>1
  listOption[['isParallel']] = isParallel
  
  # return
  return(listOption)
}

#==============================================================================
#' infer.mat.omega
#'
#' find omega correlation leading to estimates cor between occurrences
#'
#' @param P.mat matrix of precipitation n x p
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param nLag order of the Markov chain
#' @param pr.state output of function \code{\link{lag.trans.proba.mat.byClass}}
#' @param nChainFit length of the simulated chains used during the fitting
#' @param isParallel logical: indicate computation in parallel or not (easier for debugging)
#'
#' @return A list with different objects 
#' \itemize{
#'   \item \strong{Qtrans.mat}: matrix nStation x n.comb of transition probabilites
#'   \item \strong{mat.comb}: matrix of possible combination n.comb x nLag
#'   \item \strong{mat.omega}: The spatial correlation matrix of occurrences \eqn{\Omega} (see Evin et al., 2018).
#' }
#'
#' @author Guillaume Evin
infer.mat.omega = function(P.mat,isPeriod,th,nLag,pr.state,nChainFit,isParallel){
  # class = Period
  isClass = isPeriod
  
  # filtered matrix of precipitation for this period (3-mon moving window by dafault)
  P.mat.class = P.mat[isClass,]
  
  # observed correlation of dry/wet states (see Eq. 6 in Evin et al. 2018)
  pi0 = dry.day.frequency(P.mat.class,th)
  pi1 = wet.day.frequency(P.mat.class,th)
  pi.occ = joint.proba.occ(P.mat.class,th)
  cor.occ.obs = cor.obs.occ(pi.occ$p00,pi0,pi1)
  
  # number of possible transitions
  n.comb = 2^nLag
  
  # a matrix of possible combination n.comb x nLag
  mat.comb = as.matrix(pr.state[[1]][,1:nLag])
  
  #####   #####   #####   #####   ##### 
  ##### prob / normal quantiles of transitions
  #####   #####   #####   #####   ##### 
  
  # retrieve the last column (pobabilities) for all stations
  # return a list with one field per station containing a vector n.comb x 1
  # of transition probabilites
  Ptrans.list = lapply(pr.state,'[',nLag+1)
  
  # transform to normal quantiless
  Qtrans.list = lapply(Ptrans.list,function(x) qnorm(unlist(x)))
  
  # reshape in a matrix nStation x n.comb
  Qtrans.mat = matrix(unlist(Qtrans.list), ncol=n.comb, byrow=T)
  
  # filter infinite values if Pr = 0 or 1
  Qtrans.mat[Qtrans.mat==-Inf] = -10^5
  Qtrans.mat[Qtrans.mat==Inf] = 10^5
  
  ############################## 
  # estimation of omega matrices
  ############################## 
  mat.omega = get.mat.omega(cor.occ.obs,Qtrans.mat,mat.comb,nLag,nChainFit,isParallel)
  mat.omega = modify.cor.matrix(mat.omega)
  
  return(list(Qtrans.mat=Qtrans.mat,mat.comb=mat.comb,mat.omega=mat.omega))
}


#==============================================================================
#' get.mat.omega
#'
#' find omega correlation leading to estimates cor between occurrences
#'
#' @param cor.obs matrix p x p of observed correlations between occurrences for all pairs of stations
#' @param Qtrans.mat transition probabilities, 2 x ncomb matrix
#' @param mat.comb matrix of logical: ncomb x nlag
#' @param nLag order of the Markov chain
#' @param nChainFit length of the simulated chains used during the fitting
#' @param isParallel logical: indicate computation in parallel or not (easier for debugging)
#'
#' @return \item{matrix}{omega correlations for all pairs of stations}
#'
#' @author Guillaume Evin
get.mat.omega = function(cor.obs,Qtrans.mat,mat.comb,nLag,nChainFit,isParallel){
  # number of stations
  p = ncol(cor.obs)
  
  # possible pairs
  vec.pairs = combn(1:p, 2)
  n.pairs = ncol(vec.pairs)
  
  # apply find.omega for each pair of stations
  if(isParallel){
    omega.paral = foreach(i.pair=1:n.pairs, .combine='cbind') %dopar% {
      i = vec.pairs[1,i.pair]
      j = vec.pairs[2,i.pair]
      return(find.omega(cor.obs[i,j],Qtrans.mat[c(i,j),],mat.comb,nLag,nChainFit))
    }
  }else{
    omega.paral = vector(length=n.pairs)
    for(i.pair in 1:n.pairs){
      i = vec.pairs[1,i.pair]
      j = vec.pairs[2,i.pair]
      
      omega.paral[i.pair] = find.omega(cor.obs[i,j],Qtrans.mat[c(i,j),],mat.comb,nLag,nChainFit)
    }
  }
  
  # prepare results
  mat.omega = matrix(0,nrow=p,ncol=p)
  diag(mat.omega) = 1
  
  # fill matrix
  for(i.pair in 1:n.pairs){
    i = vec.pairs[1,i.pair]
    j = vec.pairs[2,i.pair]
    mat.omega[i,j] = omega.paral[i.pair]
    mat.omega[j,i] = omega.paral[i.pair]
  }
  
  return(mat.omega)
}




#==============================================================================
#' find.omega
#'
#' finds the correlation between normal variates leading to correlation between occurrences
#'
#' @param rho.emp target correlation between occurences
#' @param Qtrans.mat transition probabilities, 2 x ncomb matrix
#' @param mat.comb matrix of logical: ncomb x nlag
#' @param nLag order of the Markov chain
#' @param nChainFit length of the simulated chains used during the fitting
#'
#' @return \item{scalar}{needed correlation}
#'
#' @author Guillaume Evin
find.omega = function(rho.emp,Qtrans.mat,mat.comb,nLag,nChainFit){
  # f.inf and f.sup represent resp. the minimum and maximum differences
  # which can be obtained with the simulated correlation and the empirical
  # correlation
  f.inf = cor.emp.occ(-1,Qtrans.mat,mat.comb,nLag,nChainFit) - rho.emp
  f.sup = cor.emp.occ(1,Qtrans.mat,mat.comb,nLag,nChainFit) - rho.emp
  
  # if f.sup<=0, it means that even with the a max omega value of 1, we cannot reach
  # the empirical correlation (negative difference between simulated cor and emp. cor)
  # we simply return a max possible omega value (less than 1 since 1 correlation leads to numerical pb)
  if(f.sup<=0){
    return(1)
  }else if(f.inf>=0){
    # if f.inf>=0, it means that even with the a min omega value of -1, we cannot reach
    # the empirical correlation (positive difference between simulated cor and emp. cor)
    # we simply return a 0 value since negative values are not physically plausible
    # it can happen when emp correlation are estimated on short series (uncertain estimate) 
    return(0)
  }else{
    if(rho.emp<0){
      # if the empirical correlation is negative, we simply return a 0 value since negative 
      # values are not physically plausible
      # it can happen when emp correlation are estimated on short series (uncertain estimate) 
      return(0)
    }else{
      # else, we find omega that leads to rho.emp (f should be zero for this omega)
      f = function(w){
        cor.emp.occ(w,Qtrans.mat,mat.comb,nLag,nChainFit) - rho.emp
      }
      return(uniroot(f,c(rho.emp,1),extendInt="upX", tol = 1e-3)$root)
    }
  }
}



#==============================================================================
#' cor.emp.occ
#'
#' Finds observed correlations between occurrences corresponding
#' to a degree of correlation of Gaussian multivariate random numbers
#' @useDynLib GWEX, .registration = TRUE
#' @noRd
#' @param w correlation of Gaussian multivariates
#' @param Qtrans.mat transition probabilities, 2 x ncomb matrix
#' @param mat.comb matrix of logical: ncomb x nlag
#' @param nLag order of the Markov chain
#' @param nChainFit number of simulated variates
#' @param myseed seed of random variates
#'
#' @return \item{scalar}{correlation between occurrences}
#'
#' @author Guillaume Evin
cor.emp.occ = function(w,Qtrans.mat,mat.comb,nLag,nChainFit,myseed=1){
  # set seed of random generation
  set.seed(myseed)
  
  # genere gaussienne multivariee
  w.mat = rbind(c(1,w),c(w,1))
  rndNorm = MASS::mvrnorm(nChainFit,rep(0,2),w.mat)
  
  # empirical correlation with the Markov process
  matComb.f = mat.comb*1
  storage.mode(matComb.f) <- "integer"
  
  cor.emp = .Fortran("corMarkovChain", rndNorm=rndNorm, QtransMat=Qtrans.mat, PACKAGE="GWEX",
                     matComb=matComb.f, n=as.integer(nChainFit), nLag=as.integer(nLag), r=double(1))$r
  
  # correlations occurrences
  return(cor.emp)
}


#==============================================================================
#' get.M0
#'
#' find matrix of correlations leading to estimates cor between intensities
#'
#' @param cor.obs matrix p x p of observed correlations between intensities for all pairs of stations
#' @param infer.mat.omega.out output of \code{\link{infer.mat.omega}}
#' @param nLag order of the Markov chain
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#' @param nChainFit integer indicating the length of simulated chains
#' @param isParallel logical: indicate computation in parallel or not (easier for debugging)
#' 
#' @return list with two items
#' \itemize{
#'   \item \strong{Xt}{long simulation of the wet/dry states according to the model}
#'   \item \strong{M0}{omega correlations for all pairs of stations}
#' }
#'
#' @author Guillaume Evin
get.M0 = function(cor.obs,infer.mat.omega.out,nLag,parMargin,typeMargin,nChainFit,isParallel){
  
  # retrieve outputs from infer.mat.omega.out
  Qtrans.mat = infer.mat.omega.out$Qtrans.mat
  mat.comb = infer.mat.omega.out$mat.comb
  mat.omega = infer.mat.omega.out$mat.omega
  
  # number of stations
  p = ncol(cor.obs)
  
  # number of combinations
  n.comb = ncol(Qtrans.mat)
  
  # simulate occurrence (once only)
  Qtrans = QtransMat2Array(nChainFit,p,Qtrans.mat)
  Xt = array(0,dim=c(nChainFit,p))
  rndNorm = MASS::mvrnorm(nChainFit,rep(0,p),mat.omega)
  for(t in (nLag+1):nChainFit){
    for(st in 1:p){
      comb = Xt[(t-nLag):(t-1),st]==1
      i.comb = which(apply(mat.comb, 1, function(x) all(x==comb)))
      qTr = Qtrans[t,st,i.comb]
      Xt[t,st] = (rndNorm[t,st]<=qTr)
    }
  }
  
  # possible pairs
  vec.pairs = combn(1:p, 2)
  n.pairs = ncol(vec.pairs)
  
  # apply find.zeta for each pair of stations
  if(isParallel){
    zeta.paral = foreach(i.pair=1:n.pairs, .combine='cbind') %dopar% {
      i = vec.pairs[1,i.pair]
      j = vec.pairs[2,i.pair]
      return(find.zeta(cor.obs[i,j],nChainFit,Xt[,c(i,j)],parMargin[c(i,j),],typeMargin))
    }
  }else{
    zeta.paral = vector(length=n.pairs)
    for(i.pair in 1:n.pairs){
      i = vec.pairs[1,i.pair]
      j = vec.pairs[2,i.pair]
      zeta.paral[i.pair] = find.zeta(cor.obs[i,j],nChainFit,Xt[,c(i,j)],parMargin[c(i,j),],typeMargin)
    }
  }
  
  # prepare results
  M0 = matrix(1,nrow=p,ncol=p)
  diag(M0) = 1
  
  # fill matrix
  for(i.pair in 1:n.pairs){
    i = vec.pairs[1,i.pair]
    j = vec.pairs[2,i.pair]
    M0[i,j] = zeta.paral[i.pair]
    M0[j,i] = zeta.paral[i.pair]
  }
  
  return(list(M0=M0,Xt=Xt))
}


#==============================================================================
#' find.zeta
#'
#' finds the correlation between normal variates leading to correlation between intensities
#'
#' @param eta.emp target correlation between intensities
#' @param nChainFit number of simulations
#' @param Xt simulated occurrences, n x 2 matrix
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#'
#' @return \item{scalar}{needed correlation}
#'
#' @author Guillaume Evin
find.zeta = function(eta.emp,nChainFit,Xt,parMargin,typeMargin){
  # f.inf and f.sup represent resp. the minimum and maximum differences
  # which can be obtained with the simulated correlation and the empirical
  # correlation
  f.inf = cor.emp.int(0,nChainFit,Xt,parMargin,typeMargin) - eta.emp
  f.sup = cor.emp.int(1,nChainFit,Xt,parMargin,typeMargin) - eta.emp
  
  # if f.sup<=0, it means that even with the a max zeta value of 1, we cannot reach
  # the empirical correlation (negative difference between simulated cor and emp. cor)
  # we simply return a max possible zeta value (less than 1 since 1 correlation leads to numerical pb)
  if(f.sup<=0){
    return(0.99999)
  }else if(f.inf>=0){
    # if f.inf>=0, it means that even with the a min zeta value of -1, we cannot reach
    # the empirical correlation (positive difference between simulated cor and emp. cor)
    # we simply return a 0 value since negative values are not physically plausible
    # it can happen when emp correlation are estimated on short series (uncertain estimate) 
    return(0)
  }else{
    if(eta.emp<0){
      # if the empirical correlation is negative, we simply return a 0 value since negative 
      # values are not physically plausible
      # it can happen when emp correlation are estimated on short series (uncertain estimate) 
      return(0)
    }else{
      # else, we find zeta that leads to eta.emp (f should be zero for this zeta)
      f = function(w){
        cor.emp.int(w,nChainFit,Xt,parMargin,typeMargin) - eta.emp
      }
      zeta = uniroot(f,c(0,1),extendInt="upX",tol = 1e-3)$root
      return(zeta)
    }
  }
}


#==============================================================================
#' cor.emp.int
#'
#' Finds observed correlations between intensities corresponding
#' to a degree of correlation of Gaussian multivariate random numbers
#' @param zeta correlation of Gaussian multivariates
#' @param nChainFit number of simulated variates
#' @param Xt simulated occurrences, n x 2 matrix
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#'
#' @return \item{scalar}{correlation between simulated intensities}
#'
#' @author Guillaume Evin
cor.emp.int = function(zeta,nChainFit,Xt,parMargin,typeMargin){
  # generate the same random numbers if zeta=1, which avoids problems at the bounds (minimzation diverving, or
  # function having the same signs at both bounds if cor.emp.int(1,...) is slightly less than 0 if zeta=1)
  set.seed(1)
  # Simulation of the spatial and temporal dependence between amounts
  Yt.Pr = pnorm(MASS::mvrnorm(n=nChainFit, mu=rep(0,2), Sigma=matrix(c(1,zeta,zeta,1),2,2)))
  # We obtain directly related intensities
  Yt = array(0,dim=c(nChainFit,2))
  for(i.st in 1:2) Yt[,i.st] = unif.to.prec(parMargin[i.st,],typeMargin,Yt.Pr[,i.st])
  # Mask with occurrences
  Yt[Xt==0] = 0
  # cor.Pearson.cor
  cor.Pearson = cor(Yt[,1],Yt[,2])
  return(cor.Pearson)
}


#==============================================================================
#' get.vec.autocor
#'
#' find rho autocorrelation leading to empirical estimates
#'
#' @param vec.ar1.obs vector of observed autocorrelations for all stations
#' @param Xt simulated occurrences given model parameters of wet/dry states
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#' @param nChainFit integer indicating the length of the simulated chains
#' @param isParallel logical: indicate computation in parallel or not (easier for debugging)
#'
#' @return \item{vector}{vector of rho parameters to simulate the MAR process}
#'
#' @author Guillaume Evin
get.vec.autocor = function(vec.ar1.obs,Xt,parMargin,typeMargin,nChainFit,isParallel){
  # number of stations
  p = length(vec.ar1.obs)
  
  # find autocor leading to obs autocor
  if(isParallel){
    iSt = NULL
    vec.autocor = foreach(iSt=1:p, .combine='cbind') %dopar% {
      return(find.autocor(vec.ar1.obs[iSt],nChainFit,Xt[,iSt],parMargin[iSt,],typeMargin))
    }
    vec.autocor = as.vector(vec.autocor)
  }else{
    vec.autocor = vector(length = p)
    for(iSt in 1:p){
      vec.autocor[iSt] = find.autocor(vec.ar1.obs[iSt],nChainFit,Xt[,iSt],parMargin[iSt,],typeMargin)
    }
  }
  
  return(vec.autocor)
}

#==============================================================================
#' find.autocor
#'
#' finds the autocorrelation leading to observed autocorrelation
#'
#' @param autocor.emp target correlation between intensities
#' @param nChainFit number of simulations
#' @param Xt simulated occurrences, nChainFit x 2 matrix
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#'
#' @return \item{scalar}{needed correlation}
#'
#' @author Guillaume Evin
find.autocor = function(autocor.emp,nChainFit,Xt,parMargin,typeMargin){
  f1 = autocor.emp.int(0.95,nChainFit,Xt,parMargin,typeMargin) - autocor.emp
  if(f1<=0){
    return(0.95)
  }else{
    f0 = autocor.emp.int(0,nChainFit,Xt,parMargin,typeMargin) - autocor.emp
    if(f0>=0){
      return(0)
    }else{
      f = function(w){
        autocor.emp.int(w,nChainFit,Xt,parMargin,typeMargin) - autocor.emp
      }
      return(uniroot(f,c(0,0.95),tol = 1e-3)$root)
    }
  }
}

#==============================================================================
#' autocor.emp.int
#'
#' Finds empirical autocorrelations (lag-1) between intensities corresponding
#' to a degree of autocorrelation of an AR(1) process
#' @param rho autocorrelation of the AR(1) process
#' @param nChainFit number of simulated variates
#' @param Xt simulated occurrences, nChainFit x 2 matrix
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#' 
#' @return \item{scalar}{correlation between simulated intensities}
#'
#' @author Guillaume Evin
autocor.emp.int = function(rho,nChainFit,Xt,parMargin,typeMargin){
  # control random seed
  set.seed(1)
  # Simulation from an AR(1) process
  if(rho==0){
     Yt.AR1=rnorm(n=nChainFit)
  }else{
    Yt.AR1=stats::arima.sim(n=nChainFit,model=list(ar=rho),rand.gen=rnorm, sd=sqrt(1-rho^2))
  }
  # to proba
  Yt.Pr = pnorm(Yt.AR1)
  # Related intensities
  Yt = unif.to.prec(parMargin,typeMargin,Yt.Pr)
  # Mask with occurrences
  Yt[Xt==0] = 0
  # Resulting cor
  cor.Pearson = cor(Yt[1:(nChainFit-1)],Yt[2:nChainFit])
  return(cor.Pearson)
}

#==============================================================================
#' QtransMat2Array
#'
#' reshape Qtrans.mat to an array
#'
#' @param n matrix of precipitation
#' @param p number of stations
#' @param Qtrans.mat  transition probabilities, 2 x ncomb matrix
#'
#' @return \item{array}{array of transition probabilities with dimension n x p x n.comb}
#'
#' @author Guillaume Evin
QtransMat2Array = function(n,p,Qtrans.mat){
  n.comb = ncol(Qtrans.mat)
  Qtrans = array(0,dim=c(n,p,n.comb))
  for(i.st in 1:p){
    for(i.comb in 1:n.comb){
      # find corresponding days for this combination of dry/wet states
      ind.mat = cbind(1:n,rep(i.st,n),i.comb)
      # fill array of transitions
      Qtrans[ind.mat] = Qtrans.mat[i.st,i.comb]
    }
  }
  return(Qtrans)
}


#==============================================================================
#' joint.proba.occ
#'
#' joint probabilities of occurrences for all pairs of stations
#'
#' @param P matrix of precipitation
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#'
#' @return \item{list}{list of joint probabilities}
#'
#' @author Guillaume Evin
joint.proba.occ = function(P,th){
  # number of stations
  p = ncol(P)
  
  # prepare results
  p00 = p10 = p01 = p11 = matrix(1,nrow=p,ncol=p)
  
  for(i in 1:(p-1)){
    ri = P[,i]
    for(j in (i+1):p){
      rj = P[,j]
      # no nan
      nz = (!is.na(ri)) & (!is.na(rj))
      
      # dry-dry
      p00ij = mean(ri[nz]<=th&rj[nz]<=th)
      p00[i,j] = p00ij
      p00[j,i] = p00ij
      
      # dry-wet
      p01ij = mean(ri[nz]<=th&rj[nz]>th)
      p01[i,j] = p01ij
      p10[j,i] = p01ij
      
      # wet-dry
      p10ij = mean(ri[nz]>th&rj[nz]<=th)
      p10[i,j] = p10ij
      p01[j,i] = p10ij
      
      # wet-wet
      p11ij = mean(ri[nz]>th&rj[nz]>th)
      p11[i,j] = p11ij
      p11[j,i] = p11ij
    }
  }
  
  # joint probabilities
  return(list(p00=p00,p01=p01,p10=p10,p11=p11))
}


#==============================================================================
#' cor.obs.occ
#'
#' provide observed correlations between occurrences for all pairs of stations
#' see Mhanna et al. (2012)
#' 
#' @references Mhanna, Muamaraldin, and Willy Bauwens. “A Stochastic Space-Time Model 
#' for the Generation of Daily Rainfall in the Gaza Strip.” International Journal of 
#' Climatology 32, no. 7 (June 15, 2012): 1098–1112. doi:10.1002/joc.2305.
#'
#' @param pi00 joint probability of having dry states
#' @param pi0 probability of having a dry state
#' @param pi1 probability of having a wet state
#'
#' @return \item{scalar}{matrix of observed correlations}
#'
#' @author Guillaume Evin
cor.obs.occ = function(pi00,pi0,pi1){
  # number of stations
  p = ncol(pi00)
  
  # prepare results
  cor.obs = matrix(1,nrow=p,ncol=p)
  diag(cor.obs) = 1
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      # Eq 6 de Mhanna (2012)
      sigi = sqrt(pi0[i]*pi1[i])
      sigj = sqrt(pi0[j]*pi1[j])
      cor.ij =  (pi00[i,j]-pi0[i]*pi0[j])/(sigi*sigj)
      
      # Assign values to matrix
      cor.obs[i,j] = cor.ij
      cor.obs[j,i] = cor.ij
    }
  }
  
  # joint probabilities
  return(cor.obs)
}


#==============================================================================
#' infer.dep.amount
#'
#' estimate parameters which control the spatial dependence between intensities using a copula
#' @param P.mat precipitation matrix
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param infer.mat.omega.out output of \code{\link{infer.mat.omega}}
#' @param nLag order of he Markov chain for the transitions between dry and wet states (=2 by default)
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param parMargin parameters of the margins 2 x 3
#' @param typeMargin 'EGPD' (Extended GPD) or 'mixExp' (Mixture of Exponentials). 'EGPD' by default
#' @param nChainFit integer, length of the runs used during the fitting procedure. =100000 by default
#' @param isMAR logical value, do we apply a Autoregressive Multivariate Autoregressive model (order 1) =TRUE by default
#' @param copulaInt 'Gaussian' or 'Student': type of dependence for amounts (='Student' by default)
#' @param isParallel logical: indicate computation in parallel or not (easier for debugging)
#' 
#' @return \item{list}{list of estimates (e.g., M0, dfStudent)}
#'
#' @author Guillaume Evin
infer.dep.amount = function(P.mat,isPeriod,infer.mat.omega.out,nLag,th,parMargin,typeMargin,nChainFit,
                            isMAR,copulaInt,isParallel){
  # Class = Period
  isClass = isPeriod
  
  # number of stations
  p = ncol(P.mat)
  n = nrow(P.mat)
  
  # filtered data for this class: replace small values below the threshold by zero
  # for the sake of comparison with simulated values
  P.mat.class = P.mat[isClass,]
  is.Zero = P.mat.class<=th&!is.na(P.mat.class)
  P.mat.class[is.Zero] = 0
  
  # direct pairwise Pearson correlations
  cor.int = cor(P.mat.class, use="pairwise.complete.obs")

  
  #### find corresponding needed zeta correlations between simulated intensities
  get.M0.out = get.M0(cor.int,infer.mat.omega.out,nLag,parMargin,typeMargin,nChainFit,isParallel)
  M0 = modify.cor.matrix(get.M0.out$M0)
  
  # find remaining parameters if necessary
  if(isMAR){
    ### MAR is a Multivariate Autoregressive Process of order 1
    # Here we model the temporal dep. by a simple diagonal matrix M1
    # where the diagonal elements are lag-1 autocorrelation
    vec.ar1.obs = vector(length=p)
    
    # on the contrary to GWEX, here I consider all days of the class (where there might
    # not be many elements) and the corr. previous days in order to increase the possible
    # days (instead of consider only days which have also previous days within the same class)
    indClassCur = which(isClass)
    if(indClassCur[1]==1) indClassCur = indClassCur[2:length(indClassCur)]
    indClassPrev = indClassCur-1
    
    # loop over the stations
    for(i.st in 1:p){
      P.st = P.mat[,i.st]
      
      # on the contrary to GWEX, here I consider all days of the class (where there might
      # not be many elements) and the corr. previous days in order to increase the possible
      # days (instead of consider only days which have also previous days within the same class)
      vec.ar1.obs[i.st] = cor(P.st[indClassPrev],P.st[indClassCur],use="pairwise.complete.obs")
    }
    
    # find corresponding ar(1) parameter
    vec.ar1 = get.vec.autocor(vec.ar1.obs,get.M0.out$Xt,parMargin,typeMargin,nChainFit,isParallel)
    
    # apply a MAR(1) process: multivariate autoregressive process, order 1
    par.dep.amount = fit.MAR1.amount(P.mat,isPeriod,th,copulaInt,M0=M0,A=diag(vec.ar1))
  }else{
    # spatial process
    par.dep.amount = fit.copula.amount(P.mat,isPeriod,th,copulaInt,M0)
  }
}


#==============================================================================
#' fit.copula.amount
#'
#' estimate parameters which control the spatial dependence between intensities using a copula
#' @param P.mat precipitation matrix
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param copulaInt type of dependence between inter-site amounts: 'Gaussian' or 'Student'
#' @param M0 Matrix containing the inter-site spatial correlations
#'
#' @return \item{list}{list of estimates (e.g., M0, dfStudent)}
#'
#' @author Guillaume Evin
fit.copula.amount = function(P.mat,isPeriod,th,copulaInt,M0)
{
  # Class = Period
  isClass = isPeriod
  
  # number of stations
  p = ncol(P.mat)
  n = nrow(P.mat)
  
  if(copulaInt=='Gaussian'){
    # For 'Gaussian', we only compute the inter-site correlations
    return(list(M0=M0))
    
  }else if(copulaInt=='Student'){
    # on the contrary to GWEX, here I consider all days of the class (where there might
    # not be many elements) and the corr. previous days in order to increase the possible
    # days (instead of consider only days which have also previous days within the same class)
    indNotClassCur = which(!isClass[2:n])
    indNotClassPrev = indNotClassCur-1
    
    # build matrix for this class and set small values <=th to NA
    P.class = P.mat[isClass,]
    is.small = P.class<=th&!is.na(P.class)
    P.class[is.small] = NA
    
    # For 'Student', we also compute the degree of freedom of the Student copula
    dfStudent = get.df.Student(P.class,M0)
    # return results
    return(list(M0=M0,dfStudent=dfStudent))
    
  }else{
    stop('fit.copula.amount: unknown value for copulaInt: must be Gaussian or Student')
  }
}


#==============================================================================
#' fit.MAR1.amount
#'
#' estimate parameters which control the dependence between intensities with a
#' MAR(1) process
#' 
#' @references Matalas, N. C. 1967. “Mathematical Assessment of Synthetic Hydrology.” Water Resources 
#' Research 3 (4): 937–45. https://doi.org/10.1029/WR003i004p00937.
#' @references Bárdossy, A., and G. G. S. Pegram. 2009. “Copula Based Multisite Model for Daily 
#' Precipitation Simulation.” Hydrology and Earth System Sciences 13 (12): 2299–2314. https://doi.org/10.5194/hess-13-2299-2009.
#' 
#' @param P.mat precipitation matrix
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param copulaInt type of dependance between inter-site amounts: 'Gaussian' or 'Student'
#' @param M0 Matrix containing the inter-site spatial correlations
#' @param A Matrix containing the autocorrelation (temporal) correlations
#'
#' @return list with the following items
#' \itemize{
#'   \item \strong{M0}{omega correlations for all pairs of stations}
#'   \item \strong{A}{omega correlations for all pairs of stations}
#'   \item \strong{covZ}{covariance matrix of the MAR(1) process}
#'   \item \strong{sdZ}{standard deviation of the diagonal elements}
#'   \item \strong{corZ}{correlation matrix of the MAR(1) process}
#'   \item \strong{dfStudent}{degrees of freedom for the Student copula if CopulaInt is equal to "Student"}
#' }
#'
#' @author Guillaume Evin
fit.MAR1.amount = function(P.mat,isPeriod,th,copulaInt,M0,A){
  # Class = Period
  isClass = isPeriod
  
  # dimensions
  p = ncol(P.mat)
  n = nrow(P.mat)
  
  # covariance matrices of the MAR(1) process (Matalas, 1967; Bardossy & Pegram, 2009)
  covZ = M0 - A%*%t(M0)%*%t(A)
  
  # standard deviation
  sdZ = sqrt(diag(covZ))
  
  # correlation matrix
  corZ = modify.cor.matrix(cov2cor(covZ))
  
  if(copulaInt=='Gaussian'){
    # For 'Gaussian', we only compute the inter-site correlations
    return(list(M0=M0,A=A,covZ=covZ,sdZ=sdZ,corZ=corZ))
  }else if(copulaInt=='Student'){
    # on the contrary to GWEX, here I consider all days of the class (where there might
    # not be many elements) and the corr. previous days in order to increase the possible
    # days (instead of consider only days which have also previous days within the same class)
    indNotClassCur = which(!isClass[2:n])
    indNotClassPrev = indNotClassCur-1
    
    # build matrix where days not concerned by this class are set to NA
    P.class.nonzero = P.mat
    P.class.nonzero[indNotClassCur,] = NA
    P.class.nonzero[indNotClassPrev,] = NA
    
    # only positive values are retained
    is.small = P.class.nonzero<=th&!is.na(P.class.nonzero)
    P.class.nonzero[is.small] = NA
    
    # transform margins
    P.cdf.emp = get.emp.cdf.matrix(P.class.nonzero)
    P.norm = qnorm(P.cdf.emp)
    P.norm.t = t(P.norm)
    
    # standardized de-correlated residuals (in time)
    inno.std = t(P.norm.t[,2:n] - A%*%P.norm.t[,1:(n-1)]) %*% diag(1 / sdZ)
    
    # estimate number of degrees of freedom
    dfStudent = get.df.Student(inno.std,corZ)
    
    # return results
    return(list(M0=M0,A=A,covZ=covZ,sdZ=sdZ,corZ=corZ,dfStudent=dfStudent))
    
  }else{
    stop('fit.copula.amount: unknown value for copulaInt: must be Gaussian or Student')
  }
}


#==============================================================================
#' unif.to.prec
#'
#' from uniform variates to precipitation variates
#' @param pI vector of three parameters of the marginal distributions
#' @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#' @param U vector of uniform variates
#'
#' @return \item{matrix}{matrix of estimates p x 3}
#'
#' @author Guillaume Evin
unif.to.prec = function(pI,typeMargin,U){
  # inverse-cdf
  if(typeMargin == 'EGPD'){
    prec.sim = qEGPD.GI(U,pI[1],pI[2],pI[3])
  }else if(typeMargin == 'mixExp'){
    prec.sim = as.vector(Renext::qmixexp2(U,pI[1],pI[2],pI[3]))
  }
  
  return(prec.sim)
}


#==============================================================================
#' fit.GWex.prec
#'
#' estimate all the parameters for the G-Wex model of precipitation
#' 
#' @references Evin, G., A.-C. Favre, and B. Hingray. 2018. 'Stochastic Generation of Multi-Site 
#' Daily Precipitation Focusing on Extreme Events.' Hydrol. Earth Syst. Sci.
#' 22 (1): 655-672. doi.org/10.5194/hess-22-655-2018.
#'
#' @param objGwexObs object of class \code{\linkS4class{GwexObs}}
#' @param parMargin if not NULL, list where each element parMargin[[iM]] 
#' corresponds to a month iM=1...12 and contains a matrix nStation x 3 of estimated parameters
#'  of the marginal distributions (EGPD or mixture of exponentials)
#' @param listOption  list with the following fields:
#' \itemize{
#'   \item \strong{th}: threshold value in mm above which precipitation observations are considered to be non-zero (=0.2 by default)
#'   \item \strong{nLag}: order of he Markov chain for the transitions between dry and wet states (=2 by default)
#'   \item \strong{typeMargin}: 'EGPD' (Extended GPD) or 'mixExp' (Mixture of Exponentials). 'EGPD' by default
#'   \item \strong{copulaInt}: 'Gaussian' or 'Student': type of dependence for amounts (='Student' by default)
#'   \item \strong{isMAR}: logical value, do we apply a Autoregressive Multivariate Autoregressive model (order 1) =TRUE by default
#'   \item \strong{is3Damount}: logical value, do we apply the model on 3D-amount. =FALSE by default
#'   \item \strong{nChainFit}: integer, length of the runs used during the fitting procedure. =100000 by default
#'   \item \strong{nCluster}: integer, number of clusters which can be used for the parallel computation
#' }
#'
#' @export
#'
#' @return a list containing the list of options \code{listOption} and the list of estimated parameters \code{listPar}. 
#' The parameters of the occurrence process are contained in \code{parOcc} and the parameters related to the  precipitation 
#' amounts are contained in \code{parInt}. Each type of parameter is a list containing the estimates for each month. In \code{parOcc}, we find:
#'
#' \itemize{
#'   \item \strong{p01}: For each station, the probability of transition from a dry state to a wet state.
#'   \item \strong{p11}: For each station, the probability of staying in a wet state.
#'   \item \strong{list.pr.state}: For each station, the probabilities of transitions for a Markov chain with lag \code{p}.
#'   \item \strong{list.mat.omega}: The spatial correlation matrix of occurrences \eqn{\Omega} (see Evin et al., 2018).}
#'
#' In \code{parInt}, we have:
#'
#' \itemize{
#'  \item \strong{parMargin}: list of matrices nStation x nPar of parameters for the marginal distributions (one element per Class).
#'   \item \strong{cor.int}: Matrices nStation x nStation \eqn{M_0}, \eqn{A}, \eqn{\Omega_Z} representing 
#'   the spatial and temporal correlations between all the stations (see Evin et al., 2018). For the 
#'   Student copula, \code{dfStudent} indicates the \eqn{\nu} parameter.
#'   }
#' @author Guillaume Evin
fit.GWex.prec = function(objGwexObs,parMargin,listOption=NULL){
  
  # get/check options
  listOption = get.listOption(listOption)
  
  ######### Retrieve observations and dates ##########
  
  if(listOption$is3Damount){
    # Precipitation matrix (cumulated amounts on 3 days)
    P.1D = objGwexObs@obs
    P.mat = agg.matrix(P.1D,3)
    colnames(P.mat) = colnames(P.1D)
    # Dates
    n = nrow(P.1D)
    vec.dates = objGwexObs@date[seq(from=1,to=n-2,by=3)]
    # Time resolution
    dayScale = 3
  }else{
    # Precipitation matrix
    P.mat = objGwexObs@obs
    # Dates
    vec.dates = objGwexObs@date
    # Time resolution
    dayScale = 1
  }
  
  
  #########  Process dates and periods   ######### 
  
  # vector of integers indicating the months (1,1,1,...)
  vec.month = as.numeric(strftime(vec.dates, "%m"))
  
  # vector of 12 x 3-char "JAN" "FEB" ....
  vec.month.char = get.list.month()
  
  # number of stations
  p = ncol(P.mat)
  
  #########   check parMargin   #########  
  if(!is.null(parMargin)){
    text.err.parMargin = "parMargin must be a list where each element parMargin[[iM]]
  corresponds to a month iM=1...12 and contains a matrix nStation x 3 of estimated
  parameters"
    if(!is.list(parMargin)){
      stop(paste0("parMargin is not a list. ",text.err.parMargin))
    }else if(length(parMargin)!=12){
      stop(paste0("parMargin does not have 12 elements. ",text.err.parMargin))
    }else{
      for(iM in 1:12){
        if(any(dim(parMargin[[iM]])!=c(p,3))){
          stop(paste0("Error in parMargin for the month ",iM,". ",text.err.parMargin))
        }
      }
    }
  }
  
  ######### prepare results ##########
  list.pr.state = list.parMargin = list.mat.omega = list.par.dep.amount = list()
  
  # prepare parallelization
  if(listOption$isParallel){
    cl <- parallel::makeCluster(listOption$nCluster)
    doParallel::registerDoParallel(cl)
  }
  
  #==============================================================================
  # estimations by season for
  # - wet/dry transition probabilites
  # - the marginal distributions
  # - spatial process for wet.dry states
  # - spatial process for positive intensities
  #==============================================================================
  
  for(iMonth in 1:12){
    # prepare lists
    list.pr.state[[iMonth]] = list()
    list.parMargin[[iMonth]] = list()
    list.mat.omega[[iMonth]] = list()
    list.par.dep.amount[[iMonth]] = list()
    
    # month
    m = vec.month.char[iMonth]
    
    # period for this month
    period.m = get.period.fitting.month(m)
    
    # dates concerned by the 3-month window
    is3monthPeriod = vec.month%in%period.m
    isMonth = vec.month==iMonth
    
    
    #-------------------------------------------
    # wet/dry transition probabilities
    #-------------------------------------------
    pr.state = lag.trans.proba.mat.byClass(P.mat,isMonth,listOption$th,
                                           listOption$nLag,listOption$dayScale)
    list.pr.state[[iMonth]] = pr.state
    
    
    #-------------------------------------------
    # parameters of the marginal distributions
    #-------------------------------------------
    if(is.null(parMargin)){
      parMargin.class = fit.margin.cdf(P.mat,isMonth,listOption$th,listOption$typeMargin)
    }else{
      parMargin.class = parMargin[[iMonth]]
    }
    list.parMargin[[iMonth]] = parMargin.class
    
    
    #-------------------------------------------
    # - spatial process for wet.dry states
    #-------------------------------------------
    infer.mat.omega.out = infer.mat.omega(P.mat,is3monthPeriod,listOption$th,listOption$nLag,
                                          pr.state,listOption$nChainFit,
                                          listOption$isParallel)
    list.mat.omega[[iMonth]] = infer.mat.omega.out
    
    
    #-------------------------------------------
    # - spatial process for positive intensities
    #-------------------------------------------
    infer.dep.amount.out = infer.dep.amount(P.mat,is3monthPeriod,infer.mat.omega.out,
                                            listOption$nLag,listOption$th,parMargin.class,listOption$typeMargin,
                                            listOption$nChainFit,listOption$isMAR,listOption$copulaInt,
                                            listOption$isParallel)
    list.par.dep.amount[[iMonth]] = infer.dep.amount.out
  }
  
  if(listOption$isParallel){
    parallel::stopCluster(cl)
  }
  
  # return a list of all the parameters
  listPar = list(parOcc=list(list.pr.state=list.pr.state,list.mat.omega=list.mat.omega),
                 parInt=list(cor.int=list.par.dep.amount,parMargin=list.parMargin),
                 p=p)
  
  # return options and estimated parameters
  return(list(listOption=listOption,listPar=listPar))
}


###===============================###===============================###
#' disag.3D.to.1D
#' @useDynLib GWEX, .registration = TRUE
#' @noRd
#' @param Yobs matrix of observed intensities at 24h: (nTobs*3) x nStation
#' @param YObsAgg matrix of observed 3-day intensities: nTobs x nStation
#' @param mObsAgg vector of season corresponding to YobsAgg
#' @param YSimAgg matrix of simulated intensities per 3-day period: nTsim x nStation
#' @param mSimAgg vector of season corresponding to the period simulated
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#'
#' @author Guillaume Evin
disag.3D.to.1D = function(Yobs, # matrix of observed intensities at 24h: (nTobs*3) x nStation,
                          YObsAgg, # matrix of observed 3-day intensities: nTobs x nStation,
                          mObsAgg, # vector of season corresponding to YobsAgg
                          YSimAgg, # matrix of simulated intensities per 3-day period: nTsim x nStation
                          mSimAgg,  # vector of season corresponding to the period simulated
                          prob.class # vector of probabilities indicating class of "similar" mean intensities
                          
){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(YObsAgg))
  nTsim = as.integer(nrow(YSimAgg))
  
  ###### number of stations
  nStat = as.integer(ncol(Yobs))
  
  ##### class: classification of precipitation events according
  # to the mean precipitation over all stations. 4 classes
  classObs = vector(length = nTobs)
  classSim =  vector(length = nTsim)
  # for each season
  for(i.s in 1:4){
    # mean obs
    iObs.s = mObsAgg==i.s
    Yobs.s = YObsAgg[iObs.s,]
    mean.s = apply(Yobs.s,1,mean,na.rm=T)
    # 4 breaks by default: small, moderate, high, extremes precipitation
    q.mean.s = quantile(mean.s, probs=prob.class)
    if(any(q.mean.s==0)){
      q.mean.s = q.mean.s[q.mean.s!=0]
    }
    # observed class
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
    classObs[iObs.s] = class.s
    # simulated class
    iSim.s = mSimAgg==i.s
    Ysim.s = YSimAgg[iSim.s,]
    mean.s = apply(Ysim.s,1,mean,na.rm=T)
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
    classSim[iSim.s] = class.s
  }
  
  ###### replace NA values by -9999 (can be processed in Fortran)
  Yobs[is.na(Yobs)] = -9999
  YObsAgg[is.na(YObsAgg)] = -9999
  
  ##### call Fortran function
  disag3Day.out = .Fortran("disag3DayGWexPrec_F",  PACKAGE="GWEX",
                           Yobs=Yobs, Y3obs=YObsAgg, mObs=as.integer(mObsAgg), cObs=as.integer(classObs),
                           Y3sim=YSimAgg, mSim=as.integer(mSimAgg), cSim=as.integer(classSim),
                           nTobs=nTobs, nStat=nStat, nTsim=nTsim, nLagScore=as.integer(1),
                           Ysim=matrix(0,nrow=nTsim*3,ncol=nStat),codeDisag=matrix(0,nrow=nTsim,ncol=nStat))
  return(list(Ysim=disag3Day.out$Ysim,codeDisag=disag3Day.out$codeDisag))
}



#==============================================================================
#' sim.GWex.occ
#'
#' generate boolean variates which describe the dependence
#' between intersite occurrence correlations and wet/dry persistence
#' @param objGwexFit object of class GwexFit
#' @param vecMonth vector n x 1 of integers indicating the months
#'
#' @return \item{matrix of logical}{occurrences simulated}
#'
#' @author Guillaume Evin
sim.GWex.occ = function(objGwexFit,vecMonth){
  # number of stations
  p = getNbStations(objGwexFit)
  n = length(vecMonth)
  
  # number of days for the transition probas
  nLag = objGwexFit@fit$listOption$nLag
  
  # prepare simulation occurrence
  Xt = rndNorm = array(0,dim=c(n,p))
  
  # number of possible transitions
  n.comb = 2^nLag
  
  # parameters for the occurrence process
  parOcc = objGwexFit@fit$listPar$parOcc
  
  # a matrix of possible combination n.comb x nLag
  mat.comb = as.matrix(parOcc$list.pr.state[[1]][[1]][,1:nLag])
  
  # initialize array of transitions (Gaussian quantiles)
  Qtrans = array(0,dim=c(n,p,n.comb))
  
  for(t in 1:n){
    # filter
    iMonth.t = vecMonth[t]
    
    # prob / normal quantiles of transitions
    Ptrans.list = lapply(parOcc$list.pr.state[[iMonth.t]],'[',nLag+1)
    Qtrans.list = lapply(Ptrans.list,function(x) qnorm(unlist(x)))
    Qtrans.mat = matrix(unlist(Qtrans.list), ncol=n.comb, byrow=T)
    for(i.st in 1:p){
      for(i.comb in 1:n.comb){
        # fill array of transitions
        Qtrans[t,i.st,i.comb] = Qtrans.mat[i.st,i.comb]
      }
    }
    
    # generate multivariate gaussian
    rndNorm[t,] = MASS::mvrnorm(1,rep(0,p),parOcc$list.mat.omega[[iMonth.t]]$mat.omega)
  }
  
  
  # Start simulation
  for(t in (nLag+1):n){
    for(st in 1:p){
      comb = Xt[(t-nLag):(t-1),st]==1
      i.comb = which(apply(mat.comb, 1, function(x) all(x==comb)))
      qTr = Qtrans[t,st,i.comb]
      Xt[t,st] = (rndNorm[t,st]<=qTr)
    }
  }
  
  return(Xt)
}


#==============================================================================
#' sim.GWex.Yt.Pr.get.param
#'
#' get relevant parameters
#' @param objGwexFit object of class GwexFit
#' @param iM integer indicating the month
#'
#' @return \item{list}{list of parameters}
#'
#' @author Guillaume Evin
sim.GWex.Yt.Pr.get.param = function(objGwexFit,iM){
  # extract relevant parameters
  fitParCorIntclass = objGwexFit@fit$listPar$parInt$cor.int[[iM]]
  
  # return results
  return(fitParCorIntclass)
}


#==============================================================================
#' sim.Zt.Spatial
#'
#' generate gaussian variates which describe the spatial dependence between the sites
#' @param PAR parameters for a class
#' @param copulaInt 'Gaussian' or 'Student'
#' @param p number of stations
#'
#' @return \item{matrix}{matrix n x p of uniform dependent variates}
#'
#' @author Guillaume Evin
sim.Zt.Spatial = function(PAR,copulaInt,p){
  # generate from the corresponding multivariate distribution
  if(copulaInt=='Gaussian'){
    # generate from a multivariate Gaussian
    return(MASS::mvrnorm(n=1, mu=rep(0,p), Sigma=PAR[['M0']]))
  }else if(copulaInt=='Student'){
    # generate from a multivariate Student
    rmvt.out = mvtnorm::rmvt(n=1, df=PAR[['dfStudent']], sigma=PAR[['M0']])
    return(qnorm(matrix(pt(rmvt.out, df=PAR[['dfStudent']]), ncol = p)))
  }
}


#==============================================================================
#' sim.Zt.MAR
#'
#' generate gaussian variates which describe the spatial and temporal dependence
#' between the sites (MAR(1) process)
#' @param PAR parameters for this class
#' @param copulaInt 'Gaussian' or 'Student'
#' @param Zprev previous Gaussian variate
#' @param p number of stations
#'
#' @return \item{matrix}{matrix n x p of uniform dependent variates}
#'
#' @author Guillaume Evin
sim.Zt.MAR = function(PAR,copulaInt,Zprev,p){
  # generate from the corresponding multivariate distribution
  # t-1
  Zt.prev = t(PAR$A%*%Zprev)
  
  if(copulaInt=='Gaussian'){
    # generate from a multivariate Gaussian
    inno = MASS::mvrnorm(n=1, mu=rep(0,p), Sigma=PAR[['corZ']]) %*% diag(PAR[['sdZ']])
  }else if(copulaInt=='Student'){
    # generate from a Student copula with normal margins
    rmvt.out = mvtnorm::rmvt(n=1, df=PAR[['dfStudent']], sigma=PAR[['corZ']])
    r.cop.t.out = matrix(pt(rmvt.out, df=PAR[['dfStudent']]), ncol = p)
    inno = qnorm(r.cop.t.out)  %*% diag(PAR[['sdZ']])
  }
  
  # return MAR(1)
  return(Zt.prev + inno)
}


#==============================================================================
#' sim.GWex.Yt.Pr
#'
#' generate uniform variates which describe the dependence between intersite amount
#' correlations
#' @param objGwexFit object of class GwexFit
#' @param vecMonth vector n x 1 of integer indicating the months
#'
#' @return \item{matrix}{matrix n x p of uniform dependent variates}
#'
#' @author Guillaume Evin
sim.GWex.Yt.Pr = function(objGwexFit,vecMonth){
  # retrieve some options
  isMAR = objGwexFit@fit$listOption$isMAR
  copulaInt = objGwexFit@fit$listOption$copulaInt
  
  # number of stations
  p = getNbStations(objGwexFit)
  
  # length of the time series generated at the end
  n = length(vecMonth)
  
  # prepare matrix with Gaussian variates
  Yt.Gau = array(0,dim=c(n,p))
  
  #_____________ t=1 _________________
  # for the first iteration, we simulate from the marginal multivariate distribution (no temporal dependence)
  
  # retrieve period and parameters
  PAR = sim.GWex.Yt.Pr.get.param(objGwexFit,vecMonth[1])
  Yt.Gau[1,] = sim.Zt.Spatial(PAR,copulaInt,p)
  
  #_____________ t=2...n _____________
  for(t in 2:n){
    # retrieve period and parameters if necessary
    if((vecMonth[t-1]!=vecMonth[t])){
      PAR = sim.GWex.Yt.Pr.get.param(objGwexFit,vecMonth[t])
    }
    
    if(isMAR){
      Yt.Gau[t,] = sim.Zt.MAR(PAR,copulaInt,Yt.Gau[t-1,],p)
      
      # Issue with the MAR(1) process: instability of the autoregressive process if A is nearly non-inversible
      # during the simulation process, Z values can become very large (>10)and leads to prob=1, and non-possible invers transform
      # with the marginal distributions
      if(any(Yt.Gau[t,]>6)){
        print(paste0('fail: t=',t))
      }
    }else{
      Yt.Gau[t,] = sim.Zt.Spatial(PAR,copulaInt,p)
    }
  }
  
  # return variates
  return(pnorm(Yt.Gau))
}


#==============================================================================
#' sim.GWex.Yt
#'
#' Inverse PIT: from the probability space to the precipitation
#' space
#'
#' @param objGwexFit object of class GwexFit
#' @param vecMonth vector of integer indicating the months
#' @param Yt.Pr uniform variates describing dependence between inter-site amounts
#'
#' @return \item{matrix}{matrix n x p of simulated non-zero precipitation intensities}
#'
#' @author Guillaume Evin
sim.GWex.Yt = function(objGwexFit,vecMonth,Yt.Pr){
  # number of stations
  p = getNbStations(objGwexFit)
  
  # length of the simulated period: necessarily related to Yt.Pr
  n = nrow(Yt.Pr)
  
  # prepare simulation precip
  Yt = array(0,dim=c(n,p))
  
  # marginal distributions
  parMargin = objGwexFit@fit$listPar$parInt$parMargin
  
  # Start simulation
  for(iM in 1:12){
    # filter
    isClass = (vecMonth == iM)
    nClass = sum(isClass)
    
    # pour chaque station
    for(st in 1:p){
      # days for this period as matrix indices
      i.mat = cbind(which(isClass),rep(st,nClass))
      # intensities
      pI = parMargin[[iM]][st,]
      # inverse-cdf
      Yt[i.mat] = unif.to.prec(pI,objGwexFit@fit$listOption$typeMargin,Yt.Pr[i.mat])
    }
  }
  
  return(Yt)
}


#==============================================================================
#' mask.GWex.Yt
#'
#' Mask intensities where there is no occurrence
#'
#' @param Xt simulated occurrences
#' @param Yt simulated intensities
#'
#' @return \item{matrix}{matrix n x p of simulated precipitations}
#'
#' @author Guillaume Evin
mask.GWex.Yt = function(Xt,Yt){
  # number of stations
  p = ncol(Xt)
  
  # length of the time series generated
  n = nrow(Xt)
  
  # prepare simulation pluies
  Pt = array(0,dim=c(n,p))
  
  # pour chaque station
  for(st in 1:p){
    # wet days of this period as matrix indices
    is.wet = (Xt[,st]==1)
    n.wet = sum(is.wet)
    i.mat = cbind(which(is.wet),rep(st,n.wet))
    # mask intensities
    Pt[i.mat] = Yt[i.mat]
  }
  
  return(Pt)
}


#==============================================================================
#' sim.GWex.prec.1it
#'
#' Simulate one scenario of precipitation from the GWex model
#'
#' @param objGwexFit object of class GwexFit
#' @param vecDates vector of continuous dates
#' @param myseed seed of the random generation, to be fixed if the results need to be replicated
#' @param objGwexObs optional: necessary if we need observations to simulate (e.g. disaggregation of 3-day periods)
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#' @export
#' @return \item{matrix}{Precipitation simulated for the dates contained in vec.Dates at the different stations}
#'
#' @author Guillaume Evin
sim.GWex.prec.1it = function(objGwexFit,vecDates,myseed,objGwexObs,prob.class){
  # set seed of random generation
  set.seed(myseed)
  
  # do we simulate 3-day period
  is3Damount = objGwexFit@fit$listOption$is3Damount
  
  # tweak vector of dates
  if(is3Damount){
    n.orig = length(vecDates)
    # complete vector of dates if not a multiple of 3
    n = ceiling(n.orig/3)*3
    vecDates = c(vecDates,vecDates[(n.orig-2):n.orig])
    # get dates for these 3-day periods
    vecDates = vecDates[seq(from=1,to=n,by=3)]
    mSimAgg = month2season(as.numeric(format(vecDates,'%m')))
  }
  
  # vector of integers indicating the months (1,1,1,...)
  vecMonth = as.numeric(strftime(vecDates, "%m"))
  
  
  # Simulation of occurrences
  Xt = sim.GWex.occ(objGwexFit,vecMonth)
  
  # Simulation of the spatial and temporal dependence between amounts
  Yt.Pr = sim.GWex.Yt.Pr(objGwexFit,vecMonth)
  
  # We obtain directly related intensities
  Yt = sim.GWex.Yt(objGwexFit,vecMonth,Yt.Pr)
  
  # if we simulate 3-day periods, we disaggregate these periods, otherwise, we just mask the intensities
  # with the simulated occurrences
  if(is3Damount){
    # Time resolution
    dayScale = 3
    
    # mask simulated intensities with the occurrences Xt
    simAgg = mask.GWex.Yt(Xt,Yt)
    
    ### we need the observations to find observed 3-day structure of precipitation
    ### we first process observations: (aggregate; get months)
    #  get obs
    if(is.null(objGwexObs)|(!is.GwexObs(objGwexObs))) stop('sim.GWex.prec.1it: we need a proper objGwexObs object')
    obs = objGwexObs@obs
    # aggregate
    obsAgg = agg.matrix(obs,dayScale)
    nObsAgg = nrow(obsAgg)
    dateObsAgg = objGwexObs@date[seq(from=1,length.out=nObsAgg,by=dayScale)]
    mObsAgg = month2season(as.numeric(format(dateObsAgg,'%m')))
    #### disagregate 3-day int. Yt into daily intensities
    Pt = disag.3D.to.1D(obs[1:(nObsAgg*3),],obsAgg,mObsAgg,simAgg,mSimAgg,prob.class)$Ysim[1:n.orig,]
  }else{
    # mask with the occurrences Xt
    Pt = mask.GWex.Yt(Xt,Yt)
  }
  
  # return results
  return(Pt)
}


#==============================================================================
#' fit.margin.cdf
#'
#' estimate parameters which control the marginal distribution of precipitation amounts
#' @param P.mat precipitation matrix
#' @param isPeriod vector of logical n x 1 indicating the days concerned by a 3-month period
#' @param th threshold above which we consider that a day is wet (e.g. 0.2 mm)
#' @param type distribution: 'EGPD' or 'mixExp'
#'
#' @return \item{matrix}{matrix of estimates p x 3}
#' @export
#' @author Guillaume Evin
fit.margin.cdf = function(P.mat,isPeriod,th,type=c('EGPD','mixExp')){
  # number of stations
  p = ncol(P.mat)
  
  # class = Period
  isClass = isPeriod
  
  # filtered matrix of precipitation for this period (3-mon moving window by dafault)
  P.mat.class = P.mat[isClass,]
  
  # prepare output
  list.out = matrix(nrow=p,ncol=3)
  
  if(type == 'EGPD'){
    #  Applies Extented GPD
    for(i.st in 1:p){
      P.st = P.mat.class[,i.st]
      is.Prec = P.st>th&!is.na(P.st)
      P.nz = P.st[is.Prec]
      list.out[i.st,] = EGPD.GI.fit.PWM(P.nz)$x
    }
  }else if(type == 'mixExp'){
    #  Applies mixture of exponentials
    for(i.st in 1:p){
      P.st = P.mat.class[,i.st]
      is.Prec = P.st>th&!is.na(P.st)
      P.nz = P.st[is.Prec]
      list.out[i.st,] = Renext::EM.mixexp(P.nz)$estimate[c(1,3,4)]
    }
  }
  
  # return matrix of estimates
  return(list.out)
}