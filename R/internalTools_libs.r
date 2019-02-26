###===============================###===============================###
### Guillaume Evin
### 22/01/2019, Grenoble
###  IRSTEA
### guillaume.evin@irstea.fr
###
### The following functions compute standard statistics for precipitation
### time series
###===============================###===============================###


#==============================================================================
wet.wet.trans.proba = function(mat.prec, # matrix of precipitation
                               vec.dates, # vector of dates
                               th, # threshold above which we consider that a day is wet
                               dayScale=1 # time resolution
){
  # wet.wet.trans.proba (computeStat_lib): wet.wet transition probability
  return(apply(mat.prec,2,wet.wet.trans.proba.vector,vec.dates,th,dayScale))
}


#==============================================================================
dry.wet.trans.proba = function(mat.prec, # matrix of precipitation
                               vec.dates, # vector of dates
                               th, # th above which we consider that a day is wet
                               dayScale=1 # time resolution
){
  # dry.wet.trans.proba (computeStat_lib): wet.wet transition probability
  return(apply(mat.prec,2,dry.wet.trans.proba.vector,vec.dates,th,dayScale))
}


#==============================================================================
# Dry day frequency DDF
#==============================================================================
dry.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # dry day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x<=th,na.rm=T),th))
}


#==============================================================================
# Wet day frequency WDF
#==============================================================================
wet.day.frequency = function(mat.prec, # matrix of precipitation
                             th # th above which we consider that a day is wet
){
  # wet day frequency (computeStat_lib)
  return(apply(mat.prec,2,function(x,th) mean(x>th,na.rm=T),th))
}


#==============================================================================
# transition probabilities on several consecutive days
#==============================================================================
lag.trans.proba.vector = function(vec.prec, # vector of precipitation
                                  vec.dates, #vector of dates
                                  th, # threshold
                                  nlag, # number of lag days
                                  dayScale # time resolution
){
  # lag.trans.proba.vector (computeStat_lib): number of consecutive days on
  # which we compute the probabilities
  
  ndays = nlag+1
  
  # filter consecutive dates
  has.all.days = diff(vec.dates,lag=nlag)==(nlag*dayScale)
  ind.has.all.days = which(has.all.days)
  
  # is.wet.lag gives the boolean values indicating if the time step 0, 1, ... ,0+nlag
  # are wet or not
  is.wet.lag = matrix(nrow=length(ind.has.all.days),ncol=ndays)
  for(i.lag in 0:nlag) is.wet.lag[,(i.lag+1)] = vec.prec[ind.has.all.days+i.lag]>th
  
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
  
  # return  wconditional probabilities
  return(comb.prT)
}


#==============================================================================
lag.trans.proba.proba = function(mat.prec, # matrix of precipitation
                                 vec.dates, # vector of dates
                                 th, # threshold above which we consider that a day is wet
                                 nlag, # number of lag days
                                 dayScale=1 # time resolution
){
  # lag.trans.proba.proba (computeStat_lib)
  return(apply(mat.prec,2,lag.trans.proba.vector,vec.dates,th,nlag,dayScale))
}


#==============================================================================
# wet-wet transition probability WWTP
#==============================================================================
wet.wet.trans.proba.vector = function(vec.prec, # vector of precipitation
                                      vec.dates, # vector of dates
                                      th, # threshold above which we consider that a day is wet
                                      dayScale # time resolution
){
  # wet.wet.trans.proba.vector (computeStat_lib): wet-wet transition probability WWTP
  
  # filter consecutive dates
  has.next.day = diff(vec.dates,lag=1)==dayScale
  
  # time t and lag time t+1
  t1 = which(has.next.day)
  t2 = t1+1
  
  # no nan
  nz = (!is.na(vec.prec[t1])) & (!is.na(vec.prec[t2]))
  t1nz = t1[nz]
  t2nz = t2[nz]
  
  # return  wet.wet transition probability
  return(mean(vec.prec[t1nz]>th&vec.prec[t2nz]>th)/mean(vec.prec[t1nz]>th))
}


#==============================================================================
# dry-wet transition probability DWTP
#==============================================================================
dry.wet.trans.proba.vector = function(vec.prec, # vector of precipitation
                                      vec.dates, # vector of dates
                                      th, # threshold above which we consider that a day is wet
                                      dayScale # time resolution
){
  # dry.wet.trans.proba.vector(computeStat_lib): dry-wet transition probability DWTP
  
  # filter consecutive dates
  has.next.day = diff(vec.dates,lag=1)==dayScale
  
  # time t and lag time t+1
  t1 = which(has.next.day)
  t2 = t1+1
  
  # no nan
  nz = (!is.na(vec.prec[t1])) & (!is.na(vec.prec[t2]))
  t1nz = t1[nz]
  t2nz = t2[nz]
  
  # return  wet.wet transition probability
  return(mean(vec.prec[t1nz]<=th&vec.prec[t2nz]>th)/mean(vec.prec[t1nz]<=th))
}


#========================================================================
# Modify a non-positive definite correlation matrix
# in order to have a positive definite matrix
#========================================================================
modify.cor.matrix = function(cor.matrix){
  # modify.cor.matrix (GWex_lib.r): Modify a non-positive definite correlation matrix
  # in order to have a positive definite matrix
  # Rousseeuw, P. J. and G. Molenberghs. 1993. Transformation of non positive semidefinite
  # correlation matrices. Communications in Statistics: Theory and Methods 22(4):965-984.
  # Rebonato, R., & Jackel, P. (2000). The most general methodology
  # to create a valid correlation matrix for risk management and
  # option pricing purposes. J. Risk, 2(2), 17-26.
  #
  # INPUT: cor.matrix: positive definite correlation matrix
  #
  # OUTPUT: modified correlation matrix
  
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


#========================================================================
# get.df.Student: this function estimates the
# nu (degrees of freedom) of the multivariate Student
# distribution
#========================================================================
get.df.Student = function(P,Sig,max.df=20){
  # get.df.Student (GWex_lib.r): estimate nu parameter (degrees of freedom) from a multivariate t distribution
  # when the correlation matrix Sig is given (McNeil et al. (2005) "Quantitative Risk Management")
  
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
# Fonction get.emp.cdf.matrix: Empirical cdf values for all columns of a matrix
#==============================================================================
get.emp.cdf.matrix = function(X,th=NULL){
  # get.emp.cdf.matrix (processData_lib): Empirical cdf values for all columns of a matrix
  #
  # X: matrix n x p
  # th (optional): threshold
  
  # number of columns
  p = ncol(X)
  
  # prepare output
  Y = matrix(NA,nrow=nrow(X),ncol=p)
  
  #
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
# Liste des mois
#==============================================================================
get.list.month = function(){
  return(c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AOU','SEP','OCT','NOV','DEC'))
}

#==============================================================================
# Liste des mois
#==============================================================================
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
# Liste des saisons
#==============================================================================
get.list.season = function(){
  return(c('DJF','MAM','JJA','SON'))
}

#==============================================================================
# transform vector of months to seasons
#==============================================================================
month2season = function(vecMonth){
  iSeason = c(1,1,2,2,2,3,3,3,4,4,4,1)
  return(iSeason[vecMonth])
}


#==============================================================================
# return a vector of boolean indicating sub-periods of vec.Dates
#==============================================================================
get.vec.dates.filter = function(vec.Dates,filter.name){
  # processData_lib
  #
  # INPUTS:
  # - vec.Dates: vector of Dates
  # - filter.name: indicate some months
  
  # vecteur des mois
  vec.month = as.numeric(strftime(vec.Dates, "%m"))
  
  if(filter.name=="ALL"){
    filt = rep(T,length(vec.month))
  }else if(filter.name%in%get.list.month()){
    # mois concerne
    i.m = which(filter.name==get.list.month())
    # filtre
    filt = vec.month==i.m
  }else if(filter.name%in%get.list.season()){
    if(filter.name=='DJF'){
      filt = vec.month%in%c(12,1,2)
    }else if(filter.name=='MAM'){
      filt = vec.month%in%c(3,4,5)
    }else if(filter.name=='JJA'){
      filt = vec.month%in%c(6,7,8)
    }else if(filter.name=='MAM'){
      filt = vec.month%in%c(9,10,11)
    }
  }else{
    stop("error in: get.vec.dates.filter, unknown filter.name value")
  }
  
  return(filt)
}




#==============================================================================
# Simple accumulation of a matrix
#==============================================================================
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