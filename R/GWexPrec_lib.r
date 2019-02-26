###===============================###===============================###
### Guillaume Evin
### 18/03/2016, Grenoble
###  IRSTEA
### guillaume.evin@irstea.fr
###
###  Provide utilities for the estimation and simulation of the GWex
### multi-site precipitation model.
###
### Wilks, D.S. (1998) "Multisite generalization of a daily stochastic
### precipitation generation model", J Hydrol, 210: 178-191
###===============================###===============================###

#==============================================================================
# get.mat.omega
#
# find omega correlation leading to estimates cor between occurrences
#
# @param cor.obs matrix p x p of observed correlations between occurrences for all pairs of stations
# @param Qtrans.mat transition probabilities, 2 x ncomb matrix
# @param mat.comb matrix of logical: ncomb x nlag
# @param nLag order of the Markov chain
# @param nChainFit length of the simulated chains used during the fitting
#
# @return \item{matrix}{omega correlations for all pairs of stations}
#
# @author Guillaume Evin
get.mat.omega = function(cor.obs,Qtrans.mat,mat.comb,nLag,nChainFit){
  # number of stations
  p = ncol(cor.obs)

  # possible pairs
  vec.pairs = combn(1:p, 2)
  n.pairs = ncol(vec.pairs)

  # apply find.omega for each pair of stations
  omega.paral = foreach(i.pair=1:n.pairs, .combine='cbind') %dopar% {
    i = vec.pairs[1,i.pair]
    j = vec.pairs[2,i.pair]
    return(find.omega(cor.obs[i,j],Qtrans.mat[c(i,j),],mat.comb,nLag,nChainFit))
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
# find.omega
#
# finds the correlation between normal variates leading to correlation between occurrences
#
# @param rho.emp target correlation between occurences
# @param Qtrans.mat transition probabilities, 2 x ncomb matrix
# @param mat.comb matrix of logical: ncomb x nlag
# @param nLag order of the Markov chain
# @param nChainFit length of the simulated chains used during the fitting
#
# @return \item{scalar}{needed correlation}
#
# @author Guillaume Evin
find.omega = function(rho.emp,Qtrans.mat,mat.comb,nLag,nChainFit){
  f1 = cor.emp.occ(1,Qtrans.mat,mat.comb,nLag,nChainFit) - rho.emp
  if(f1<=0){
    return(1)
  }else{
    f = function(w){
      cor.emp.occ(w,Qtrans.mat,mat.comb,nLag,nChainFit) - rho.emp
    }
    return(uniroot(f,c(rho.emp,1),extendInt="upX", tol = 1e-3)$root)
  }
}



#==============================================================================
#' cor.emp.occ
#'
#' Finds observed correlations between occurrences corresponding
#' to a degree of correlation of Gaussian multivariate random numbers
#' @useDynLib GWEX
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
# get.mat.zeta
#
# find omega correlation leading to estimates cor between intensities
#
# @param cor.obs matrix p x p of observed correlations between intensities for all pairs of stations
# @param Qtrans.mat transition probabilities, 2 x ncomb matrix
# @param mat.comb matrix of logical: ncomb x nlag
# @param mat.omega omega correlations for all pairs of stations
# @param nLag order of the Markov chain
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
# @param nChainFit integer indicating the length of simulated chains
# 
# @return \item{matrix}{omega correlations for all pairs of stations}
#
# @author Guillaume Evin
get.mat.zeta = function(cor.obs,Qtrans.mat,mat.comb,mat.omega,nLag,par.margin,typeMargin,nChainFit){
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
  zeta.paral = foreach(i.pair=1:n.pairs, .combine='cbind') %dopar% {
    i = vec.pairs[1,i.pair]
    j = vec.pairs[2,i.pair]
    return(find.zeta(cor.obs[i,j],nChainFit,Xt[,c(i,j)],par.margin[c(i,j),],typeMargin))
  }

  # prepare results
  mat.zeta = matrix(1,nrow=p,ncol=p)
  diag(mat.zeta) = 1

  # fill matrix
  for(i.pair in 1:n.pairs){
    i = vec.pairs[1,i.pair]
    j = vec.pairs[2,i.pair]
    mat.zeta[i,j] = zeta.paral[i.pair]
    mat.zeta[j,i] = zeta.paral[i.pair]
  }

  return(mat.zeta)
}


#==============================================================================
# find.zeta
#
# finds the correlation between normal variates leading to correlation between intensities
#
# @param eta.emp target correlation between intensities
# @param nChainFit number of simulations
# @param Xt simulated occurrences, n x 2 matrix
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#
# @return \item{scalar}{needed correlation}
#
# @author Guillaume Evin
find.zeta = function(eta.emp,nChainFit,Xt,par.margin,typeMargin){
  f1 = cor.emp.int(1,nChainFit,Xt,par.margin,typeMargin) - eta.emp
  if(f1<=0){
    return(0.99999)
  }else{
    f = function(w){
      cor.emp.int(w,nChainFit,Xt,par.margin,typeMargin) - eta.emp
    }
    zeta = uniroot(f,c(eta.emp,1),extendInt="upX",tol = 1e-3)$root
    return(zeta)
  }
}


#==============================================================================
# cor.emp.int
#
# Finds observed correlations between intensities corresponding
# to a degree of correlation of Gaussian multivariate random numbers
# @param zeta correlation of Gaussian multivariates
# @param nChainFit number of simulated variates
# @param Xt simulated occurrences, n x 2 matrix
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
# @return \item{scalar}{correlation between simulated intensities}
#
# @author Guillaume Evin
cor.emp.int = function(zeta,nChainFit,Xt,par.margin,typeMargin){
  # generate the same random numbers if zeta=1, which avoids problems at the bounds (minimzation diverving, or
  # function having the same signs at both bounds if cor.emp.int(1,...) is slightly less than 0 if zeta=1)
  set.seed(1)
  # Simulation of the spatial and temporal dependence between amounts
  Yt.Pr = pnorm(MASS::mvrnorm(n=nChainFit, mu=rep(0,2), Sigma=matrix(c(1,zeta,zeta,1),2,2)))
  # We obtain directly related intensities
  Yt = array(0,dim=c(nChainFit,2))
  for(i.st in 1:2) Yt[,i.st] = unif.to.prec(par.margin[i.st,],Yt.Pr[,i.st],typeMargin)
  # Mask with occurrences
  Yt[Xt==0] = 0
  # cor.Pearson.cor
  cor.Pearson = cor(Yt[,1],Yt[,2])
  return(cor.Pearson)
}


#==============================================================================
# get.vec.autocor
#
# find rho autocorrelation leading to empirical estimates
#
# @param autocor.obs matrix p x p of observed correlations between intensities for all pairs of stations
# @param Qtrans.mat transition probabilities, 2 x ncomb matrix
# @param mat.comb matrix of logical: ncomb x nlag
# @param mat.omega omega correlations for all pairs of stations
# @param nLag order of the Markov chain
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
# @param nChainFit integer indicating the length of the simulated chains
#
# @return \item{matrix}{omega correlations for all pairs of stations}
#
# @author Guillaume Evin
get.vec.autocor = function(autocor.obs,Qtrans.mat,mat.comb,mat.omega,nLag,par.margin,typeMargin,nChainFit){
  # number of stations
  p = length(autocor.obs)

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

  # find autocor leading to obs autocor
  iSt = NULL
  vec.autocor = foreach(iSt=1:p, .combine='cbind') %dopar% {
     return(find.autocor(autocor.obs[iSt],nChainFit,Xt[,iSt],par.margin[iSt,],typeMargin))
  }

  return(as.vector(vec.autocor))
}


#==============================================================================
# find.autocor
#
# finds the autocorrelation leading to observed autocorrelation
#
# @param autocor.emp target correlation between intensities
# @param nChainFit number of simulations
# @param Xt simulated occurrences, nChainFit x 2 matrix
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#
# @return \item{scalar}{needed correlation}
#
# @author Guillaume Evin
find.autocor = function(autocor.emp,nChainFit,Xt,par.margin,typeMargin){
  f1 = autocor.emp.int(0.95,nChainFit,Xt,par.margin,typeMargin) - autocor.emp
  if(f1<=0){
    return(0.95)
  }else{
    f0 = autocor.emp.int(0,nChainFit,Xt,par.margin,typeMargin) - autocor.emp
    if(f0>=0){
      return(0)
    }else{
      f = function(w){
        autocor.emp.int(w,nChainFit,Xt,par.margin,typeMargin) - autocor.emp
      }
      return(uniroot(f,c(0,0.95),tol = 1e-3)$root)
    }
  }
}

#==============================================================================
# autocor.emp.int
#
# Finds empirical autocorrelations (lag-1) between intensities corresponding
# to a degree of autocorrelation of an AR(1) process
# @param rho autocorrelation of the AR(1) process
# @param nChainFit number of simulated variates
# @param Xt simulated occurrences, nChainFit x 2 matrix
# @param par.margin parameters of the margins 2 x 3
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
# @return \item{scalar}{correlation between simulated intensities}
#
# @author Guillaume Evin
autocor.emp.int = function(rho,nChainFit,Xt,par.margin,typeMargin){
  # control random seed
  set.seed(1)
  # Simulation from an AR(1) process
  Yt.AR1 = vector(length=nChainFit)
  Yt.AR1[1] = rnorm(1)
  for(t in 2:nChainFit) Yt.AR1[t] = rho*Yt.AR1[t-1] + rnorm(1,sd = sqrt(1-rho^2))
  # to proba
  Yt.Pr = pnorm(Yt.AR1)
  # Related intensities
  Yt = unif.to.prec(par.margin,Yt.Pr,typeMargin)
  # Mask with occurrences
  Yt[Xt==0] = 0
  # Resulting cor
  cor.Pearson = cor(Yt[1:(nChainFit-1)],Yt[2:nChainFit])
  return(cor.Pearson)
}

#==============================================================================
# QtransMat2Array
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
# joint.proba.occ
#
# joint probabilities of occurrences for all pairs of stations
#
# @param P matrix of precipitation
# @param th threshold
#
# @return \item{list}{list of joint probabilities}
#
# @author Guillaume Evin
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
# cor.obs.occ
#
# provide observed correlations between occurrences for all pairs of stations
# see Mhanna et al. (2012)
# @references Mhanna, Muamaraldin, and Willy Bauwens. “A Stochastic Space-Time Model for the Generation of Daily Rainfall in the Gaza Strip.” International Journal of Climatology 32, no. 7 (June 15, 2012): 1098–1112. doi:10.1002/joc.2305.
#
# @param pi00 joint probability of having dry states
# @param pi0 probability of having a dry state
# @param pi1 probability of having a wet state
#
# @return \item{scalar}{matrix of observed correlations}
#
# @author Guillaume Evin
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
# fit.copula.amount
#
# estimate parameters which control the spatial dependence between intensities using a copula
# @param P precipitation matrix
# @param copulaInt type of dependance between inter-site amounts: 'Gaussian' or 'Student'
# @param M0 Matrix containing the inter-site spatial correlations
#
# @return \item{list}{list of estimates (e.g., M0, dfStudent)}
#
# @author Guillaume Evin
fit.copula.amount = function(P,copulaInt,M0)
{
  if(copulaInt=='Gaussian'){
    # For 'Gaussian', we only compute the inter-site correlations
    return(list(M0=M0))

  }else if(copulaInt=='Student'){
    # For 'Student', we also compute the degree of freedom of the Student copula
    dfStudent = get.df.Student(P,M0)
    # return results
    return(list(M0=M0,dfStudent=dfStudent))

  }else{
    stop('fit.copula.amount: unknown value for copulaInt: must be Gaussian or Student')
  }
}


#==============================================================================
# fit.MAR1.amount
#
# estimate parameters which control the dependence between intensities with a
# MAR(1) process
# @param P precipitation matrix
# @param copulaInt type of dependance between inter-site amounts: 'Gaussian' or 'Student'
# @param M0 Matrix containing the inter-site spatial correlations
# @param A Matrix containing the autocorrelation (temporal) correlations
#
# @return \item{list}{list of estimates (e.g., M0, dfStudent)}
#
# @author Guillaume Evin
fit.MAR1.amount = function(P,copulaInt,M0,A){
  # dimensions
  p = ncol(P)
  n = nrow(P)

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
    # transform margins
    P.cdf.emp = get.emp.cdf.matrix(P)
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
# fit.margin.cdf
#
# estimate parameters which control the marginal distribution of precipitation amounts
# @param P precipitation matrix
# @param type distribution: 'EGPD' or 'mixExp'
# @param xiHat xi estimates if provided
#
# @return \item{matrix}{matrix of estimates p x 3}
#
# @author Guillaume Evin
fit.margin.cdf = function(P,type=c('EGPD','mixExp'),xiHat){
  # number of stations
  p = ncol(P)

  # prepare output
  list.out = matrix(nrow=p,ncol=3)

  if(type == 'EGPD'){
    #  Applies Extented GPD
    for(i.st in 1:p){
      P.st = P[,i.st]
      P.nz = P.st[!is.na(P.st)]
      xiHat.st = ifelse(is.null(xiHat),0.05,xiHat[i.st])
      list.out[i.st,] = NHRH.fit.PWM.cst(P.nz,xiHat.st)$x
    }
  }else if(type == 'mixExp'){
    #  Applies mixture of exponentials
    for(i.st in 1:p){
      P.st = P[,i.st]
      P.nz = P.st[!is.na(P.st)]
      list.out[i.st,] = Renext::EM.mixexp(P.nz)$estimate[c(1,3,4)]
    }
  }

  # return matrix of estimates
  return(list.out)
}


#==============================================================================
# unif.to.prec
#
# from uniform variates to precipitation variates
# @param pI vector of parameters
# @param U vector of uniform variates
# @param typeMargin type of marginal distribution: 'EGPD' or 'mixExp'
#
# @return \item{matrix}{matrix of estimates p x 3}
#
# @author Guillaume Evin
unif.to.prec = function(pI,U,typeMargin){
  # inverse-cdf
  if(typeMargin == 'EGPD'){
    prec.sim = qNHRH.GI(U,pI[1],pI[2],pI[3])
  }else if(typeMargin == 'mixExp'){
    prec.sim = as.vector(Renext::qmixexp2(U,pI[1],pI[2],pI[3]))
  }

  return(prec.sim)
}



#==============================================================================
# fit.GWex.prec
#
# estimate all the parameters for the G-Wex model of precipitation
# @references Evin, G., A.-C. Favre, and B. Hingray. 2018. “Stochastic Generation of Multi-Site Daily Precipitation Focusing on Extreme Events.” Hydrol. Earth Syst. Sci.
# 22 (1): 655–672. doi.org/10.5194/hess-22-655-2018.
#
# @param objGwexObs object of class \code{\linkS4class{GwexObs}}
# @param listOption  list with the following fields:
# \itemize{
#   \item \strong{th}: threshold value in mm above which precipitation observations are considered to be non-zero (=0.2 by default)
#   \item \strong{nLag}: order of he Markov chain for the transitions between dry and wet states (=2 by default)
#   \item \strong{typeMargin}: 'EGPD' (Extended GPD) or 'mixExp' (Mixture of Exponentials). 'EGPD' by default
#   \item \strong{xiHat}: pre-determined values for the xi parameters of the EGPD distribution on prec. amounts
#   \item \strong{copulaInt}: 'Gaussian' or 'Student': type of dependence for amounts (='Student' by default)
#   \item \strong{isMAR}: logical value, do we apply a Autoregressive Multivariate Autoregressive model (order 1) =TRUE by default
#   \item \strong{is3Damount}: logical value, do we apply the model on 3D-amount. =FALSE by default
#   \item \strong{nChainFit}: integer, length of the runs used during the fitting procedure. =100000 by default
#   \item \strong{nCluster}: integer, number of clusters which can be used for the parallel computation
# }
#
# @export
#
# @return a list containing the list of options \code{listOption} and the list of estimated parameters \code{listPar}. The parameters of the occurrence process are contained in \code{parOcc} and the parameters related to the  precipitation amounts are contained in \code{parInt}. Each type of parameter is a list containing the estimates for each month. In \code{parOcc}, we find:
#
# \itemize{
#   \item \strong{p01}: For each station, the probability of transition from a dry state to a wet state.
#   \item \strong{p11}: For each station, the probability of staying in a wet state.
#   \item \strong{list.pr.state}: For each station, the probabilities of transitions for a Markov chain with lag \code{p}.
#   \item \strong{list.mat.omega}: The spatial correlation matrix of occurrences \eqn{\Omega} (see Evin et al., 2018).}
#
# In \code{parInt}, we have:
#
# \itemize{
#   \item \strong{par.margin}: Matrices nStation x nPar of parameters for the marginal distributions.
#   \item \strong{cor.int}: Matrices nStation x nStation \eqn{M_0}, \eqn{A}, \eqn{\Omega_Z} representing the spatial and temporal correlations between all the stations (see Evin et al., 2018). For the Student copula, \code{dfStudent} indicates the \eqn{\nu} parameter.}
# @author Guillaume Evin
fit.GWex.prec = function(objGwexObs,listOption=NULL){
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
    if(!(typeMargin %in% c('mixExp','EGPD'))) stop('wrong value for typeMargin')
  }else{
    typeMargin = 'EGPD'
    listOption[['typeMargin']] = typeMargin
  }


  # xiHat
  if('xiHat' %in% names(listOption)){
    xiHat = listOption[['xiHat']]
    xiVec = as.vector(as.matrix(xiHat[,2:13]))
    if(!all(is.numeric(xiVec)&(xiVec>=0))) stop('wrong value for xiHat')
  }else{
    xiHat = NULL
    listOption[['xiHat']] = xiHat
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


  ######### Retrieve observations and dates ##########

  if(is3Damount){
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

  # Months
  vec.month = as.numeric(strftime(vec.dates, "%m"))

  # number of stations
  p = ncol(P.mat)

  # liste des mois
  vec.month.char = get.list.month()


  ######### prepare results ##########
  list.p01 = list.p11 = list.pr.state = list.mat.omega = list.par.dep.amount = list.par.margin = list()


  # number of months, should be 12...
  n.m = length(vec.month.char)

  # prepare parallelization
  cl <- parallel::makeCluster(nCluster)
  doParallel::registerDoParallel(cl)
  
 
  # for the progress bar
  pb <- txtProgressBar()

  for(m in vec.month.char){
    ########### fit on each month #########
    # mois concerne
    i.m = which(m==vec.month.char)
    is.month = vec.month==i.m
    # donnees filtrees
    P.mat.m = P.mat[is.month,]
    dates.m = vec.dates[is.month]

    ########### parameters for the occurrences #########
    # proba transition etats sec et humide
    p01 = dry.wet.trans.proba(P.mat.m,dates.m,th,dayScale)
    p11 = wet.wet.trans.proba(P.mat.m,dates.m,th,dayScale)
    list.p01[[m]] = p01
    list.p11[[m]] = p11

    # proba jointes etats sec et humide sur plusieurs jours
    pr.state = lag.trans.proba.proba(P.mat.m,dates.m,th,nLag,dayScale)
    list.pr.state[[m]] = pr.state

    # observed correlation of occurrences
    pi0 = dry.day.frequency(P.mat.m,th)
    pi1 = wet.day.frequency(P.mat.m,th)
    pi.occ = joint.proba.occ(P.mat.m,th=th)
    cor.occ.obs = cor.obs.occ(pi.occ$p00,pi0,pi1)

    # a matrix of possible combination n.comb x nLag
    mat.comb = as.matrix(pr.state[[1]][,1:nLag])

    # number of possible transitions
    n.comb = 2^nLag

    # prob / normal quantiles of transitions
    Ptrans.list = lapply(pr.state,'[',nLag+1)
    Qtrans.list = lapply(Ptrans.list,function(x) qnorm(unlist(x)))
    Qtrans.mat = matrix(unlist(Qtrans.list), ncol=n.comb, byrow=T)

    # filter infinite values if Pr = 0 or 1
    Qtrans.mat[Qtrans.mat==-Inf] = -10^5
    Qtrans.mat[Qtrans.mat==Inf] = 10^5

    # estimation of omega matrices
    mat.omega = get.mat.omega(cor.occ.obs,Qtrans.mat,mat.comb,nLag,nChainFit)
    mat.omega = modify.cor.matrix(mat.omega)
    list.mat.omega[[m]] = mat.omega

    ########### fit on 3-month periods (moving window) #########
    # period for this month
    per.m = get.period.fitting.month(m)

    # filtered data
    P.mat.per = P.mat[vec.month%in%per.m,]
    is.Zero = P.mat.per<=th&!is.na(P.mat.per)
    P.mat.per[is.Zero] = 0

    # number of obs. for this period
    n.per = nrow(P.mat.per)

    # filter non-zero amounts
    P.int = matrix(NA,nrow=n.per,ncol=p)
    is.Prec = P.mat.per>th&!is.na(P.mat.per)
    P.int[is.Prec] = P.mat.per[is.Prec]


    ########### marginal distributions of amounts #########
    if(is.null(xiHat)){
      xiHat.m = NULL
    }else{
      xiHat.m = xiHat[[m]]
    }
    par.margin = fit.margin.cdf(P.int,typeMargin,xiHat.m)
    list.par.margin[[m]] = par.margin

    ############## Inter-site dep. ############
    ### Saptial dep.: M0
    # observed correlations between intensities (zeo and positive)
    cor.int = cor(P.mat.per, use="pairwise.complete.obs")
    # find corresponding needed zeta correlations between simulated intensities
    mat.zeta = get.mat.zeta(cor.int,Qtrans.mat,mat.comb,mat.omega,nLag,par.margin,typeMargin,nChainFit)
    mat.zeta = modify.cor.matrix(mat.zeta)

    # find remaining parameters if necessary
    if(isMAR){
      ### Temporal dep.: M1
      vec.ar1.obs = vector(length=p)
      n = nrow(P.mat)
      for(i.st in 1:p){
        P.st = P.mat[,i.st]
        P.st[!(vec.month%in%per.m)] = NA
        vec.ar1.obs[i.st] = cor(P.st[1:(n-1)],P.st[2:n],use="pairwise.complete.obs")
      }

      # find corresponding ar(1) parameter
      vec.ar1 = get.vec.autocor(vec.ar1.obs,Qtrans.mat,mat.comb,mat.omega,nLag,par.margin,typeMargin,nChainFit)

      # apply a MAR(1) process: multivariate autoregressive process, order 1
      par.dep.amount = fit.MAR1.amount(P.int,copulaInt,M0=mat.zeta,A=diag(vec.ar1))
    }else{
      # spatial process
      par.dep.amount = fit.copula.amount(P.int,copulaInt,M0=mat.zeta)
    }

    # save list of parameters
    list.par.dep.amount[[m]] = par.dep.amount

    # progress bar
    setTxtProgressBar(pb, i.m/n.m)
  }
  close(pb)
  parallel::stopCluster(cl)

  # return a list of all the parameters
  listPar = list(parOcc=list(p01=list.p01,p11=list.p11,list.pr.state=list.pr.state,list.mat.omega=list.mat.omega),
                 parInt=list(cor.int=list.par.dep.amount,par.margin=list.par.margin),
                 p=p,periods=vec.month.char)

  # return options and estimated parameters
  return(list(listOption=listOption,listPar=listPar))
}

###===============================###===============================###
#' disag.3D.to.1D
#' @useDynLib GWEX
#' @noRd
#' @param Yobs matrix of observed intensities at 24h: (nTobs*3) x nStation
#' @param YObsAgg matrix of observed 3-day intensities: nTobs x nStation
#' @param mObsAgg vector of month corresponding to YobsAgg
#' @param YSimAgg matrix of simulated intensities per 3-day period: nTsim x nStation
#' @param mSimAgg vector of month corresponding to the period simulated
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#'
#' @author Guillaume Evin
#
# @author Guillaume Evin
disag.3D.to.1D = function(Yobs, # matrix of observed intensities at 24h: (nTobs*3) x nStation,
                          YObsAgg, # matrix of observed 3-day intensities: nTobs x nStation,
                          mObsAgg, # vector of month corresponding to YobsAgg
                          YSimAgg, # matrix of simulated intensities per 3-day period: nTsim x nStation
                          mSimAgg,  # vector of month corresponding to the period simulated
                          prob.class # vector of probabilities indicating class of "similar" mean intensities
                          
){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(YObsAgg))
  nTsim = as.integer(nrow(YSimAgg))
  
  ###### number of stations
  nStat = as.integer(ncol(Yobs))
  
  ##### class: classification of precipitation events according
  # to the mean precipitation over the Aare catchment. 4 classes
  classObs = vector(length = nTobs)
  classSim =  vector(length = nTsim)
  # for each season
  for(i.s in 1:4){
    # mean obs
    iObs.s = mObsAgg==i.s
    Yobs.s = YObsAgg[iObs.s,]
    mean.s = apply(Yobs.s,1,mean,na.rm=T)
    # 4 breaks: small, moderate, high, extremes precipitation
    q.mean.s = quantile(mean.s, probs=prob.class)
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
# sim.GWex.occ
#
# generate boolean variates which describe the dependence
# between intersite occurrence correlations and wet/dry persistence
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of objects of class 'Date'
#
# @return \item{matrix of logical}{occurrences simulated}
#
# @author Guillaume Evin
sim.GWex.occ = function(objGwexFit,vec.Dates){
  # number of stations
  p = getNbStations(objGwexFit)

  # number of days for the transition probas
  nLag = objGwexFit@fit$listOption$nLag

  # length of the time series generated
  n = length(vec.Dates)

  # prepare simulation occurrence
  Xt = rndNorm = array(0,dim=c(n,p))

  # number of possible transitions
  n.comb = 2^nLag

  # parameters for the occurrence process
  parOcc = objGwexFit@fit$listPar$parOcc

  # a matrix of possible combination n.comb x nLag
  mat.comb = as.matrix(parOcc$list.pr.state[[1]][[1]][,1:nLag])

  # initialise array of transitions (Gaussian quantiles)
  Qtrans = array(0,dim=c(n,p,n.comb))

  for(per in objGwexFit@fit$listPar$periods){
    # filter
    ff = get.vec.dates.filter(vec.Dates,per)
    n.per = sum(ff)

    # prob / normal quantiles of transitions
    Ptrans.list = lapply(parOcc$list.pr.state[[per]],'[',nLag+1)
    Qtrans.list = lapply(Ptrans.list,function(x) qnorm(unlist(x)))
    Qtrans.mat = matrix(unlist(Qtrans.list), ncol=n.comb, byrow=T)
    for(i.st in 1:p){
      for(i.comb in 1:n.comb){
        # find corresponding days for this combination of dry/wet states
        ind.mat = cbind(which(ff),rep(i.st,n.per),i.comb)
        # fill array of transitions
        Qtrans[ind.mat] = Qtrans.mat[i.st,i.comb]
      }
    }

    # generate multivariate gaussian
    rndNorm[ff,] = MASS::mvrnorm(n.per,rep(0,p),parOcc$list.mat.omega[[per]])
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
# sim.GWex.Yt.Pr.get.per
#
# get fitting periods
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of objects of class 'Date'
#
# @return \item{vector of char}{vector of char indicating periods ('ALL' or months ('JAN','FEB',...,'DEC'))}
#
# @author Guillaume Evin
sim.GWex.Yt.Pr.get.per = function(objGwexFit,vec.Dates){
  # get periods
  periods = objGwexFit@fit$listPar$periods

  #!!!! only works for a vector of months. If periods are 3-month, this needs to be changed
  if(length(periods)==1){
    vec.per = rep('ALL',length(vec.Dates))
  }else{
    vec.month = as.numeric(strftime(vec.Dates, "%m"))
    vec.per = get.list.month()[vec.month]
  }

  # return results
  return(vec.per)
}


#==============================================================================
# sim.GWex.Yt.Pr.get.update
#
# get days when parameters change
# @param vec.per vector of periods
#
# @return \item{vector of logical}{A vector logicals indicating if we must update the parameters}
#
# @author Guillaume Evin
sim.GWex.Yt.Pr.get.update = function(vec.per){
  # length of simulations
  n = length(vec.per)

  # difference of periods
  do.update = c(TRUE,vec.per[2:n]!=vec.per[1:(n-1)])

  # return results
  return(do.update)
}


#==============================================================================
# sim.GWex.Yt.Pr.get.param
#
# get relevant parameters
# @param objGwexFit object of class GwexFit
# @param per period (e.g. month)
#
# @return \item{list}{list of parameters}
#
# @author Guillaume Evin
sim.GWex.Yt.Pr.get.param = function(objGwexFit,per){
  # extract relevant parameters
  fitParCorIntPer = objGwexFit@fit$listPar$parInt$cor.int[[per]]

  # return results
  return(fitParCorIntPer)
}


#==============================================================================
# sim.Zt.Spatial
#
# generate gaussian variates which describe the spatial dependence between the sites
# @param PAR parameters for this period
# @param copulaInt 'Gaussian' or 'Student'
# @param p number of stations
#
# @return \item{matrix}{matrix n x p of uniform dependent variates}
#
# @author Guillaume Evin
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
# sim.Zt.MAR
#
# generate gaussian variates which describe the spatial and temporal dependence
# between the sites (MAR(1) process)
# @param PAR parameters for this period
# @param copulaInt 'Gaussian' or 'Student'
# @param Zprev previous Gaussian variate
# @param p number of stations
#
# @return \item{matrix}{matrix n x p of uniform dependent variates}
#
# @author Guillaume Evin
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
# sim.GWex.Yt.Pr
#
# generate uniform variates which describe the dependence between intersite amount
# correlations
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of objects of class 'Date'
#
# @return \item{matrix}{matrix n x p of uniform dependent variates}
#
# @author Guillaume Evin
sim.GWex.Yt.Pr = function(objGwexFit,vec.Dates){
  # retrieve some options
  isMAR = objGwexFit@fit$listOption$isMAR
  copulaInt = objGwexFit@fit$listOption$copulaInt

  # number of stations
  p = getNbStations(objGwexFit)

  # length of the time series generated at the end
  n = length(vec.Dates)

  # get periods
  vec.per = sim.GWex.Yt.Pr.get.per(objGwexFit,vec.Dates)

  # prepare matrix with Gaussian variates
  Yt.Gau = array(0,dim=c(n,p))

  # when do we update the parameters
  do.update = sim.GWex.Yt.Pr.get.update(vec.per)


  #_____________ t=1 _________________
  # for the first iteration, we simulate from the marginal multivariate distribution (no temporal dependence)

  # retrieve period and parameters
  per = vec.per[1]
  PAR = sim.GWex.Yt.Pr.get.param(objGwexFit,per)

  Yt.Gau[1,] = sim.Zt.Spatial(PAR,copulaInt,p)

  #_____________ t=2...n _____________
  for(t in 2:n){
    # retrieve period and parameters if necessary
    if(do.update[t]){
      PAR = sim.GWex.Yt.Pr.get.param(objGwexFit,vec.per[t])
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
# sim.GWex.Yt
#
# Inverse PIT: from the probability space to the precipitation
# space
#
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of objects of class 'Date'
# @param Yt.Pr niform variates describing dependence between inter-site amounts
#
# @return \item{matrix}{matrix n x p of simulated non-zero precipitation intensities}
#
# @author Guillaume Evin
sim.GWex.Yt = function(objGwexFit,vec.Dates,Yt.Pr){
  # number of stations
  p = getNbStations(objGwexFit)

  # length of the simulated period: necessarily related to Yt.Pr
  n = nrow(Yt.Pr)

  # prepare simulation pluies
  Yt = array(0,dim=c(n,p))

  # marginal distributions
  typeMargin = objGwexFit@fit$listOption$typeMargin

  # Start simulation
  for(per in objGwexFit@fit$listPar$periods){
    # filter
    is.per = get.vec.dates.filter(vec.Dates,per)
    n.per = sum(is.per)

    # pour chaque station
    for(st in 1:p){
      # days for this period as matrix indices
      i.mat = cbind(which(is.per),rep(st,n.per))
      # intensities
      pI = objGwexFit@fit$listPar$parInt$par.margin[[per]][st,]
      # inverse-cdf
      Yt[i.mat] = unif.to.prec(pI,Yt.Pr[i.mat],typeMargin)
    }
  }

  return(Yt)
}


#==============================================================================
# mask.GWex.Yt
#
# Mask intensities where there is no occurrence
#
# @param Xt simulated occurrences
# @param Yt simulated intensities
#
# @return \item{matrix}{matrix n x p of simulated precipitations}
#
# @author Guillaume Evin
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
# sim.GWex.prec.1it
#
# Simulate one scenario of precipitation from the GWex model
#
# @param objGwexFit object of class GwexFit
# @param vec.Dates vector of dates
# @param myseed seed of the random generation, to be fixed if the results need to be replicated
# @param objGwexObs optional: necessary if we need observations to simulate (e.g. disaggregation of 3-day periods)
# @param prob.class vector of probabilities indicating class of "similar" mean intensities
# @export
# @return \item{matrix}{Precipitation simulated for the dates contained in vec.Dates at the different stations}
#
# @author Guillaume Evin
sim.GWex.prec.1it = function(objGwexFit,vec.Dates,myseed,objGwexObs,prob.class){
  # set seed of random generation
  set.seed(myseed)

  # do we simulate 3-day period
  is3Damount = objGwexFit@fit$listOption$is3Damount

  # tweak vector of dates
  if(is3Damount){
    n.orig = length(vec.Dates)
    # complete vector of dates if not a multiple of 3
    n = ceiling(n.orig/3)*3
    vec.Dates = c(vec.Dates,vec.Dates[(n.orig-2):n.orig])
    # get dates for these 3-day periods
    vec.Dates = vec.Dates[seq(from=1,to=n,by=3)]
    mSimAgg = month2season(as.numeric(format(vec.Dates,'%m')))
  }

  # Simulation of occurrences
  Xt = sim.GWex.occ(objGwexFit,vec.Dates)

  # Simulation of the spatial and temporal dependence between amounts
  Yt.Pr = sim.GWex.Yt.Pr(objGwexFit,vec.Dates)

  # We obtain directly related intensities
  Yt = sim.GWex.Yt(objGwexFit,vec.Dates,Yt.Pr)

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
