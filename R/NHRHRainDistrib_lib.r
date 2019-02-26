###===============================###===============================###
### Guillaume Evin
### 18/03/2016, Grenoble
###  LTHE
### guillaume_evin@yahoo.fr
###
###  Provide utilities for the estimation and simulation of the unified
### distributions for positive rainfall intensities
###
### Naveau, P., Huser, R., Ribereau, P. and Hannart, A. (2016) "Modeling
### jointly low, moderate and heavy rainfall intensities without a
### threshold selection", Water Resources Research, in press
###===============================###===============================###


###===============================###===============================###
# First parametric family for G(v): distribution,density and quantile function
###===============================###===============================###
NHRH.pGI = function(v,kappa) return(v^kappa)
NHRH.dGI = function(v,kappa) return(kappa*v^(kappa-1))
NHRH.qGI = function(p,kappa) return(p^(1/kappa))


###===============================###===============================###
# Distribution function of the unified distribution
###===============================###===============================###
pNHRH.GI = function(x,G.par,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  return(NHRH.pGI(pH,G.par))
}


###===============================###===============================###
# Quantile function of the unified distribution
###===============================###===============================###
qNHRH.GI = function(p,G.par,sig,xi){
  qG = NHRH.qGI(p,G.par)
  qH = Renext::qGPD(p=qG,scale=sig,shape=xi)
  return(qH)
}

###===============================###===============================###
# Random generation from the unified distribution
###===============================###===============================###
rNHRH.GI = function(n,G.par,sig,xi){
  u = runif(n=n)
  rNHRH = qNHRH.GI(u,G.par,sig,xi)
  return(rNHRH)
}


###===============================###===============================###
# Density function of the unified distribution
###===============================###===============================###
dNHRH.GI = function(x,G.par,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  dH = Renext::dGPD(x=x,scale=sig,shape=xi)
  dNHRH = dH*NHRH.dGI(pH,G.par)
  return(dNHRH)
}


###===============================###===============================###
# Probability Weighted Moments
###===============================###===============================###
NHRHI.mu0 = function(kappa,sig,xi){
  mu0 = (sig/xi)*(kappa*beta(kappa,1-xi)-1)
  return(mu0)
}
NHRHI.mu1 = function(kappa,sig,xi){
  mu1 = (sig/xi)*(kappa*(beta(kappa,1-xi)-beta(2*kappa,1-xi))-1/2)
  return(mu1)
}
NHRHI.mu2 = function(kappa,sig,xi){
  mu2 = (sig/xi)*(kappa*(beta(kappa,1-xi)-2*beta(2*kappa,1-xi)+beta(3*kappa,1-xi))-1/3)
  return(mu2)
}


###===============================###===============================###
# Parameter estimation with the PWM method
###===============================###===============================###
# set of equations which have to be equal to zero
NHRH.fPWM = function(par,pwm){
  kappa = par[1]
  sig = par[2]
  xi = par[3]

  y = numeric(3)
  y[1] = NHRHI.mu0(kappa,sig,xi) - pwm[1]
  y[2] = NHRHI.mu1(kappa,sig,xi) - pwm[2]
  y[3] = NHRHI.mu2(kappa,sig,xi) - pwm[3]

  return(y)
}

# numerical solver of the system of nonlinear equations
NHRH.fit.PWM = function(x){
  sam.pwm = c(EnvStats::pwMoment(x,k=0),
              EnvStats::pwMoment(x,k=1),
              EnvStats::pwMoment(x,k=2))

  return(nleqslv::nleqslv(c(2,sd(x),0.1),NHRH.fPWM,jac=NULL,pwm=sam.pwm))
}

###===============================###===============================###
# Parameter estimation with the PWM method: Xi=0
###===============================###===============================###
# set of equations which have to be equal to zero
NHRH.fPWM.cst = function(par,pwm,xi.val){
  kappa = par[1]
  sig = par[2]
  xi = xi.val

  y = numeric(2)
  y[1] = NHRHI.mu0(kappa,sig,xi) - pwm[1]
  y[2] = NHRHI.mu1(kappa,sig,xi) - pwm[2]

  return(y)
}

# numerical solver of the system of nonlinear equations
NHRH.fit.PWM.cst = function(x,xi.val){
  sam.pwm = c(EnvStats::pwMoment(x,k=0),
              EnvStats::pwMoment(x,k=1))
  if(xi.val<0.001) xi.val=0.001
  fit.out = nleqslv::nleqslv(c(2,sd(x)),NHRH.fPWM.cst,jac=NULL,pwm=sam.pwm,xi.val=xi.val)
  fit.out$x = c(fit.out$x,xi.val)
  return(fit.out)
}
