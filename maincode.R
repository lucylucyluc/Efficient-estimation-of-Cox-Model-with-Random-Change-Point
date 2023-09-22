### The Cox model with covariate-free random change point
RCPFREE = function(dat,par0,k,intebound,trunnorm){
  #dat: data
  #par0: the initial value of parameters (par0 = c(omega0,alpha0,beta0,mu0,sigma0,xi0))
  #k:number of inner knots
  #intebound: the low bound and up bound of integral
  #trunnorm : including the lower and upper bounds for truncated normal distribution
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  set.seed(0)
  n = dim(dat)[1]
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  z = dat[[3]]
  X = dat[[4]]   ##covariates 
  U = dat[[5]]   ## variable with threshold effect
  px = dim(X)[2]
  pz = dim(Z)[2]
  
  
  ######----- M,I spline 
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  
  dbasc=deriv(basc)     
  valdc=evaluate(dbasc,Y)
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  
  m=200
  emin = intebound[1]
  emax = intebound[2]
  tmin = trunnorm[1]
  tmax = trunnorm[2]
  grid=seq(emin,emax,length.out =m)
  pz = dim(Z)[2]
  px = dim(X)[2] 
  loglik=function(par){
    omega = par[1:pz]
    alpha  = par[(pz+1):(pz+px)]
    beta   = par[(pz+px+1):(pz+2*px)]
    muu    = par[pz+2*px+1]
    sigma  = par[pz+2*px+2]
    xi     = (par[(pz+2*px+3):(pz+2*px+2+od+k)])^2
    tt = matrix(rep(grid,n),nrow=m)
    t0 = matrix(rep(Z%*%omega,m),nrow=n)
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)
    t5 = t4>t(tt)
    t6 = exp(t0+t2+(t3*t5))
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)
    t9 = -(t(tt)-muu)^2/(2*sigma^2)
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = exp(t9)/(sqrt(2*pi)*sigma)#lam*exp(-lam*(t(tt)))
    term4 = pnorm((tmax-muu)/sigma) - pnorm((tmin-muu)/sigma)
    term5 = term3/term4
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term5*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term5*Gd,1,sum)
    lint  = log(int1+int2)
    loglf = -sum(lint)
    return(loglf)
  }
  par00   = c(par0,xi0) ### initial value
  result = tryCatch({ optim(par00,loglik,hessian = TRUE,control = list(maxit = 5000,abstol=10^-6,reltol=10^-6))},error = function(e){99} )
  
  if (is.double(result)==1){res1=NULL}else{
    est = result$par[1:(pz+2*px+2)]
    hes = result$hessian
    inhes = ginv(hes)
    sdinhes1 = sqrt(diag(inhes)[1:(pz+2*px+2)])
    pp = (1-pnorm(abs(est/sdinhes1),0,1))*2         
    res = rbind(est,sdinhes1,pp)
    BIC  = result$value + log(n)*(length(par0))
    res1 = list(res,result$convergence,BIC)}
  return(res1)
}


### The Cox model with covariate-dependent random change point
RCPDEPENDENT = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,Z,X,U)
  ## par0: initial value of parameters (par0 = c(omega0,alpha0,beta0,mu0,sigma0,xi0))
  ## k :  number of inner knots of spline
  ## intebound: including the Upper and lower bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  set.seed(0)
  n = dim(dat)[1]
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  z = dat[[3]]
  X = dat[[4]]   ##covariates 
  W = dat[[5]]
  U = dat[[6]]   ## variable with threshold effect
  px = dim(X)[2]
  pz = dim(Z)[2]
  pw = dim(W)[2]
  
  
  #####------ M,I spline 
  od <- 3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  
  dbasc=deriv(basc)     
  valdc=evaluate(dbasc,Y)
  
  
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  ############## M-spline
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  ## sieve loglikelihood function
  m=200
  inteup = intebound[1]
  intelow = intebound[2]
  grid=seq(intelow,inteup,length.out =m)
  seieveloglikehood1 = function(par){
    omega  = par[1:pz]
    alpha  = par[(pz+1):(pz+px)]
    beta   = par[(pz+px+1):(pz+2*px)]
    gamma  = par[(pz+2*px+1):(pz+2*px+pw)]
    sigma  = par[pz+2*px+pw+1]
    xi     = (par[(pz+2*px+pw+2):(pz+2*px+pw+p+1)])^2  
    tt = matrix(rep(grid,n),nrow=m)#m*n
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)
    t0 = matrix(rep(Z%*%omega,m),nrow=n)
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)
    t5 = t4>t(tt)
    t6 = exp(t0+t2+t3*t5)
    
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)
    t8 = matrix(rep(W%*%gamma,m),nrow=n)
    t9 = -(t(tt)-t8)^2/(2*sigma^2)
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = exp(t9)/(sqrt(2*pi)*sigma)
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd    = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term3*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term3*Gd,1,sum)
    lint  = int1+int2
    loglf = -sum(log(lint))
    return(loglf)
  }
  
  par00   = c(par0,xi0) 
  result = tryCatch({optim(par00,seieveloglikehood1,hessian = TRUE,control = list(maxit=5000,abstol=10^-6,reltol=10^-6))},error = function(e){99} )
  if (is.double(result)==1){res1=NULL}else{
    est = result$par[1:(pz+2*px+pw+1)]
    hes = result$hessian   
    inhes = ginv(hes)     
    sdinhes1 = sqrt(diag(inhes)[1:(pz+2*px+pw+1)]) 
    pp = (1-pnorm(abs(est/sdinhes1),0,1))*2         
    res = rbind(est,sdinhes1,pp)
    BIC  = result$value + log(n)*(length(par0))
    res1 = list(res,result$convergence,BIC)}
  return(res1)
}
