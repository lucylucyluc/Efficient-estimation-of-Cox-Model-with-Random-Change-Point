# Efficient-estimation-of-Cox-Model-with-Random-Change-Point

In the maincode, there are six codes for the paper "Efficient estimation of Cox Model with Random Change Point".

In the first code, "set1_exp",  we assume that the random change point is from a exponential distribution under the Cox Model with Covariate-free Random Change Point  without covariate Z. 
while the random change point is assumed to follow a normal distribution with mean mu and variance sigma^2 in code "set_norm". And the form of cox model is the same as set1_exp.

The code "set2" is for setting 2 in this paper, where the random change point also form a normal distribution. The difference between this code and set1_norm is that the form of the cox model contain covariate Z.

The code "set3" is for setting 3 in this paper, where the random change point has a linear structure with an error. And the error form a normal distribution. In the cox model, covariate Z = ∅.

The code "realdata1" is for the breast dataset under Cox Model with Covariate-free Random Change Point. And the error follows a truncated normal distribution on [25,65].
The code "realdata2" is for the breast dataset under Cox Model with Covariate-dependent Random Change Point. 

In all code, the "dat" is data, which is a list. The "par0" is the initial value that needs to be input, "k" is the number of inner knots of spline, and "intebound" include the lower and upper bounds of integral.

Example
Here, we will provide an example to illustrate how to run the code.

## set1_norm
###generating data process 
n=500
alpha0 = c(-0.3,1.2,0.6)
beta0 = c(1.4,-0.5,0.2)
mu0 = 0
lam0  = 1 
X1 = rnorm(n,1,2)     ###NORMIAL DISTRIBUTION
X2 = rbinom(n,1,0.6)    ### BINOMIAL DISTRIBUTION
X3 = runif(n,-2,2)
X = cbind(X1,X2,X3)   #### covariates

U = runif(n,-2,2)      
eta0 = rnorm(n,mu0,lam0)  
sigmaU = sd(U)

############# Cox model data
V = runif(n,0,1)
H = exp(X%*%alpha0 + X%*%beta0*(1*(U>eta0)))
TT = sqrt(-(log(1-V)/H))   ### lambda0=2*t  event time
C=runif(n,1,3)             ### noninformative censoring time 
delta = 1*(TT<=C)
Y = pmin(TT,C)             ### survival time
Y = c(Y)

xi0 = c()   # The intial value of nonparametric part
par0 = c(alpha0,beta0,mu0,lam0,xi0) #The initial value
