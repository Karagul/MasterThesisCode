# User-Defined Functions Called in the Rest of the Code Files -------------
# Last updated: Aug 23, 2017

# Before running the rest of the code files to replicate graphs in each chapter, 
# we need to first run this file containing user-defined functions to bring 
# all the user-defined functions into the local environment.


# to simulate square-root diffusion using Euler-Maruyama scheme -----------

EM_cir = function(n,kappa,sigma,theta,start,t0,t) {
  y=c()
  y[1] = start
  
  for (i in 1:n) {
    y[i+1] = y[i] + kappa*(t-t0)/n*(theta - y[i]) + sigma*(t-t0)/sqrt(n)*sqrt(y[i])*rnorm(1)
  }
  
  return(y)
}


# to simulate bond price paths by decreasing time to maturity ------------

C_tT.fun = function(t_i,a1,kappa,sigma,theta,t){
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  C_tT = -2*(exp(gamma*(t-t_i))-1)/
    (gamma-kappa+(gamma+kappa)*exp(gamma*(t-t_i)))
  
  return(C_tT)}

A_tT.fun = function(t_i,a1,kappa,sigma,theta,t){
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  A_tT = 2*kappa*theta/sigma^2 * 
    log(2*gamma*exp((gamma+kappa)*(t-t_i)/2)/
          (gamma-kappa+(gamma+kappa)*exp(gamma*(t-t_i))))
  
  return(A_tT)}


# to simulate bond price paths by increasing time to maturity -------------

C_tT.fun2 = function(t_i,a1,kappa,sigma,theta,gamma){
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  C_tT = -2*a1*(exp(gamma*t_i)-1)/
    (gamma-kappa+(gamma+kappa)*exp(gamma*t_i))
  
  return(C_tT)}

A_tT.fun2 = function(t_i,a1,kappa,sigma,theta,gamma){
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  A_tT = 2*kappa*theta/sigma^2 * 
    log(2*gamma*exp((gamma+kappa)*t_i/2)/
          (gamma-kappa+(gamma+kappa)*exp(gamma*t_i)))
  return(A_tT)}


# The bond pricing formula from p. 167 Lamberton-Lapeyre ------------------

phi.fun = function(q,a1,kappa,sigma,theta,gamma,t) {
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  phi = 2*kappa*theta/sigma^2 * 
    log(2*gamma*exp(t*(gamma+kappa)/2)/
          (gamma-kappa+(gamma+kappa)*exp(gamma*t)+
             sigma^2*q*(exp(gamma*t)-1)))
  
  return(phi)
}

psi.fun = function(q,a1,kappa,sigma,theta,gamma,t) {
  gamma = sqrt(kappa^2+2*sigma^2*a1)
  
  psi = - (q*(gamma+kappa+exp(gamma*t)*(gamma-kappa))+2*a1*(exp(gamma*t)-1))/
    (sigma^2*q*(exp(gamma*t)-1)+gamma-kappa+exp(gamma*t)*(gamma+kappa))
  
  return(psi)
}  

# reflected EM scheme -----------------------------------------------------

refEM_cir = function(n,kappa,sigma,theta,start,t0,t) {
  y=c()
  y[1] = start
  
  for (i in 1:n) {
    y[i+1] = abs(y[i] + kappa*(t-t0)/n*(theta - y[i]) + sigma*(t-t0)/sqrt(n)*sqrt(y[i])*rnorm(1))
  }
  
  return(y)
}


# direct sampling from non-central Chi-squared distribution ----------------

chi_sampling = function(n,kappa,sigma,theta,start,t0,t){
  cir = c()
  
  c = sigma^2*(1-exp(-kappa*(t-t0)/n))/(4*kappa)
  dof = 4*kappa*theta/sigma^2
  
  cir[1] = start
  
  for (i in 1:n) {
    lambda = 4*kappa*exp(-kappa*(t-t0)/n)/(sigma^2*(1-exp(-kappa*(t-t0)/n)))*cir[i]
    noncentralChi = rchisq(1, df=dof, ncp = lambda)
    cir[i+1] = c * noncentralChi
  }
  return(cir)
}

# end of functions code