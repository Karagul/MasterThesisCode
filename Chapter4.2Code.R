# Chapter 4.2 code --------------------------------------------------------
# Last updated: Aug 24, 2017

set.seed(19930225)

npath = 100

# simulate bond price sample paths using the closed-form solution

numTime = 50
t0 = 0
t = 1
kappa2 = 0.6 # parameter values from Duffie and Garleanu (2001)
sigma2 = 0.14
theta2 = 0.02
lambda0 = 0.01
a1 = 1

cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa2,
                        sigma = sigma2,
                        theta = theta2,
                        start = lambda0,
                        t0 = t0,
                        t = t))  

timestep = seq(from=t0,to=t,by=(t-t0)/numTime)

C_tT = C_tT.fun(
  t_i = timestep,
  a1 = a1,
  kappa = kappa2,
  sigma = sigma2,
  theta = theta2,
  t = t)

A_tT = A_tT.fun(
  t_i = timestep,
  a1 = a1,
  kappa = kappa2,
  sigma = sigma2,
  theta = theta2,
  t = t)

Bond_01_closedform <- exp(cir_paths*C_tT + A_tT)

# plot figure 4.2
yran <- range(Bond_01_closedform[,1:100],na.rm=T) 

plot(x=0:numTime, 
     y=Bond_01_closedform[,1], # plot the first sample path
     type = "l",
     ylim = yran,
     xlab = "Time", ylab = "Bond price p(t,1); t ranges from 0 to 1 in 50 time steps",
     main = expression(paste("Bond Price Sample Paths; ", 
                             lambda[t],
                             " Follows a Square-Root Diffusion")))
for(k in 2:npath) # plot the other sample paths
  lines(0:numTime, Bond_01_closedform[,k], col = adjustcolor("black", alpha.f=1/(1+100/14)))

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa, " = 0.6, ",
  bar(theta), " = 0.02, ", 
  sigma, " = 0.14, npath = 100")),
  side = 3)
# end of plotting figure 4.2

# end of Chapter 4.2 code