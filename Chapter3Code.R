# Chapter 3 code ----------------------------------------------------------
# Last updated: Aug 23, 2017

# To simulate 1 default time which admits a stochastic hazard rate that 
# follows a square-root diffusion. We use the Euler-Maruyama discretization scheme.

set.seed(19930225)

numTime = 50 # number of time points
t0 = 0
t = 1
npath = 1
lambda0 = 0.01 # we start the cir process at 0.01
kappa1 = 0.1
sigma1 = 1
theta1 = 50

npath = 1 # number of sample paths to replicate
cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa1,
                        sigma = sigma1,
                        theta = theta1,
                        start = lambda0,
                        t0 = t0,
                        t = t))

hazardProcess = cumsum(cir_paths[2:(numTime+1),]*((t-t0)/numTime)+lambda0)
threshold = rexp(npath, rate = 1) #
defaultTime = min(which(hazardProcess > threshold))

# plot figure 3.1
yran <- range(hazardProcess) # plot range
plot(1:numTime, 
     hazardProcess, 
     type = "l", 
     ylim = yran,
     xlab = "Time", 
     ylab = expression(paste(Gamma[t], " follows a square-root diffusion")),
     main = expression(paste("1 Hazard Process Path; Default Occurs when ", Gamma[t]," Exceeds Threshold")))

abline(h=threshold, col="red")
text(t(c(0,threshold)), "threshold level", pos=4) 
abline(v=defaultTime, col="red")
text(t(c(defaultTime,0)), "default time", pos=3) 

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa, " = 0.1, ",
  bar(theta), " = 50, ",
  sigma, " = 1, npath = 1")),
  side = 3)
# end of plotting figure 3.1

# end of Chapter 3 code