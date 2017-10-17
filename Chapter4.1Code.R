# Chapter 4.1 code --------------------------------------------------------
# Last Updated: Aug 24, 2017

set.seed(19930225)

numTime = 50 # number of time points
t0 = 0
t = 1

kappa2 = 0.6 # parameter values from Duffie and Garleanu (2001)
sigma2 = 0.14
theta2 = 0.02
lambda0 = 0.01

npath = 100 # number of sample paths

# Running the code below will produce 6 null values (NaN) because 
# EM scheme cannot guarantee that the CIR values would not reach 
# go negative. However, since only 6% (6 out of 100) of the total 
# number of sample paths eventually reach 0, and the EM scheme 
# produces the closest result to the theoretical mean compared to 
# other discretization schemes, we still use the EM scheme to
# simulate the square-root diffusion.
cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa2,
                        sigma = sigma2,
                        theta = theta2,
                        start = lambda0,
                        t0 = t0,
                        t = t))  

pct_neg = sum(colSums(is.na(cir_paths)) > 0)/npath
pct_neg # find the percentage of cir paths that at one point become negative

# plot Figure 4.1
yran <- c(0,max(cir_paths,na.rm=T)) # plot range

plot(0:numTime, # plot the first sample path
     cir_paths[,1],
     type = "l", 
     ylim = yran, 
     xlab = "Time", 
     ylab = expression(paste(lambda[t], " follows a square-root diffusion")),
     main = "Square-Root Diffusion Sample Paths (EM Scheme)",
     yaxt="n") # remove y axis
for(k in 2:npath) # plot the other sample paths
  lines(0:numTime, cir_paths[,k], col = adjustcolor("black", alpha.f = 1/(1+100/14)))

axis(2, at=seq(0,0.09,by=0.01),labels=seq(0,0.09,by=0.01), las=2)

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa, " = 0.6, ", 
  bar(theta), " = 0.02, ", 
  sigma, " = 0.14, npath = 100")),
  side = 3)
# end of plotting Figure 4.1

# end of Chapter 4.1 code 

