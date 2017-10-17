# Chapter 2 code ----------------------------------------------------------
# Last updated: Aug 23, 2017

# To simulate 20 sample paths of homogenous Poisson processes 

# Main issue: on a given time interval we don't know in advance 
# the number of arrivals. So we simulate more than we need.
# We find this technique from the following link:
# https://stats.stackexchange.com/questions/148997/poisson-process-in-r-from-exponential-distribution

set.seed(19930225)

lambda = 0.2 # constant hazard rate
t = 100 # number of time points
npath = 20 # number of sample paths
n = qpois(1 - 1e-8, # Pr{N(t) <= n} <= 1 - 10^-8
          lambda = lambda * t) # maximum number of exponential rvs to generate

# simulate exponential interarrivals
X = matrix(
  rexp(n * npath, rate = lambda), # generate n*npath exponential rvs
  ncol = npath) # each column represents a sample path
S = apply(X, 2, cumsum) # compute interarrival times
S = rbind("T0" = rep(0, npath), S) # add a row for t0 = 0

# plot figure 2.1
matplot(x = S, 
        y = 0:n, 
        type="S", 
        col=c(1:10),
        xlim = c(0, 100),
        ylim= c(0,40),
        main = "Homogeneous Poisson Process Paths", 
        xlab = "time", 
        ylab = "N(t)")
mtext(expression(paste(gamma, " = 0.2, npath = 20")),
      side = 3)
# end of plotting figure 2.1

# end of Chapter 2 code
