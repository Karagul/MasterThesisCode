
# Last updated: Jul 5 2017 ------------------------------------------------

### 1 Simulating CIR SDE Sample Paths ###

x0 = runif(1) # starting point of y_t
t = 1000 # number of time points

my.fun = function() {
  n = t
  
  k = runif(1)
  nu = runif(1)
  beta = runif(1)
  
  y[1] = x0
  
  for (i in 1:(n-1)) {
    y[i+1] = abs(y[i] + k/n*(beta - y[i]) + nu/sqrt(n)*sqrt(y[i])*rnorm(1))
  }
  return(y)
}


npath = 400 # number of sample paths
W.t.paths = replicate(npath, my.fun())

n = 1:t

### Plot of the paths ###
yran <- range(W.t.paths) # plot range
plot(n, W.t.paths[,1], type = "l", ylim = yran, # plot the first path
     xlab = "time t", ylab = expression("Y(t)"),
     main = "CIR Process paths")
for(k in 2:npath) # plot the other paths
  lines(n, W.t.paths[,k], col = adjustcolor("black", alpha.f = 1/(1+npath/14)))

### 2 Simulation of Poisson processes ###

lambda <- 0.05
tMax <- 100

## find the number 'n' of exponential r.vs required by imposing that
## Pr{N(t) <= n} <= 1 - eps for a small 'eps'
n <- qpois(1 - 1e-8, lambda = lambda * tMax)

## simulate exponential interarrivals
X <- rexp(n = n, rate = lambda)
S <- c(0, cumsum(X))
plot(x = S, y = 0:n, type = "s", xlim = c(0, tMax)) 

## several paths?
npath <- 50
## simulate exponential interarrivals
X <- matrix(rexp(n * npath, rate = lambda), ncol = npath,
            dimnames = list(paste("S", 1:n, sep = ""), paste("samp", 1:npath)))
## compute arrivals, and add a fictive arrival 'T0' for t = 0
S <- apply(X, 2, cumsum)
S <- rbind("T0" = rep(0, npath), S)
head(S)
## plot using steps
matplot(x = S, y = 0:n, type = "s", col = "darkgray",
        xlim = c(0, tMax),
        main = "Homogeneous Poisson Process paths", xlab = "time t", 
        ylab = "N(t)")