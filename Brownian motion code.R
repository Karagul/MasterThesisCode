
# Last updated: Jul 3 2017 ------------------------------------------------

### 1 Simulating CIR SDE Sample Paths ###

x0 = runif(1)
my.fun = function() {
  n = 500
  
  k = runif(1)
  nu = runif(1)
  beta = runif(1)
  
  y[1] = x0
  
  for (i in 1:(n-1)) {
    y[i+1] = abs(y[i] + k/n*(beta - y[i]) + nu/sqrt(n)*sqrt(y[i])*rnorm(1))
  }
  return(y)
}


npath <- 400
W.t.paths <- replicate(npath, my.fun())

n = 1:500

### Plot of the paths
layout(cbind(1, 2), width = c(1, 1/2)) # plot layout
yran <- range(W.t.paths) # plot range
opar <- par(mar = c(4.5, 4.2, 1, 0.5)) # reduce space around plot 1, especially to the right
plot(n, W.t.paths[,1], type = "l", ylim = yran, # plot the first path
     xlab = "Time steps", ylab = expression("Yt that follows a CIR process"))
for(k in 2:npath) # plot the other paths
  lines(n, W.t.paths[,k], col = adjustcolor("black", alpha.f = 1/(1+npath/14)))
