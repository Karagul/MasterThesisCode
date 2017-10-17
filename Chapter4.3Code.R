# Chapter 4.3 code --------------------------------------------------------
# Last updated: Aug 24, 2017

set.seed(19930225)

numTime = 50 # number of time points
t0 = 0
t = 1

kappa2 = 0.6 # parameter values from Duffie and Garleanu (2001)
sigma2 = 0.14
theta2 = 0.02
lambda0 = 0.01

npath = 1000 # number of sample paths
cir_paths = replicate(npath, 
                      EM_cir(
                       n = numTime,
                       kappa = kappa2,
                       sigma = sigma2,
                       theta = theta2,
                       start = lambda0,
                       t0 = t0,
                       t = t))

t2 = 2

theta3 = 0.8 # Increase theta to 0.8

cir_paths_paramchange = t(matrix(unlist(lapply(
  cir_paths[dim(cir_paths)[1],], EM_cir,
  n = numTime+1,
  kappa = kappa2,
  sigma = sigma2,
  theta = theta3,
  t0 = t,
  t = t2)),
  ncol=numTime+2, 
  byrow=T))

cir_paths2 = rbind(cir_paths, cir_paths_paramchange[2:(numTime+1),])

t3 = 3

cir_paths_paramchange2 = t(matrix(unlist(lapply(
  cir_paths2[dim(cir_paths2)[1],], EM_cir,
  n = numTime+1,
  kappa = kappa2,
  sigma = sigma2,
  theta = theta2, # Change theta back to 0.02
  t0 = t2,
  t = t3)),
  ncol=numTime+2, 
  byrow=T))

cir_paths3 = rbind(cir_paths2, cir_paths_paramchange2[2:(numTime+1),])

# plot figure 4.3
yran <- range(cir_paths3,na.rm=T) # plot range

plot(0:(numTime*3), 
     cir_paths3[,1], 
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time",
     ylab = expression(paste(lambda[t], " is a square-root diffusion")),
     main = expression(paste(
       "Square-Root Diffusion; Parameter ", 
       bar(theta), 
       " Changes Value at Time 50 and 100")))

for(k in 2:100) # plot the other paths
  lines(0:(numTime*3), cir_paths3[,k], col = adjustcolor("black", alpha.f = 1/(1+100/14)))

mtext(expression(paste(lambda[0], " = 0.01, ",
                       kappa, " = 0.6, ", 
                       bar(theta)[0-50], " = 0.02, ", 
                       bar(theta)[50-100], " = 0.8, ",
                       bar(theta)[100-150], " = 0.02, ",
                       sigma, " = 0.14, npath = 100")),
      side = 3)

abline(v=50, col="red")
abline(v=100, col="red")
# end of plotting figure 4.3

# use calculate the default probability using closed-form solution
# first, between time 0 to 50 

a1 = 1
timestep = seq(from=t0,to=t,by=(t-t0)/numTime)

C_tT1 = C_tT.fun2(t_i = timestep,
                  a1 = a1,
                  kappa = kappa2,
                  theta = theta2,
                  sigma = sigma2)
A_tT1 = A_tT.fun2(t_i = timestep,
                  a1 = a1,
                  kappa = kappa2,
                  theta = theta2,
                  sigma = sigma2)

DefaultProb_01 <- 1-exp(lambda0*C_tT1 + A_tT1)

# second, between time 51 to 100

C_tT2 = C_tT.fun2(t_i = timestep,
                  a1 = a1,
                  kappa = kappa2,
                  theta = theta3,
                  sigma = sigma2)
A_tT2 = A_tT.fun2(t_i = timestep,
                  a1 = a1,
                  kappa = kappa2,
                  theta = theta3,
                  sigma = sigma2)

phi = phi.fun(q = -C_tT2,
              a1 = a1,
              t = 1,
              kappa = kappa2,
              theta = theta2,
              sigma = sigma2)
psi = psi.fun(q = -C_tT2,
              a1 = a1,
              t = 1,
              kappa = kappa2,
              theta = theta2,
              sigma = sigma2)

DefaultProb_12 <- 1-exp(A_tT2)*exp(lambda0*psi + phi)
DefaultProb_02 = c(DefaultProb_01[1:numTime],DefaultProb_12)

# use calculate the default probability by averaging simulations

Bond_tT_Euler_path2 = exp(apply(
  -cir_paths2[1:(numTime*2),]*((t-t0)/numTime),2,cumsum))
Bond_tT_Euler_path_avg2 = apply(Bond_tT_Euler_path2,1,mean,na.rm=T)
DefaultProb_02_icir = 1-Bond_tT_Euler_path_avg2
DefaultProb_02_icir = c(0,DefaultProb_02_icir)

# plot figure 4.4 
yran <- range(DefaultProb_02) # plot range
plot(0:(numTime*2), 
     DefaultProb_02, 
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time", 
     ylab = expression("Probability of Default"),
     main = expression(paste(
       "Probability of Default; Square-Root Diffusion Parameter ", 
       bar(theta), 
       " Changes Value at Time 50")))
lines(x=0:(numTime*2), y=DefaultProb_02_icir[1:101],col="green")

abline(v=50, col="red")

mtext(expression(paste(lambda[0], " = 0.01, ",
                       kappa, " = 0.6, ", 
                       bar(theta), " = 0.02, ", 
                       bar(theta)["new"], " = 0.8, ",
                       sigma, " = 0.14, npath = 1000")),
      side = 3)

legend("topleft",inset=0.01,legend=c("Closed-form solution", "Mean of simulations"),
       col=c("black","green"), lty=1:1,cex=0.9,box.lty=0, y.intersp=0.8)
# end of plotting figure 4.4
