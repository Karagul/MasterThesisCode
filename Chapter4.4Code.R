# Chapter 4.4 code --------------------------------------------------------
# Last updated: Aug 29, 2017

set.seed(19930225)

npath = 1000

# part I
# uncorrelated square-root diffusion

numTime = 50 # number of time points
t0 = 0
t = 1

kappa2 = 0.6 # parameter values from Duffie and Garleanu (2001)
sigma2 = 0.14
theta2 = 0.02
lambda0 = 0.01

cir_paths1 = replicate(
  npath, 
  EM_cir(
    n = numTime,
    kappa = kappa2,
    sigma = sigma2,
    theta = theta2,
    start = lambda0,
    t0 = t0,
    t = t))

icir_path1_Fig4.5 = exp(apply(
  -cir_paths1[1:numTime,]*((t-t0)/numTime), 
  2,
  cumsum)) 

# first, we vary the theta parameter in lambda_2t, holding kappa and sigma
# the same as ones in lambda_1t

DefaultProb_icir_Fig4.5 = c()

for (theta in seq(0.02,0.92,by=0.1)) { 
  cir_paths4 = replicate(
    npath, 
    EM_cir(
      n = numTime,
      kappa = kappa2,
      sigma = sigma2,
      theta = theta,
      start = lambda0,
      t = 1,
      t0 = 0))
  icir_path2_Fig4.5 = exp(apply(
    -cir_paths4[1:numTime,]*((t-t0)/numTime), 
    2, 
    cumsum)) 
  twoIcir = icir_path1_Fig4.5 * icir_path2_Fig4.5
  icir_path_avg_Fig4.5 = apply(twoIcir,1,mean,na.rm=T)
  
  DefaultProb_icir_Fig4.5 = cbind(
    DefaultProb_icir_Fig4.5, 
    c(0,1-icir_path_avg_Fig4.5))
}

# plot figure 4.5 
colfunc <- colorRampPalette(c("black", "red"))
yran <- range(DefaultProb_icir_Fig4.5,na.rm=T) # plot range

plot(0:numTime, 
     DefaultProb_icir_Fig4.5[,1], 
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time", ylab = expression("Default Probability"),
     main = "Probability of at least 1 Default; 2 Uncorrelated Square-Root Diffusion",
     col = colfunc(10)[1])
for(i in 2:10) # plot the other paths
  lines(0:numTime, DefaultProb_icir_Fig4.5[,i], col = colfunc(10)[i])

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa[1], " = ", kappa[2], " = 0.6, ", 
  bar(theta)[1], " = 0.02, ",
  bar(theta)[2], " = 0.02 ~ 0.9, ",
  sigma[1], " = ", sigma[2], " = 0.14, npath = 1000")),
  side = 3)

legend("topleft",
       inset=0.02,
       legend = c(expression(paste(bar(theta)[2], " = 0.02")), seq(0.12,0.92,by=0.1)),
       lty = 1,
       box.lty=0,
       col = colfunc(10))
# end of plotting figure 4.5

# second, we vary the kappa parameter in lambda_2t, holding theta and sigma
# the same as the theta and sigma in lambda_1t

DefaultProb_icir_Fig4.5_kappa = c()

for (kappa in seq(0.3,1.2,by=0.1)) {
  cir_paths = replicate(
    npath,
    EM_cir(
      n = numTime,
      kappa = kappa,
      sigma = sigma2,
      theta = theta2,
      start = lambda0,
      t = 1,
      t0 = 0))
  icir_path2_Fig4.5_kappa = exp(apply(
    -cir_paths[1:numTime,]*((t-t0)/numTime), 
    2,
    cumsum)) 
  twoIcir = icir_path1_Fig4.5 * icir_path2_Fig4.5_kappa
  icir_path_avg_Fig4.5_kappa = apply(twoIcir,1,mean,na.rm=T)
  
  DefaultProb_icir_Fig4.5_kappa = cbind(
    DefaultProb_icir_Fig4.5_kappa, 
    c(0,1-icir_path_avg_Fig4.5_kappa))
}

# plot figure 4.6
colfunc <- colorRampPalette(c("black", "red"))
yran <- range(DefaultProb_icir_Fig4.5_kappa,na.rm=T) # plot range

plot(0:numTime, 
     DefaultProb_icir_Fig4.5_kappa[,1], 
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time", 
     ylab = expression("Default Probability"),
     main = "Probability of at least 1 Default; 2 Uncorrelated Square-Root Diffusion",
     col = colfunc(10)[1])
for(i in 2:10) # plot the other paths
  lines(0:numTime, DefaultProb_icir_Fig4.5_kappa[,i], col = colfunc(10)[i])

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa[1], " = 0.6, ", 
  kappa[2],"= 0.3~1.2, ",
  bar(theta)[1], " = ", bar(theta)[2], " = 0.02, ",
  sigma[1], " = ", sigma[2]," = 0.14, npath = 1000")),
  side = 3)

legend("topleft",
       inset=0.02,
       legend = c(expression(paste(kappa[2], " = 0.3")), seq(0.4,1.2,by=0.1)),
       lty = 1,
       box.lty=0,
       col = colfunc(10))
# end of plotting figure 4.6

# third, we vary the sigma parameter in lambda_2t, holding kappa and theta
# the same as the kappa and theta in lambda_1t

DefaultProb_icir_Fig4.5_sigma = c()

for (sigma in seq(0.12,0.3,by=0.02)) {
  cir_paths4 = replicate(
    npath, 
    EM_cir(
      n = numTime,
      kappa = kappa2,
      sigma = sigma,
      theta = theta2,
      start = lambda0,
      t = 1,
      t0 = 0))
  
  icir_path2_Fig4.5_sigma = exp(apply(
    -cir_paths4[1:numTime,]*((t-t0)/numTime), 
    2,
    cumsum)) 
  twoIcir = icir_path1_Fig4.5 * icir_path2_Fig4.5_sigma
  icir_path_avg_Fig4.5_sigma = apply(twoIcir,1,mean,na.rm=T)
  
  DefaultProb_icir_Fig4.5_sigma = cbind(
    DefaultProb_icir_Fig4.5_sigma, 
    c(0,1-icir_path_avg_Fig4.5_sigma))
}

# plot figure 4.7
colfunc <- colorRampPalette(c("black", "red"))
yran <- range(DefaultProb_icir_Fig4.5_sigma,na.rm=T) # plot range

plot(0:numTime, 
     DefaultProb_icir_Fig4.5_sigma[,1], 
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time", 
     ylab = expression("Default Probability"),
     main = "Probability of at least 1 Default; 2 Uncorrelated Square-Root Diffusion",
     col = colfunc(10)[1])
for(i in 2:10) # plot the other paths
  lines(0:numTime, DefaultProb_icir_Fig4.5_sigma[,i], col = colfunc(10)[i])

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa[1]," = ",kappa[2]," = 0.6, ",
  bar(theta)[1]," = ", bar(theta[2]), " = 0.02, ",
  sigma[1], " = 0.14, ",
  sigma[2], " = 0.12~0.3, npath = 1000")),
  side = 3)

legend("topleft",
       inset=0.02,
       legend = c(expression(paste(sigma[2], " = 0.12")), seq(0.14,0.3,by=0.02)),
       lty = 1,
       box.lty=0,
       col = colfunc(10))
# end of plotting figure 4.7

# part II
# correlated square-root diffusion
# probability of at least 1 default
set.seed(19930225)
npath = 1000
a11 = 1
a12 = 0.5
a21 = 0.5
a22 = 1
theta3 = 0.8

cir_paths_Fig4.8 = replicate(npath, EM_cir(n = numTime,
                                    kappa = kappa2,
                                    sigma = sigma2,
                                    theta = theta2,
                                    start = lambda0,
                                    t = 1,
                                    t0 = 0))
cir_paths2_Fig4.8= replicate(npath, EM_cir(n = numTime,
                                     kappa = kappa2,
                                     sigma = sigma2,
                                     theta = theta3,
                                     start = lambda0,
                                     t = 1,
                                     t0 = 0))


Bond_tT_Euler_path1_Fig4.8 = exp(apply(
  -((a11+a21)*cir_paths_Fig4.8[1:numTime,])*((t-t0)/numTime), 
  2, 
  cumsum)) 
Bond_tT_Euler_path2_Fig4.8 = exp(apply(
  -((a12+a22)*cir_paths2_Fig4.8[1:numTime,])*((t-t0)/numTime), 
  2,
  cumsum))

TwoBonds = Bond_tT_Euler_path1_Fig4.8 * Bond_tT_Euler_path2_Fig4.8
Bond_tT_Euler_path_avg2_Fig4.8 = apply(TwoBonds,1,mean,na.rm=T)
DefaultProb_02_icir_Fig4.8 = c(0,1-Bond_tT_Euler_path_avg2_Fig4.8)

timestep = seq(from=t0,to=t,by=(t-t0)/numTime)

C_tT1 = C_tT.fun2(t_i = timestep,
                  a1 = a11+a21,
                  kappa = kappa2,
                  theta = theta2,
                  sigma = sigma2)
A_tT1 = A_tT.fun2(t_i = timestep,
                  a1 = a11+a21,
                  kappa = kappa2,
                  theta = theta2,
                  sigma = sigma2)
C_tT2 = C_tT.fun2(t_i = timestep,
                  a1 = a12+a22,
                  kappa = kappa2,
                  theta = theta3,
                  sigma = sigma2)
A_tT2 = A_tT.fun2(t_i = timestep,
                  a1 = a12+a22,
                  kappa = kappa2,
                  theta = theta3,
                  sigma = sigma2)


Prob_lambda1t <- exp(lambda0*C_tT1 + A_tT1)
Prob_lambda2t <- exp(lambda0*C_tT2 + A_tT2)
NeitherDefaultProb_01 <- Prob_lambda1t * Prob_lambda2t
OneDefault <- 1- NeitherDefaultProb_01

# plot figure 4.6 
yran <- range(OneDefault) # plot range

plot(0:numTime, 
     OneDefault ,
     type = "l", 
     ylim = yran, # plot the first path
     xlab = "Time", 
     ylab = expression("Default Probability"),
     main = "Probability of at least 1 Default by Time t; 2 Correlated Square-Root Diffusions")
lines(x=0:numTime, y=DefaultProb_02_icir_Fig4.8, col="green")

mtext(expression(paste(
  lambda[0], " = 0.01, ",
  kappa[1], " = ", kappa[2], " = 0.6, ", 
  bar(theta)[1], " = 0.02, ", 
  bar(theta)[2], " = 0.8, ",
  sigma[1], " = ", sigma[2], " = 0.14, npath = 1000")),
  side = 3)
legend("topleft",inset=0.01,legend=c("Closed-form solution", "Mean of exp(int CIR)"),
       col=c("black","green"), lty=1:1,cex=0.9,box.lty=0, y.intersp=0.8)
# end of plotting figure 4.6

# end of Chapter 4.4 code