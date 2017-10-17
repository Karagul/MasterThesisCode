# Chapter 5 code ----------------------------------------------------------
# Last updated: Aug 29, 2017

set.seed(19930225)

numTime = 50
npath = 1000 # number of sample paths

# test performance of the EM scheme using 4 sets of parameters
lambda0_1 = 0.3
kappa1 = 0.1
sigma1 = 2
theta1 = 0.4
t0 = 0
t = 1

cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa1,
                        sigma = sigma1,
                        theta = theta1,
                        start = lambda0_1,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)
pct_neg = sum(colSums(is.na(cir_paths)) > 0)/npath
pct_neg # find the percentage of cir paths that at one point become negative

lambda0_2 = 0.1
kappa2 = 0.2
sigma2 = 1.2
theta2 = 0.2

cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa2,
                        sigma = sigma2,
                        theta = theta2,
                        start = lambda0_2,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)
pct_neg = sum(colSums(is.na(cir_paths)) > 0)/npath
pct_neg # find the percentage of cir paths that at one point become negative

lambda0_3 = 0.05
kappa3 = 0.4
sigma3 = 1
theta3 = 0.1

cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa3,
                        sigma = sigma3,
                        theta = theta3,
                        start = lambda0_3,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)
pct_neg = sum(colSums(is.na(cir_paths)) > 0)/npath
pct_neg # find the percentage of cir paths that at one point become negative

lambda0_4 = 0.01 
kappa4 = 0.6
sigma4 = 0.14
theta4 = 0.02

cir_paths = replicate(npath, 
                      EM_cir(
                        n = numTime,
                        kappa = kappa4,
                        sigma = sigma4,
                        theta = theta4,
                        start = lambda0_4,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)
pct_neg = sum(colSums(is.na(cir_paths)) > 0)/npath
pct_neg # find the percentage of cir paths that at one point become negative

# test performance of the reflected EM scheme

cir_paths = replicate(npath, 
                      refEM_cir(
                        n = numTime,
                        kappa = kappa1,
                        sigma = sigma1,
                        theta = theta1,
                        start = lambda0_1,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

cir_paths = replicate(npath, 
                      refEM_cir(
                        n = numTime,
                        kappa = kappa2,
                        sigma = sigma2,
                        theta = theta2,
                        start = lambda0_2,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)


cir_paths = replicate(npath, 
                      refEM_cir(
                        n = numTime,
                        kappa = kappa3,
                        sigma = sigma3,
                        theta = theta3,
                        start = lambda0_3,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

cir_paths = replicate(npath, 
                      refEM_cir(
                        n = numTime,
                        kappa = kappa4,
                        sigma = sigma4,
                        theta = theta4,
                        start = lambda0_4,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

# test performance of the step-by-step direct sampling

cir_paths = replicate(npath, 
                      chi_sampling(
                        n = numTime,
                        kappa = kappa1,
                        sigma = sigma1,
                        theta = theta1,
                        start = lambda0_1,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

cir_paths = replicate(npath, 
                      chi_sampling(
                        n = numTime,
                        kappa = kappa2,
                        sigma = sigma2,
                        theta = theta2,
                        start = lambda0_2,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)


cir_paths = replicate(npath, 
                      chi_sampling(
                        n = numTime,
                        kappa = kappa3,
                        sigma = sigma3,
                        theta = theta3,
                        start = lambda0_3,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

cir_paths = replicate(npath, 
                      chi_sampling(
                        n = numTime,
                        kappa = kappa4,
                        sigma = sigma4,
                        theta = theta4,
                        start = lambda0_4,
                        t0 = t0,
                        t = t))
mean(cir_paths[(numTime+1),1:npath],na.rm=T)

# end of chapter 5 code