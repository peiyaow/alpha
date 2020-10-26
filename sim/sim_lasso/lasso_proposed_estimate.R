# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------

library(glmnet)
require(methods)
library(rlist)

# path = "~/Documents/GitHub/alpha/"
path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "/function/main_function.R"))
source(paste0(path, "/function/sim_function.R"))
source(paste0(path, "/function/glmnet_POET.R"))

# set.seed(706)
p = 200
n = 100
K = 3

p_share = 10 # ridge 80 lasso 5
n_label = 3
si = 2

# control parameter for beta
h = 1 # share

ds = 10 # sqrt(p) multiplier parameter for spike L 
rho = 1/40 # determine spike smaller than 1/K
rrho = 0.1 # correlation for the core sigma matrix
d = .05 # uniform range

spike_factor1 = c(1, sapply(2:K, function(k) 1/gamma(k + 1)))*.7
spike_factor2 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.25)))
spike_factor3 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.5)))*1.3

spike1 = spike_factor1 * ds
spike2 = spike_factor2 * ds
spike3 = spike_factor3 * ds

du = (spike1[K]*rho+spike2[K]*rho + spike3[K]*rho)/3

para1 = FactorModelPara(n, p, K, spike = spike1, d, du = du, rrho)
para2 = FactorModelPara(n, p, K, spike = spike2, d, du = du, rrho)
para3 = FactorModelPara(n, p, K, spike = spike3, d, du = du, rrho)

u1 = 1
u2 = 2
u3 = 3

L1 = para1$L
L2 = para2$L
L3 = para3$L

sigma.vec = rep(si, n_label)

n.train.vec = c(n, n, n)
n.test.vec = c(n, n, n)
ix.vec = c(0, cumsum(n.train.vec))
label.train = as.factor(c(rep(1, n.train.vec[1]), rep(2, n.train.vec[2]), rep(3, n.train.vec[3])))
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

for (s in seq(0, 5, by = 0.1)[1:2]){
  gamma1 = c(1, 1, 2)*s
  gamma2 = c(1, 2, 1)*s
  gamma3 = c(2, 1, 1)*s
  beta_share = c(rep(h, p_share), rep(0, p/2-p_share), rep(-h, p_share), rep(0, p/2-p_share))
  
  #  prefix = paste0('s=', s)
  DIFF.gamma = list()
  DIFF.beta = list()
  MSE.train = list()
  MSE.train.alpha = list()
  MSE.test = list()
  SS = list()
  
  for (ii in 1:2){
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    F1.train = F1
    F2.train = F2
    F3.train = F3
    
    U1 = mvrnorm(n, rep(0, p), para1$SigmaU)
    U2 = mvrnorm(n, rep(0, p), para2$SigmaU)
    U3 = mvrnorm(n, rep(0, p), para3$SigmaU)
    
    U1.train = U1
    U2.train = U2
    U3.train = U3
    
    X1 = F1%*%L1 + U1
    X2 = F2%*%L2 + U2
    X3 = F3%*%L3 + U3
    
    X.train.list = list(X1, X2, X3)
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    
    Y1 = u1 + F1%*%gamma1 + U1%*%beta_share + eps1
    Y2 = u2 + F2%*%gamma2 + U2%*%beta_share + eps2
    Y3 = u3 + F3%*%gamma3 + U3%*%beta_share + eps3
    
    Y.train.list = list(Y1, Y2, Y3)
    
    # testing
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    U1 = mvrnorm(n, rep(0, p), para1$SigmaU)
    U2 = mvrnorm(n, rep(0, p), para2$SigmaU)
    U3 = mvrnorm(n, rep(0, p), para3$SigmaU)
    
    X1 = F1%*%L1 + U1
    X2 = F2%*%L2 + U2
    X3 = F3%*%L3 + U3
    
    X.test.list = list(X1, X2, X3)
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    
    Y1 = u1 + F1%*%gamma1 + U1%*%beta_share + eps1
    Y2 = u2 + F2%*%gamma2 + U2%*%beta_share + eps2
    Y3 = u3 + F3%*%gamma3 + U3%*%beta_share + eps3
    
    Y.test.list = list(Y1, Y2, Y3)
    
    X.train = do.call(rbind, X.train.list)
    X.test = do.call(rbind, X.test.list)
    Y.train = do.call(c, Y.train.list)
    Y.test = do.call(c, Y.test.list)
    
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
    
    # global model
    ml.global = cv.glmnet(x = X.train, y = Y.train, alpha = 1, standardize = F)
    
    Yhat.global.train = predict(ml.global, s=ml.global$lambda.min, newx = X.train)
    mse.global.train.vec = sapply(label.level, function(l) 
      mean((Yhat.global.train[label.train==l] - Y.train[label.train==l])^2))
    mse.global.train = sum(mse.global.train.vec*n.train.vec)/sum(n.train.vec)
    
    Yhat.global.test = predict(ml.global, s=ml.global$lambda.min, newx = X.test)
    mse.global.test.vec = sapply(label.level, function(l) 
      mean((Yhat.global.test[label.test==l] - Y.test[label.test==l])^2))
    mse.global.test = sum(mse.global.test.vec*n.test.vec)/sum(n.test.vec)
    
    # class lasso
    ml.class = lapply(1:n_label, function(ix) 
      cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], standardize = F, alpha = 1))
    Yhat.class.train = lapply(1:n_label, function(ix) 
      predict(ml.class[[ix]], s=ml.class[[ix]]$lambda.min, newx = X.train.list[[ix]]))
    mse.class.train.vec = sapply(1:n_label, function(ix) 
      mean(((Yhat.class.train[[ix]])-(Y.train.list[[ix]]))^2))
    mse.class.train = sum(mse.class.train.vec*n.train.vec)/sum(n.train.vec)
    
    Yhat.class.test = lapply(1:n_label, function(ix) 
      predict(ml.class[[ix]], s=ml.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
    mse.class.test.vec = sapply(1:n_label, function(ix) 
      mean(((Yhat.class.test[[ix]])-(Y.test.list[[ix]]))^2))
    mse.class.test = sum(mse.class.test.vec*n.test.vec)/sum(n.test.vec)
    
    # --------------------------------------- ALPHA ---------------------------------------------
    X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], K = 5, plot = F))
    H.list = lapply(X2U.list, function(list) list$H)
    K.list = lapply(X2U.list, function(list) list$K)
    P.list = lapply(X2U.list, function(list) list$P)
    L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 
    
    F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
    U.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)
    
    FnU.test.list = lapply(1:n_label, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
    F.test.list = lapply(FnU.test.list, function(list) list$F_) 
    U.test.list = lapply(FnU.test.list, function(list) list$U) 
    
    # OLS.F
    data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                                  F.train.list[[ix]][,-1]))
    ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
    
    D1 = diag(eigen(X.train.list[[1]]%*%t(X.train.list[[1]]))$values[1:K.list[[1]]])
    H1 = solve(D1)%*%t(F.train.list[[1]][,-1])%*%F1.train%*%L1%*%t(L1)
    
    D2 = diag(eigen(X.train.list[[2]]%*%t(X.train.list[[2]]))$values[1:K.list[[2]]])
    H2 = solve(D2)%*%t(F.train.list[[2]][,-1])%*%F2.train%*%L2%*%t(L2)
    
    D3 = diag(eigen(X.train.list[[3]]%*%t(X.train.list[[3]]))$values[1:K.list[[3]]])
    H3 = solve(D3)%*%t(F.train.list[[3]][,-1])%*%F3.train%*%L3%*%t(L3)
    
    gammas.oracle = cbind(H1%*%gamma1, H2%*%gamma2, H3%*%gamma3)
    gammas.hat = cbind(ml.lm.F.list[[1]]$coefficients[-1], 
                       ml.lm.F.list[[2]]$coefficients[-1],
                       ml.lm.F.list[[3]]$coefficients[-1])
    
    diff.gamma = apply(gammas.oracle - gammas.hat, 2, function(diff) mean(diff^2))
    DIFF.gamma = list.append(DIFF.gamma, diff.gamma)
    
    # OLS.U
    U.train = do.call(rbind, U.train.list)
    HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
    HY.train = do.call(c, HY.train.list)
    
    # original
    ml.alpha = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = F)  
    beta.alpha = coef(ml.alpha, s=ml.alpha$lambda.min)[-1]
    beta1.group = coef(ml.class[[1]], s=ml.class[[1]]$lambda.min)[-1]
    beta2.group = coef(ml.class[[2]], s=ml.class[[2]]$lambda.min)[-1]
    beta3.group = coef(ml.class[[3]], s=ml.class[[3]]$lambda.min)[-1]
    beta.global = coef(ml.global, s=ml.global$lambda.min)[-1]
    betas = cbind(beta.global, beta1.group, beta2.group, beta3.group, beta.alpha)
    diff.beta = apply(betas - beta_share, 2, function(diff) mean(diff^2))
    DIFF.beta = list.append(DIFF.beta, diff.beta)
    
    PYhat.alpha.train.list = lapply(1:n_label, function(ix) ml.lm.F.list[[ix]]$fitted.values)
    HYhat.alpha.train.list = lapply(1:n_label, function(ix) 
      predict(ml.alpha, s=ml.alpha$lambda.min, U.train.list[[ix]]))
    Yhat.alpha.train.list = lapply(1:n_label, function(ix) PYhat.alpha.train.list[[ix]] + HYhat.alpha.train.list[[ix]])
    mse.alpha.train.list = compute.mse(Y.train.list, Yhat.alpha.train.list)
    
    SSF = sapply(1:n_label, function(ix) sum(PYhat.alpha.train.list[[ix]]^2))
    SSU = sapply(1:n_label, function(ix) sum(HYhat.alpha.train.list[[ix]]^2))
    SS = list.append(SS, rbind(SSF, SSU))
    
    # training data values
    F_.train.list = list(F1.train, F2.train, F3.train)
    U_.train.list = list(U1.train, U2.train, U3.train)
    gamma.list = list(gamma1, gamma2, gamma3)
    u.list = list(u1, u2, u3)
    
    mse1.alpha.train.vec = sapply(1:n_label, function(ix) 
      mean((u.list[[ix]]+F_.train.list[[ix]]%*%gamma.list[[ix]] - PYhat.alpha.train.list[[ix]])^2))
    mse2.alpha.train.vec = sapply(1:n_label, function(ix) 
      mean((U_.train.list[[ix]]%*%beta_share - HYhat.alpha.train.list[[ix]])^2))
    
    # predicting testing data
    PYhat.alpha.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    HYhat.alpha.test.list = lapply(1:n_label, function(ix) 
      predict(ml.alpha, s=ml.alpha$lambda.min, U.test.list[[ix]]))
    Yhat.alpha.test.list = lapply(1:n_label, function(ix) 
      PYhat.alpha.test.list[[ix]] + HYhat.alpha.test.list[[ix]])
    mse.alpha.test.list = compute.mse(Y.test.list, Yhat.alpha.test.list)  
    
    #------------------------------------------ALPHA-0----------------------------------------------------
    X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = 0, plot = F))
    H.list = lapply(X2U.list, function(list) list$H)
    K.list = lapply(X2U.list, function(list) list$K)
    P.list = lapply(X2U.list, function(list) list$P)
    L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 
    
    F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
    U.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)
    
    FnU.test.list = lapply(1:n_label, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
    F.test.list = lapply(FnU.test.list, function(list) list$F_) 
    U.test.list = lapply(FnU.test.list, function(list) list$U) 
    
    # OLS.F
    data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                                  F.train.list[[ix]][,-1]))
    ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
    
    # OLS.U
    U.train = do.call(rbind, U.train.list)
    HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
    HY.train = do.call(c, HY.train.list)
    
    ml.alpha = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = F)
    
    PYhat.alpha.train.list = lapply(1:n_label, function(ix) ml.lm.F.list[[ix]]$fitted.values)
    HYhat.alpha.train.list = lapply(1:n_label, function(ix) 
      predict(ml.alpha, s=ml.alpha$lambda.min, U.train.list[[ix]]))
    Yhat.alpha.train.list = lapply(1:n_label, function(ix) PYhat.alpha.train.list[[ix]] + HYhat.alpha.train.list[[ix]])
    mse.alpha0.train.list = compute.mse(Y.train.list, Yhat.alpha.train.list)
    
    # testing
    PYhat.alpha.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    HYhat.alpha.test.list = lapply(1:n_label, function(ix) 
      predict(ml.alpha, s=ml.alpha$lambda.min, U.test.list[[ix]]))
    Yhat.alpha.test.list = lapply(1:n_label, function(ix) 
      PYhat.alpha.test.list[[ix]] + HYhat.alpha.test.list[[ix]])
    mse.alpha0.test.list = compute.mse(Y.test.list, Yhat.alpha.test.list) 
    
    mse.train = rbind(mse.global.train.vec, mse.class.train.vec, mse.alpha0.train.list[[1]], mse.alpha.train.list[[1]])
    mse.train.alpha = rbind(mse1.alpha.train.vec, mse2.alpha.train.vec)
    mse.test = rbind(mse.global.test.vec, mse.class.test.vec, mse.alpha0.test.list[[1]], mse.alpha.test.list[[1]])
    
    row.names(mse.train) = NULL
    row.names(mse.train.alpha) = NULL
    row.names(mse.test) = NULL
    
    MSE.train = list.append(MSE.train, mse.train)
    MSE.train.alpha = list.append(MSE.train.alpha, mse.train.alpha)
    MSE.test = list.append(MSE.test, mse.test)
  }
  file.name = paste0("result_lasso_estimate_s=", format(s, nsmall = 1), ".RData")
  save(DIFF.gamma, DIFF.beta, MSE.train, MSE.train.alpha, MSE.test, SS, file = file.name)
}






