# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------

library(glmnet)
require(methods)

# path = "~/Documents/GitHub/alpha/"
path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "/function/main_function.R"))
source(paste0(path, "/function/sim_function.R"))
source(paste0(path, "/function/glmnet_POET.R"))

set.seed(myseed)

p = 200
n = 100
K = 3

n_C = 21
n_lambda = 50

p_share = 10 # ridge 80 lasso 10
n_label = 3
si = 2

# control parameter for beta
h = 2 # share

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
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

for (s in c(0, 2, 4)){
  gamma1 = c(1, 1, 2)*s
  gamma2 = c(1, 2, 1)*s
  gamma3 = c(2, 1, 1)*s
  
  beta_share = c(rep(h, p_share), rep(0, p/2-p_share), rep(-h, p_share), rep(0, p/2-p_share))
  
  F1 = mvrnorm(n, rep(0, K), diag(K))
  F2 = mvrnorm(n, rep(0, K), diag(K))
  F3 = mvrnorm(n, rep(0, K), diag(K))
  
  U1 = mvrnorm(n, rep(0, p), para1$SigmaU)
  U2 = mvrnorm(n, rep(0, p), para2$SigmaU)
  U3 = mvrnorm(n, rep(0, p), para3$SigmaU)
  
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
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
  X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  
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
  
  # OLS.U
  U.train = do.call(rbind, U.train.list)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  HY.train = do.call(c, HY.train.list)
  
  # initialize lambda
  lambda.max = glmnet.lambda_max(U.train, HY.train, 1)
  epsilon <- .01
  lambdapath <- exp(seq(log(lambda.max), log(lambda.max*epsilon), 
                        length.out = n_lambda)) # 50
  C_path = seq(0, 1, length.out = n_C) # 21
  
  MSE = matrix(, nrow = 0, ncol = 3)
  for (lambda in lambdapath){
    for (C in C_path){
      SigmaU = POET(t(U.train), K = 0, C = C, thres = "soft", matrix = "vad")$SigmaU
      scout.beta = scout_POET(U.train, HY.train, SigmaU, lambda = lambda, alpha = 1)
      HYhat.test.OLS.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%scout.beta)
      PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
      Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
      mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
      MSE = rbind(MSE, mse.OLS.list$mse.vec)
    }
  }
  file.name = paste0("result_lasso_C_s=", format(s, nsmall = 1), ".csv")
  write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

  # scout POET
  cv.ml = cv.scout_POET(U.train, HY.train, 
                        C = C_path, 
                        lambda = lambdapath, 
                        alpha = 1, nfolds = 10)
  MSE.scout_POET = MSE[cv.ml$best.C.ix + n_C*(cv.ml$best.lambda.ix-1), ]
  
  # fix C = 0 scout
  cv.ml_C_0 = cv.scout_POET(U.train, HY.train, 
                            C = c(0), 
                            lambda = lambdapath, 
                            alpha = 1, nfolds = 10)
  MSE.scout_POET_C_0 = MSE[1 + n_C*(cv.ml_C_0$best.lambda.ix-1), ]
  
  # cv glmnet
  lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = F,
                          lambda = lambdapath)
  
  HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  
  Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
  
  mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
  
  MSE.cv = rbind(MSE.scout_POET, MSE.scout_POET_C_0, mse.lasso.OLS.list$mse.vec)
  
  file.name = paste0("result_lassocv_C_s=", format(s, nsmall = 1), ".csv")
  write.table(t(c(cv.ml$best.C.ix, cv.ml$best.lambda.ix, as.vector(MSE.cv))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
}
  