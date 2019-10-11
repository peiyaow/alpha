# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
# cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(caret)
library(glmnet)
require(methods)
# library(POET)

#path = "~/Documents/GitHub/alpha/"
path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "functions.R"))
source(paste0(path, "sim_func.R"))
load(paste0(path, "/sim/parameter/seeds.RData"))
load(paste0(path, "/sim/parameter/para_333.RData"))
load(paste0(path, "/sim/parameter/para_grid.RData"))

para = para_grid.list[[ii]]

dir_name = paste0(para$beta, "_", para$sigma.vec[3], "_", sprintf('%03d', para$p))
system(paste0("mkdir ", dir_name))
setwd(dir_name)

for (kk in 1: 50){
  set.seed(seeds[kk])
  n = 200
  p = para$p
  n_label = 3
  beta_F.list = lapply(1:n_label, function(ix) c(rep(ix, 5), rep(0, p-5)))
  # beta_F.list = lapply(1:n_label, function(ix) c(rep(1, 5), rep(0, p-5)))
  
  sigma.vec = para$sigma.vec
  
  beta_U = c(c(rep(para$beta, 5), rep(0, p-5)))
  # beta_U = c(rep(0, 5), rep(para$beta, 5), rep(0, p-10)))
  
  n.train.vec = c(n, n, n)
  n.test.vec = c(n, n, n)
  ix.vec = c(0, cumsum(n.train.vec))
  label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
  label.level = levels(label.test)
  
  #sim_FnU.list = lapply(1:n_label, function(ix) simulateFnU0(n, p, K.list[[ix]], mu_B.list[[ix]], Sigma_B.list[[ix]], Sigma_U))
  sim_FnU.list = lapply(1:n_label, function(ix) simulateFnU(n, p, K0.list[[ix]], L0.list[[ix]], Sigma_U))
  
  
  # true loading matrix
  L.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$L)
  gamma.list = lapply(1:n_label, function(ix) L.list[[ix]]%*%beta_F.list[[ix]])
  
  F.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$F_)
  U.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$U)
  
  
  F.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$F.test)
  U.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$U.test)
  X.train.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$X)
  X.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$X.test)
  tY.train.list = lapply(1:n_label, function(ix) F.list[[ix]]%*%gamma.list[[ix]] + U.list[[ix]]%*%beta_U)
  tY.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%gamma.list[[ix]] + U.test.list[[ix]]%*%beta_U)
  Y.train.list = lapply(1:n_label, function(ix) tY.train.list[[ix]] + rnorm(n, sd = sigma.vec[ix]))
  Y.test.list = lapply(1:n_label, function(ix) tY.test.list[[ix]] + rnorm(n, sd = sigma.vec[ix]))
  
  # -------------------------------- diagnostic ------------------------------------
  # sapply(1:n_label, function(ix) var(F.list[[ix]]%*%gamma.list[[ix]]))
  # sapply(1:n_label, function(ix) var(U.list[[ix]]%*%beta_U))
  # sapply(1:n_label, function(ix) var(tY.train.list[[ix]]))
  # sapply(1:n_label, function(ix) var(Y.train.list[[ix]]))
  # 
  # sapply(1:n_label, function(ix) var(F.test.list[[ix]]%*%gamma.list[[ix]]))
  # sapply(1:n_label, function(ix) var(U.test.list[[ix]]%*%beta_U))
  # sapply(1:n_label, function(ix) var(tY.test.list[[ix]]))
  # sapply(1:n_label, function(ix) var(Y.test.list[[ix]]))
  # ---------------------------------------------------------------------------------
  
  X.train = do.call(rbind, X.train.list)
  X.test = do.call(rbind, X.test.list)
  Y.train = do.call(c, Y.train.list)
  Y.test = do.call(c, Y.test.list)
  
  # standardize X (subtract mean)
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  data.train.list = lapply(1:n_label, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
  data.test.list = lapply(1:n_label, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
  
  # ------------------------------------ global model ---------------------------------------
  # global lm
  data.train = data.frame(Y=Y.train, X.train)
  ml.lm.global = lm(Y~., data = data.train)
  Yhat.lm.global.test = predict(ml.lm.global, new = data.frame(X.test))
  mse.lm.global.vec = sapply(label.level, function(l) mean((Yhat.lm.global.test[label.test==l] - Y.test[label.test==l])^2))
  mse.lm.global = sum(mse.lm.global.vec*n.test.vec)/sum(n.test.vec)
  
  # global ridge
  ml.ridge.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0)
  Yhat.ridge.global.test = predict(ml.ridge.global, s=ml.ridge.global$lambda.min, newx = X.test)
  mse.ridge.global.vec = sapply(label.level, function(l) mean((Yhat.ridge.global.test[label.test==l] - Y.test[label.test==l])^2))
  mse.ridge.global = sum(mse.ridge.global.vec*n.test.vec)/sum(n.test.vec)
  
  # global EN
  ml.EN.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0.5)
  Yhat.EN.global.test = predict(ml.EN.global, s=ml.EN.global$lambda.min, newx = X.test)
  mse.EN.global.vec = sapply(label.level, function(l) mean((Yhat.EN.global.test[label.test==l] - Y.test[label.test==l])^2))
  mse.EN.global = sum(mse.EN.global.vec*n.test.vec)/sum(n.test.vec)
  
  # global lasso
  ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 1)
  Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
  mse.lasso.global.vec = sapply(label.level, function(l) mean((Yhat.lasso.global.test[label.test==l] - Y.test[label.test==l])^2))
  mse.lasso.global = sum(mse.lasso.global.vec*n.test.vec)/sum(n.test.vec)
  # ------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------- WLS ---------------------------------------------
  Y.train.mean = lapply(Y.train.list, mean)
  Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]] - Y.train.mean[[l]]))
  X.train.WLS = do.call(rbind, X.train.list)
  data.train.WLS = data.frame(Y = Y.train.WLS, X.train.WLS)
  ml.lm.WLS = lm(Y~., data = data.train.WLS)
  sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.WLS$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
  w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
  
  # WLS lm
  ml.lm.WLS = lm(Y~., data = data.train.WLS, weights = w)
  Yhat.lm.WLS.test = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.test.list[[ix]])))
  mse.lm.WLS.vec = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
  mse.lm.WLS = sum(mse.lm.WLS.vec*n.test.vec)/sum(n.test.vec)
  
  # WLS ridge
  ml.ridge.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 0, weights = w)
  Yhat.ridge.WLS.test = lapply(1:n_label, function(ix) predict(ml.ridge.WLS, s=ml.ridge.WLS$lambda.min, newx = X.test.list[[ix]]))
  mse.ridge.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
  mse.ridge.WLS = sum(mse.ridge.WLS.vec*n.test.vec)/sum(n.test.vec)
  
  # WLS EN
  ml.EN.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 0.5, weights = w)
  Yhat.EN.WLS.test = lapply(1:n_label, function(ix) predict(ml.EN.WLS, s=ml.EN.WLS$lambda.min, newx = X.test.list[[ix]]))
  mse.EN.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
  mse.EN.WLS = sum(mse.EN.WLS.vec*n.test.vec)/sum(n.test.vec)
  
  # WLS lasso
  ml.lasso.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 1, weights = w)
  Yhat.lasso.WLS.test = lapply(1:n_label, function(ix) predict(ml.lasso.WLS, s=ml.lasso.WLS$lambda.min, newx = X.test.list[[ix]]))
  mse.lasso.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
  mse.lasso.WLS = sum(mse.lasso.WLS.vec*n.test.vec)/sum(n.test.vec)
  # -------------------------------------------------------------------------------------------
  
  # -------------------------------------- classwise ------------------------------------------
  # class lm
  ml.OLS.X.class = lapply(1:n_label, function(ix) lm(Y~., data = data.train.list[[ix]]))
  Yhat.OLS.X.class.test = lapply(1:n_label, function(ix) predict(ml.OLS.X.class[[ix]], new = data.frame(X.test.list[[ix]])))
  mse.OLS.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.OLS.X.class.test[[ix]]-(Y.test.list[[ix]]))^2)))
  mse.OLS.X.class = sum(mse.OLS.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # class ridge
  ml.ridge.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
  Yhat.ridge.X.class.test = lapply(1:n_label, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
  mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # class EN
  ml.EN.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0.5))
  Yhat.EN.X.class.test = lapply(1:n_label, function(ix) predict(ml.EN.X.class[[ix]], s=ml.EN.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.EN.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
  mse.EN.X.class = sum(mse.EN.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # class lasso
  ml.lasso.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 1))
  Yhat.lasso.X.class.test = lapply(1:n_label, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.lasso.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
  mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)
  # ---------------------------------------------------------------------------------------------
  
  # --------------------------------------- ALPHA ---------------------------------------------
  X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))
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
  
  POET.res = POET(t(U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
  SigmaU_hat = POET.res$SigmaU
  
  beta.OLS.U = solve(SigmaU_hat*sum(n.train.vec))%*%t(U.train)%*%HY.train
  HYhat.OLS.U = U.train%*%beta.OLS.U
  
  ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0)
  EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5)
  lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1)
  
  HYhat.test.OLS.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.OLS.U)
  HYhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
  HYhat.test.EN.OLS.list = lapply(1:n_label, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
  HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  
  Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
  Yhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
  Yhat.test.EN.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
  Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
  
  mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
  mse.ridge.OLS.list = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
  mse.EN.OLS.list = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
  mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
  
  
  # WLS
  # LD estimation for sigma
  sigma2.LD = sapply(1:n_label, function(ix) mean((HY.train.list[[ix]]-HYhat.OLS.U[(ix.vec[ix]+1):ix.vec[ix+1]])^2))
  w.LD = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2.LD[l]), n.train.vec[l])))
  
  # HD estimation for sigma
  var_moment_est.res = var_moment_est(HY.train, U.train, SigmaU_hat)
  tao2 = var_moment_est.res$B - var_moment_est.res$A
  sigma2.HD = sapply(1:n_label, function(ix) mean(HY.train.list[[ix]]^2)-tao2)
  w.HD = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2.HD[l]), n.train.vec[l])))
  
  #w.LD = w.LD/sum(w.LD)*sum(n.train.vec)
  #w.HD = w.HD/sum(w.HD)*sum(n.train.vec)
  
  W.LD = diag(w.LD)
  W.HD = diag(w.HD)
  
  # POET.res.LD = POET(t(sqrt(W.LD)%*%U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
  # SigmaU_hat.LD = POET.res.LD$SigmaU
  # 
  # POET.res.HD = POET(t(sqrt(W.HD)%*%U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
  # SigmaU_hat.HD = POET.res.HD$SigmaU
  # 
  # beta.WLS.U.LD = solve(SigmaU_hat.LD*sum(n.train.vec))%*%t(U.train)%*%W.LD%*%HY.train
  # beta.WLS.U.HD = solve(SigmaU_hat.HD*sum(n.train.vec))%*%t(U.train)%*%W.HD%*%HY.train
  
  beta.WLS.U.LD = solve((SigmaU_hat*sum(W.LD)))%*%t(U.train)%*%W.LD%*%HY.train
  beta.WLS.U.HD = solve((SigmaU_hat*sum(W.HD)))%*%t(U.train)%*%W.HD%*%HY.train
  
  # low-d penalized
  ridge.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 0)
  EN.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 0.5)
  lasso.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 1)
  
  HYhat.test.WLS.LD.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.WLS.U.LD)
  HYhat.test.ridge.WLS.LD.list = lapply(1:n_label, function(ix) predict(ridge.WLS.U.LD, s=ridge.WLS.U.LD$lambda.min, U.test.list[[ix]]))
  HYhat.test.EN.WLS.LD.list = lapply(1:n_label, function(ix) predict(EN.WLS.U.LD, s=EN.WLS.U.LD$lambda.min, U.test.list[[ix]]))
  HYhat.test.lasso.WLS.LD.list = lapply(1:n_label, function(ix) predict(lasso.WLS.U.LD, s=lasso.WLS.U.LD$lambda.min, U.test.list[[ix]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  
  Yhat.test.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.WLS.LD.list[[ix]])
  Yhat.test.ridge.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.WLS.LD.list[[ix]])
  Yhat.test.EN.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.WLS.LD.list[[ix]])
  Yhat.test.lasso.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.WLS.LD.list[[ix]])
  
  mse.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.WLS.LD.list)
  mse.ridge.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.ridge.WLS.LD.list)
  mse.EN.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.EN.WLS.LD.list)
  mse.lasso.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.lasso.WLS.LD.list)
  
  
  # high-d penalized
  if (n*3 > p){
    mse.WLS.HD.list = compute.mse(list(NA, NA, NA), list(NA, NA, NA))
    mse.ridge.WLS.HD.list = compute.mse(list(NA, NA, NA), list(NA, NA, NA))
    mse.EN.WLS.HD.list = compute.mse(list(NA, NA, NA), list(NA, NA, NA))
    mse.lasso.WLS.HD.list = compute.mse(list(NA, NA, NA), list(NA, NA, NA))
  }else{
    ridge.WLS.U.HD = cv.glmnet(x = U.train, y = HY.train, weight = w.HD, alpha = 0)
    EN.WLS.U.HD = cv.glmnet(x = U.train, y = HY.train, weight = w.HD, alpha = 0.5)
    lasso.WLS.U.HD = cv.glmnet(x = U.train, y = HY.train, weight = w.HD, alpha = 1)
    
    HYhat.test.WLS.HD.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.WLS.U.HD)
    HYhat.test.ridge.WLS.HD.list = lapply(1:n_label, function(ix) predict(ridge.WLS.U.HD, s=ridge.WLS.U.HD$lambda.min, U.test.list[[ix]]))
    HYhat.test.EN.WLS.HD.list = lapply(1:n_label, function(ix) predict(EN.WLS.U.HD, s=EN.WLS.U.HD$lambda.min, U.test.list[[ix]]))
    HYhat.test.lasso.WLS.HD.list = lapply(1:n_label, function(ix) predict(lasso.WLS.U.HD, s=lasso.WLS.U.HD$lambda.min, U.test.list[[ix]]))
    PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    
    Yhat.test.WLS.HD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.WLS.HD.list[[ix]])
    Yhat.test.ridge.WLS.HD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.WLS.HD.list[[ix]])
    Yhat.test.EN.WLS.HD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.WLS.HD.list[[ix]])
    Yhat.test.lasso.WLS.HD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.WLS.HD.list[[ix]])
    
    mse.WLS.HD.list = compute.mse(Y.test.list, Yhat.test.WLS.HD.list)
    mse.ridge.WLS.HD.list = compute.mse(Y.test.list, Yhat.test.ridge.WLS.HD.list)
    mse.EN.WLS.HD.list = compute.mse(Y.test.list, Yhat.test.EN.WLS.HD.list)
    mse.lasso.WLS.HD.list = compute.mse(Y.test.list, Yhat.test.lasso.WLS.HD.list)
  }
  
  
  
  # ------------------------------------------------------------------------------------------------
  
  file.name = "result_overall.csv"
  write.table(t(c(mse.lm.global, mse.ridge.global, mse.EN.global, mse.lasso.global, 
                  mse.lm.WLS, mse.ridge.WLS, mse.EN.WLS, mse.lasso.WLS, 
                  mse.OLS.X.class, mse.ridge.X.class, mse.EN.X.class, mse.lasso.X.class, 
                  mse.OLS.list[[2]], mse.ridge.OLS.list[[2]], mse.EN.OLS.list[[2]], mse.lasso.OLS.list[[2]],
                  mse.WLS.LD.list[[2]], mse.ridge.WLS.LD.list[[2]], mse.EN.WLS.LD.list[[2]], mse.lasso.WLS.LD.list[[2]],
                  mse.WLS.HD.list[[2]], mse.ridge.WLS.HD.list[[2]], mse.EN.WLS.HD.list[[2]], mse.lasso.WLS.HD.list[[2]])),
              file = file.name, sep = ',', append = T, col.names = F, row.names = F)
  
  file.name = "result_class.csv"
  write.table(t(c(mse.lm.global.vec, mse.ridge.global.vec, mse.EN.global.vec, mse.lasso.global.vec, 
                  mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.EN.WLS.vec, mse.lasso.WLS.vec, 
                  mse.OLS.X.class.vec, mse.ridge.X.class.vec, mse.EN.X.class.vec, mse.lasso.X.class.vec, 
                  mse.OLS.list[[1]], mse.ridge.OLS.list[[1]], mse.EN.OLS.list[[1]], mse.lasso.OLS.list[[1]],
                  mse.WLS.LD.list[[1]], mse.ridge.WLS.LD.list[[1]], mse.EN.WLS.LD.list[[1]], mse.lasso.WLS.LD.list[[1]],
                  mse.WLS.HD.list[[1]], mse.ridge.WLS.HD.list[[1]], mse.EN.WLS.HD.list[[1]], mse.lasso.WLS.HD.list[[1]])), 
              file = file.name, sep = ',', append = T, col.names = F, row.names = F)
}
