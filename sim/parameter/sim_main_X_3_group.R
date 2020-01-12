library(caret)
library(glmnet)
require(methods)

path = "~/Documents/GitHub/alpha/"

source(paste0(path, "functions.R"))
source(paste0(path, "sim_func.R"))
load(paste0(path, "/sim/parameter/seeds.RData"))

p = 200
n = 100
AA = matrix(, nrow = 0, ncol = 2)

for (h in seq(0, 2, by = 0.5)){
  for (ii in 1:5){
    p_share = 20
    mu_X = rep(0, p)
    Sigma_X1 = diag(c(c(40, 20, 4), exp(seq(log(0.3), log(0.01), length.out = p - 3))))/2
    Sigma_X2 = diag(c(c(30, 15, 5), exp(seq(log(0.4), log(0.01), length.out = p - 3))))/2
    Sigma_X3 = diag(c(c(50, 25, 5)/10*4, exp(seq(log(0.15), log(0.01), length.out = p - 3))))/2
    
    sigma.vec = c(3, 3, 3)
    
    beta1_unique = c(1, 1, 0)
    beta2_unique = c(1, 0, 1)
    beta3_unique = c(0, 1, 1)
    beta_share = rep(h, p_share)
    
    beta1 = c(beta1_unique, beta_share, rep(0, p-p_share-3))
    beta2 = c(beta2_unique, beta_share, rep(0, p-p_share-3))
    beta3 = c(beta3_unique, beta_share, rep(0, p-p_share-3))
    
    n_label = 3
    n.train.vec = c(n, n, n)
    n.test.vec = c(n, n, n)
    ix.vec = c(0, cumsum(n.train.vec))
    label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
    label.level = levels(label.test)
    
    X1 = mvrnorm(n, mu_X, Sigma_X1)
    X2 = mvrnorm(n, mu_X, Sigma_X2)
    X3 = mvrnorm(n, mu_X, Sigma_X3)
    
    X.train.list = list(X1, X2, X3)
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = X1%*%beta1+ eps1
    Y2 = X2%*%beta2+ eps2
    Y3 = X3%*%beta3+ eps3
    
    Y.train.list = list(Y1, Y2, Y3)
    
    X1 = mvrnorm(n, mu_X, Sigma_X1)
    X2 = mvrnorm(n, mu_X, Sigma_X2)
    X3 = mvrnorm(n, mu_X, Sigma_X3)
    
    X.test.list = list(X1, X2, X3)
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    
    Y1 = X1%*%beta1+ eps1
    Y2 = X2%*%beta2+ eps2
    Y3 = X3%*%beta3+ eps3
    
    Y.test.list = list(Y1, Y2, Y3)
    
    X.train = do.call(rbind, X.train.list)
    X.test = do.call(rbind, X.test.list)
    Y.train = do.call(c, Y.train.list)
    Y.test = do.call(c, Y.test.list)
    
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
    
    # global lasso
    ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 1)
    Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
    mse.lasso.global.vec = sapply(label.level, function(l) mean((Yhat.lasso.global.test[label.test==l] - Y.test[label.test==l])^2))
    mse.lasso.global = sum(mse.lasso.global.vec*n.test.vec)/sum(n.test.vec)
    
    # class lasso
    ml.lasso.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 1))
    # Yhat.lasso.X.class.test = lapply(1:n_label, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.train.list[[ix]]))
    # mse.lasso.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.train.list[[ix]]))^2))
    # 
    Yhat.lasso.X.class.test = lapply(1:n_label, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
    mse.lasso.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
    
    mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)
    
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
    
    # # WLS
    # # LD estimation for sigma
    # sigma2.LD = sapply(1:n_label, function(ix) mean((HY.train.list[[ix]]-HYhat.OLS.U[(ix.vec[ix]+1):ix.vec[ix+1]])^2))
    # w.LD = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2.LD[l]), n.train.vec[l])))
    # 
    # # HD estimation for sigma
    # var_moment_est.res = var_moment_est(HY.train, U.train, SigmaU_hat)
    # tao2 = var_moment_est.res$B - var_moment_est.res$A
    # sigma2.HD = sapply(1:n_label, function(ix) mean(HY.train.list[[ix]]^2)-tao2)
    # w.HD = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2.HD[l]), n.train.vec[l])))
    # 
    # #w.LD = w.LD/sum(w.LD)*sum(n.train.vec)
    # #w.HD = w.HD/sum(w.HD)*sum(n.train.vec)
    # 
    # W.LD = diag(w.LD)
    # W.HD = diag(w.HD)
    # 
    # # POET.res.LD = POET(t(sqrt(W.LD)%*%U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
    # # SigmaU_hat.LD = POET.res.LD$SigmaU
    # # 
    # # POET.res.HD = POET(t(sqrt(W.HD)%*%U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
    # # SigmaU_hat.HD = POET.res.HD$SigmaU
    # # 
    # # beta.WLS.U.LD = solve(SigmaU_hat.LD*sum(n.train.vec))%*%t(U.train)%*%W.LD%*%HY.train
    # # beta.WLS.U.HD = solve(SigmaU_hat.HD*sum(n.train.vec))%*%t(U.train)%*%W.HD%*%HY.train
    # 
    # beta.WLS.U.LD = solve((SigmaU_hat*sum(W.LD)))%*%t(U.train)%*%W.LD%*%HY.train
    # beta.WLS.U.HD = solve((SigmaU_hat*sum(W.HD)))%*%t(U.train)%*%W.HD%*%HY.train
    # 
    # # low-d penalized
    # ridge.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 0)
    # EN.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 0.5)
    # lasso.WLS.U.LD = cv.glmnet(x = U.train, y = HY.train, weight = w.LD, alpha = 1)
    # 
    # HYhat.test.WLS.LD.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.WLS.U.LD)
    # HYhat.test.ridge.WLS.LD.list = lapply(1:n_label, function(ix) predict(ridge.WLS.U.LD, s=ridge.WLS.U.LD$lambda.min, U.test.list[[ix]]))
    # HYhat.test.EN.WLS.LD.list = lapply(1:n_label, function(ix) predict(EN.WLS.U.LD, s=EN.WLS.U.LD$lambda.min, U.test.list[[ix]]))
    # HYhat.test.lasso.WLS.LD.list = lapply(1:n_label, function(ix) predict(lasso.WLS.U.LD, s=lasso.WLS.U.LD$lambda.min, U.test.list[[ix]]))
    # PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    # 
    # Yhat.test.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.WLS.LD.list[[ix]])
    # Yhat.test.ridge.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.WLS.LD.list[[ix]])
    # Yhat.test.EN.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.WLS.LD.list[[ix]])
    # Yhat.test.lasso.WLS.LD.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.WLS.LD.list[[ix]])
    # 
    # mse.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.WLS.LD.list)
    # mse.ridge.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.ridge.WLS.LD.list)
    # mse.EN.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.EN.WLS.LD.list)
    # mse.lasso.WLS.LD.list = compute.mse(Y.test.list, Yhat.test.lasso.WLS.LD.list)
    
    # K.list
    # mse.lasso.global.vec
    # mse.lasso.X.class.vec
    # mse.lasso.X.class
    # mse.lasso.OLS.list
    # mse.lasso.WLS.LD.list
    print(c(mse.lasso.X.class, mse.lasso.OLS.list$mse))
    AA = rbind(AA, c(mse.lasso.X.class, mse.lasso.OLS.list$mse))
  }
}





