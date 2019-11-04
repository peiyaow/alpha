library(caret)
library(glmnet)
require(methods)
library(pracma)

# path = "~/Documents/GitHub/alpha/"
path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "functions.R"))
source(paste0(path, "sim_func.R"))
# load(paste0(path, "/sim/parameter/seeds.RData"))

set.seed(706)

p = 200
n = 100
p_share = 10
K = 3
n_label = 3
si = 5

# control parameter for beta
# s = 6 # unique
h = 3 # share

ds = .2 # dilute parameter for spike
rho = 1/4
d = 0.1

spike1 = c(20, 10, 10/3)*ds
spike2 = c(15, 7.5, 2.5)*ds
spike3 = c(10, 5, 5/3)*ds

sigma.vec = rep(si, n_label)

DD = matrix(, nrow = 0, ncol = 3)
for (s in seq(0, 6, by = 0.5)){
  print(s)
  for (ii in 1:50){
    beta1_unique = c(1, 1, 0)*s
    beta2_unique = c(1, 0, 1)*s
    beta3_unique = c(0, 1, 1)*s
    beta_share = rep(h, p_share)
    
    beta1 = c(beta1_unique, beta_share, rep(0, p-p_share-3))
    beta2 = c(beta2_unique, beta_share, rep(0, p-p_share-3))
    beta3 = c(beta3_unique, beta_share, rep(0, p-p_share-3))
    
    n.train.vec = c(n, n, n)
    n.test.vec = c(n, n, n)
    ix.vec = c(0, cumsum(n.train.vec))
    label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
    label.level = levels(label.test)
    
    X1 = generateX(K, p, spike = spike1, d, du = spike1[K]*rho)$X
    X2 = generateX(K, p, spike = spike2, d, du = spike2[K]*rho)$X
    X3 = generateX(K, p, spike = spike3, d, du = spike3[K]*rho)$X
    
    X.train.list = list(X1, X2, X3)
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = X1%*%beta1+ eps1
    Y2 = X2%*%beta2+ eps2
    Y3 = X3%*%beta3+ eps3
    
    Y.train.list = list(Y1, Y2, Y3)
    
    X1 = generateX(K, p, spike = spike1, d*10/3, du = 1/3*rho)$X
    X2 = generateX(K, p, spike = spike2, d*2.5, du = 0.25*rho)$X
    X3 = generateX(K, p, spike = spike3, d*5/3, du = 1/6*rho)$X
    
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
    X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
    
    # global lasso
    ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 1, standardize = T)
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
    
    print(c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list$mse))
    DD = rbind(DD, c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list[[2]]))
  }
}

save(DD, file = "DD.RData")
