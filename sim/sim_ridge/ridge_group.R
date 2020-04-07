library(glmnet)
require(methods)

# path = "~/Documents/GitHub/alpha/"
path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "/function/main_function.R"))
source(paste0(path, "/function/sim_function.R"))

set.seed(706)

p = 200
n = 100
K = 3

p_share = 80 # ridge 80 lasso 5
si = 2 # lm error sd

h = 1 # share

ds = 10 # sqrt(p) multiplier parameter for spike L 
rho = 1/40 # determine spike smaller than 1/K

L.list = lapply(1:K, function(k) randortho(p, type = "orthonormal")[1:K,])
spike1 = c(1, sapply(2:K, function(k) 1/gamma(k + 1)))*.7*ds
spike2 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.25)))*ds
spike3 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.5)))*1.3*ds

L1 = L.list[[1]]*sqrt(spike1)
L2 = L.list[[2]]*sqrt(spike2)
L3 = L.list[[3]]*sqrt(spike3)

du = (spike1[K]*rho+spike2[K]*rho + spike3[K]*rho)/3
########
rrho = 0.1
R = matrix(0, nrow = p, ncol = p)
for (i in 1:p){
  for (j in seq(max(1,i-2), min(i+2, p), by = 1)){
    R[i,j] = rrho^(abs(i-j))
  }
}
########
SigmaU = du*R
# SigmaU = du*diag(p)

u1 = 1
u2 = 2
u3 = 3

sigma.vec = rep(si, K)

n.train.vec = c(n, n, n)
n.test.vec = c(n, n, n)
ix.vec = c(0, cumsum(n.train.vec))
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

DD = matrix(, nrow = 0, ncol = 4)
for (s in seq(0, 6, length.out = 50)){
  print(s)
  for (ii in 1:50){
    # beta1 = c(c(1, 1, 2)*s, rep(h, p_share), rep(0, p/2-p_share-K), c(1, 1, 2)*s, rep(-h, p_share), rep(0, p/2-p_share-K))
    # beta2 = c(c(1, 2, 1)*s, rep(h, p_share), rep(0, p/2-p_share-K), c(1, 2, 1)*s, rep(-h, p_share), rep(0, p/2-p_share-K))
    # beta3 = c(c(2, 1, 1)*s, rep(h, p_share), rep(0, p/2-p_share-K), c(2, 1, 1)*s, rep(-h, p_share), rep(0, p/2-p_share-K))
    
    a = 10
    beta1 = c(c(a, a, -a)*s, rep(h, p_share), rep(0, p-p_share-K))/sqrt(p)
    beta2 = c(c(a, -a, a)*s, rep(h, p_share), rep(0, p-p_share-K))/sqrt(p)
    beta3 = c(c(-a, a, a)*s, rep(h, p_share), rep(0, p-p_share-K))/sqrt(p)
    
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    U1 = mvrnorm(n, rep(0, p), SigmaU)
    U2 = mvrnorm(n, rep(0, p), SigmaU)
    U3 = mvrnorm(n, rep(0, p), SigmaU)
    
    X1 = F1%*%L1 + U1
    X2 = F2%*%L2 + U2
    X3 = F3%*%L3 + U3
    
    X.train.list = list(X1, X2, X3)
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.list = lapply(1:K, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = u1 + X1%*%beta1 + eps1
    Y2 = u2 + X2%*%beta2 + eps2
    Y3 = u3 + X3%*%beta3 + eps3
    
    Y.train.list = list(Y1, Y2, Y3)
    
    # testing
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    U1 = mvrnorm(n, rep(0, p), SigmaU)
    U2 = mvrnorm(n, rep(0, p), SigmaU)
    U3 = mvrnorm(n, rep(0, p), SigmaU)
    
    X1 = F1%*%L1 + U1
    X2 = F2%*%L2 + U2
    X3 = F3%*%L3 + U3
    
    X.test.list = list(X1, X2, X3)
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = u1 + X1%*%beta1 + eps1
    Y2 = u2 + X2%*%beta2 + eps2
    Y3 = u3 + X3%*%beta3 + eps3
    
    Y.test.list = list(Y1, Y2, Y3)
    
    X.train = do.call(rbind, X.train.list)
    X.test = do.call(rbind, X.test.list)
    Y.train = do.call(c, Y.train.list)
    Y.test = do.call(c, Y.test.list)
    
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
    X.train.list = lapply(1:K, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    X.test.list = lapply(1:K, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
    
    # global lasso
    ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0, intercept = T, standardize = F, lambda = exp(seq(log(p*2), log(.01), length.out = 100)))
    # ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0, standardize = T, intercept = T, lambda = exp(log(seq(p*2, 1))))
    Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
    mse.lasso.global.vec = sapply(label.level, function(l) mean((Yhat.lasso.global.test[label.test==l] - Y.test[label.test==l])^2))
    mse.lasso.global = sum(mse.lasso.global.vec*n.test.vec)/sum(n.test.vec)
    
    # class lasso
    ml.lasso.X.class = lapply(1:K, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0, standardize = F, lambda = exp(seq(log(p*2), log(.01), length.out = 100))))
    Yhat.lasso.X.class.test = lapply(1:K, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
    mse.lasso.X.class.vec = sapply(1:K, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
    mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)
    
    # --------------------------------------- ALPHA ---------------------------------------------
    X2U.list = lapply(1:K, function(ix) X2U1(X.train.list[[ix]], K = 5, plot = F))
    H.list = lapply(X2U.list, function(list) list$H)
    K.list = lapply(X2U.list, function(list) list$K)
    markK = K.list
    P.list = lapply(X2U.list, function(list) list$P)
    L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 
    
    F.train.list = lapply(1:K, function(ix) X2U.list[[ix]]$F_)
    U.train.list = lapply(1:K, function(ix) X2U.list[[ix]]$U)
    
    FnU.test.list = lapply(1:K, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
    F.test.list = lapply(FnU.test.list, function(list) list$F_) 
    U.test.list = lapply(FnU.test.list, function(list) list$U) 
    
    # OLS.F
    data.F.train.list = lapply(1:K, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                            F.train.list[[ix]][,-1]))
    ml.lm.F.list = lapply(1:K, function(l) lm(Y~., data = data.F.train.list[[l]]))
    
    # OLS.U
    U.train = do.call(rbind, U.train.list)
    HY.train.list = lapply(1:K, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
    HY.train = do.call(c, HY.train.list)
    
    # POET.res = POET(t(U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
    # SigmaU_hat = POET.res$SigmaU
    
    lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = F, lambda = exp(seq(log(p*2), log(.01), length.out = 100)))
    ml = lasso.OLS.U
    
    HYhat.test.lasso.OLS.list = lapply(1:K, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
    PYhat.test.list = lapply(1:K, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    
    Yhat.test.lasso.OLS.list = lapply(1:K, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
    
    mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
    
    #------------------------------------------ALPHA-0----------------------------------------------------
    X2U.list = lapply(1:K, function(ix) X2U2(X.train.list[[ix]], K = 0, plot = F))
    H.list = lapply(X2U.list, function(list) list$H)
    K.list = lapply(X2U.list, function(list) list$K)
    P.list = lapply(X2U.list, function(list) list$P)
    L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 
    
    F.train.list = lapply(1:K, function(ix) X2U.list[[ix]]$F_)
    U.train.list = lapply(1:K, function(ix) X2U.list[[ix]]$U)
    
    FnU.test.list = lapply(1:K, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
    F.test.list = lapply(FnU.test.list, function(list) list$F_) 
    U.test.list = lapply(FnU.test.list, function(list) list$U) 
    
    # OLS.F
    data.F.train.list = lapply(1:K, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                            F.train.list[[ix]][,-1]))
    ml.lm.F.list = lapply(1:K, function(l) lm(Y~., data = data.F.train.list[[l]]))
    
    # OLS.U
    U.train = do.call(rbind, U.train.list)
    HY.train.list = lapply(1:K, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
    HY.train = do.call(c, HY.train.list)
    
    lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = F, lambda = exp(seq(log(p*2), log(.01), length.out = 100)))
    ml0 = lasso.OLS.U
    
    HYhat.test.lasso.OLS.list = lapply(1:K, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
    PYhat.test.list = lapply(1:K, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    
    Yhat.test.lasso.OLS.list = lapply(1:K, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
    
    mse.lasso.OLS.list0 = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
    print(c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list0$mse, mse.lasso.OLS.list$mse))
    
    DD = rbind(DD, c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list0[[2]], mse.lasso.OLS.list[[2]]))
  }
}

save(DD, file = "DD.RData")