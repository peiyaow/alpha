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
n_label = 3
si = 2 # lm error sd

# control parameter for beta
# s = 1 # unique
h = 1 # share

ds = sqrt(p) # sqrt(p) multiplier parameter for spike L 
rho = 1/10 # determine spike smaller than 1/K
rrho = 0.1 # correlation for the core sigma matrix
d = .05 # uniform range

spike_factor1 = c(1, sapply(2:K, function(k) 1/gamma(k + 1)))*.7
spike_factor2 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.25)))
spike_factor3 = c(1, sapply(2:K, function(k) 1/gamma(k + 1.5)))*1.3

# plot(spike_factor2[1:(K-1)]/spike_factor2[2:K])

spike1 = spike_factor1 * ds
spike2 = spike_factor2 * ds
spike3 = spike_factor3 * ds

du = (spike1[K]*rho+spike2[K]*rho + spike3[K]*rho)/3

# para1 = FactorModelPara(n, p, K, spike = spike1, d, du = spike1[K]*rho, rrho)
# para2 = FactorModelPara(n, p, K, spike = spike2, d, du = spike2[K]*rho, rrho)
# para3 = FactorModelPara(n, p, K, spike = spike3, d, du = spike3[K]*rho, rrho)

para1 = FactorModelPara(n, p, K, spike = spike1, d, du = du, rrho)
para2 = FactorModelPara(n, p, K, spike = spike2, d, du = du, rrho)
para3 = FactorModelPara(n, p, K, spike = spike3, d, du = du, rrho)

u1 = 1
u2 = 2
u3 = 3

sigma.vec = rep(si, n_label)

n.train.vec = c(n, n, n)
n.test.vec = c(n, n, n)
ix.vec = c(0, cumsum(n.train.vec))
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

# # Lbeta norm
# sqrt(mean((para1$L%*%beta1)^2))
# sqrt(mean((para2$L%*%beta2)^2))
# sqrt(mean((para3$L%*%beta3)^2))

# generate data
DD = matrix(, nrow = 0, ncol = 4)
for (s in seq(0, 2.5, by = 0.1)){
  print(s)
  for (ii in 1:50){
    beta1_unique = c(1, 1, -1)*s
    beta2_unique = c(1, -1, 1)*s
    beta3_unique = c(-1, 1, 1)*s
    
    beta_share = c(rep(h, p_share), rep(-h, p_share))
    
    beta1 = c(beta1_unique, beta_share, rep(0, p-2*p_share-3))
    beta2 = c(beta2_unique, beta_share, rep(0, p-2*p_share-3))
    beta3 = c(beta3_unique, beta_share, rep(0, p-2*p_share-3))
    
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    U1 = mvrnorm(n, rep(0, p), para1$SigmaU)
    U2 = mvrnorm(n, rep(0, p), para2$SigmaU)
    U3 = mvrnorm(n, rep(0, p), para3$SigmaU)
    
    X1 = F1%*%para1$L + U1
    X2 = F2%*%para2$L + U2
    X3 = F3%*%para3$L + U3
    
    X.train.list = list(X1, X2, X3)
    X.train.mean = lapply(X.train.list, colMeans)
    X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = u1 + X1%*%beta1+ eps1
    Y2 = u2 + X2%*%beta2+ eps2
    Y3 = u3 + X3%*%beta3+ eps3
    
    Y.train.list = list(Y1, Y2, Y3)
    
    # testing
    F1 = mvrnorm(n, rep(0, K), diag(K))
    F2 = mvrnorm(n, rep(0, K), diag(K))
    F3 = mvrnorm(n, rep(0, K), diag(K))
    
    U1 = mvrnorm(n, rep(0, p), para1$SigmaU)
    U2 = mvrnorm(n, rep(0, p), para2$SigmaU)
    U3 = mvrnorm(n, rep(0, p), para3$SigmaU)
    
    X1 = F1%*%para1$L + U1
    X2 = F2%*%para2$L + U2
    X3 = F3%*%para3$L + U3
    
    X.test.list = list(X1, X2, X3)
    
    eps1 = rnorm(n, sd = sigma.vec[1])
    eps2 = rnorm(n, sd = sigma.vec[2])
    eps3 = rnorm(n, sd = sigma.vec[3])
    Y1 = u1 + X1%*%beta1+ eps1
    Y2 = u2 + X2%*%beta2+ eps2
    Y3 = u3 + X3%*%beta3+ eps3
    
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
    # ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0, standardize = F, intercept = T)
    ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0, standardize = F, intercept = T, lambda = exp(log(seq(p*2, 1))))
    Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
    mse.lasso.global.vec = sapply(label.level, function(l) mean((Yhat.lasso.global.test[label.test==l] - Y.test[label.test==l])^2))
    mse.lasso.global = sum(mse.lasso.global.vec*n.test.vec)/sum(n.test.vec)
    
    # class lasso
    ml.lasso.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], standardize = F, alpha = 0, lambda = exp(log(seq(p*2, 1)))))
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
    
    # POET.res = POET(t(U.train), K = 0, C = 0.5, thres = "soft", matrix = "vad")
    # SigmaU_hat = POET.res$SigmaU
    
    # beta.OLS.U = solve(SigmaU_hat*sum(n.train.vec))%*%t(U.train)%*%HY.train
    # HYhat.OLS.U = U.train%*%beta.OLS.U
    
    ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0)
    EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5)
    lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1)
    
    # HYhat.test.OLS.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.OLS.U)
    HYhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
    HYhat.test.EN.OLS.list = lapply(1:n_label, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
    HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
    PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    
    # Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
    Yhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
    Yhat.test.EN.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
    Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
    
    # mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
    mse.ridge.OLS.list = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
    mse.EN.OLS.list = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
    mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
    
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
    
    ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0)
    EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5)
    lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1)
    
    # HYhat.test.OLS.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.OLS.U)
    HYhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
    HYhat.test.EN.OLS.list = lapply(1:n_label, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
    HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
    PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
    
    # Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
    Yhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
    Yhat.test.EN.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
    Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])
    
    # mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
    mse.ridge.OLS.list0 = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
    mse.EN.OLS.list0 = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
    mse.lasso.OLS.list0 = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)
    
    print(c(mse.lasso.global, mse.lasso.X.class, mse.ridge.OLS.list0$mse, mse.ridge.OLS.list$mse))
    
    DD = rbind(DD, c(mse.lasso.global, mse.lasso.X.class, mse.ridge.OLS.list0[[2]], mse.ridge.OLS.list[[2]]))
  }
}

save(DD, file = "DD.RData")

