library(POET)
# OLS.U
POET.list = lapply(1:n_label, function(ix) POET(t(X.train.list[[ix]]), K = K.list[[ix]], C = 0.5, thres = "soft", matrix = "vad"))
Sigma_U.list = lapply(POET.list, function(list) list$SigmaU)
Sigma_U.list = lapply(1:n_label, function(ix) Sigma_U.list[[ix]]*n.train.vec[ix])
SigmaU_hat = Reduce('+', Sigma_U.list)/sum(n.train.vec)
HYhat = U.train%*%solve(SigmaU_hat*n*n_label)%*%t(U.train)%*%HY.train
# compute weights from OLS.F and OLS.U 
Yhat.lm.U.train.list = lapply(1:n_label, function(ix) HYhat[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l]))) 


W = diag(w)
sW.list = lapply(1:n_label, function(ix) sqrt(W)[(ix.vec[ix]+1):ix.vec[ix+1], (ix.vec[ix]+1):ix.vec[ix+1]]) #sqrt W
POET.list = lapply(1:n_label, function(ix) POET(t(sW.list[[ix]]%*%U.train.list[[ix]]), K = 0, C = 0.5, thres = "soft", matrix = "vad"))
Sigma_U.list = lapply(POET.list, function(list) list$SigmaU)
Sigma_U.list = lapply(1:n_label, function(ix) Sigma_U.list[[ix]]*n.train.vec[ix])
SigmaU_hat = Reduce('+', Sigma_U.list)/sum(n.train.vec)
HYhat = U.train%*%solve(SigmaU_hat*n*n_label)%*%t(U.train)%*%W%*%HY.train

