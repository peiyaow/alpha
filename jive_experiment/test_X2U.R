source("function_opt.R")

K = 3
X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = K, Plot = F))
L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 

F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_[,-1])
FnU.test.list = lapply(1:n_label, function(ix) FnU.svd2(X.test.list[[ix]], L.list[[ix]])) 
F.test.list = lapply(FnU.test.list, function(list) list$F) 
U.test.list = lapply(FnU.test.list, function(list) list$U)







Y.train.mean = lapply(Y.train.list, mean)
gamma_new = matrix(0, ncol = n_label, nrow = K)
Fgamma.train.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%gamma_new[,ix])
mse.lm.WLS.train.vec_prev = mse.lm.WLS.train.vec
#ix_cur = 1

mse.new.list = list()

Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]]- Fgamma.train.list[[l]] - Y.train.mean[[l]]))
plot(Y.train.WLS)
data.train.WLS = data.frame(Y = Y.train.WLS, X.train)
ml.lm.WLS = lm(Y~., data = data.train.WLS)
Yhat.lm.WLS.train = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.train.list[[ix]])))

# WLS lm
ix.vec = c(0,cumsum(n.train.vec))
sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.WLS$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
ml.lm.WLS = lm(Y~., data = data.train.WLS, weights = w)
Yhat.lm.WLS.train = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.train.list[[ix]])))

res.Y.train.list = lapply(1:n_label, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]] - Yhat.lm.WLS.train[[ix]])
plot(do.call(c, res.Y.train.list))
mse.lm.WLS.train.vec = sapply(1:n_label, function(ix) mean((res.Y.train.list[[ix]]-Fgamma.train.list[[ix]])^2))
mse.new.list = list.append(mse.new.list, mse.lm.WLS.train.vec)

#gamma = return_gamma_mtx_lasso(res.Y.train.list, F.train.list, .001)$gamma
gamma = gamma_new
#mse_diff = mse.lm.WLS.train.vec - mse.lm.WLS.train.vec_prev
#weight = as.numeric(mse_diff>0)
#weight = exp(mse_diff*10)/sum(exp(mse_diff*10))
#weight = c(1/4,1/4,1/4,1/4)
#weight = c(1,1,1,1)
# gamma_new = return_gamma_mtx_ridge(res.Y.train.list, F.train.list, .001)$gamma
gamma_new = return_gamma_mtx_ridge(res.Y.train.list, F.train.list, weight = 1/mse.lm.WLS.train.vec_prev, 0.1)$gamma
print(gamma_new)

Fgamma.train.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%gamma_new[,ix])
mse.lm.WLS.train.vec_prev = sapply(1:n_label, function(ix) mean((res.Y.train.list[[ix]]-Fgamma.train.list[[ix]])^2))
mse.new.list = list.append(mse.new.list, mse.lm.WLS.train.vec_prev)





Yhat.lm.WLS.test = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.test.list[[ix]])))
mse.lm.WLS.vec11 = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS11 = sum(mse.lm.WLS.vec11*n.test.vec)/sum(n.test.vec)
mse.new.list = list.append(mse.new.list, mse.lm.WLS.vec11)

Fgamma_hat.test = lapply(1:n_label, function(ix) (F.test.list[[ix]])%*%gamma[,ix])

mse.lm.WLS.vec1 = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]+Fgamma_hat.test[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS1 = sum(mse.lm.WLS.vec1*n.test.vec)/sum(n.test.vec)
mse.new.list = list.append(mse.new.list, mse.lm.WLS.vec1)
mse.new.list

