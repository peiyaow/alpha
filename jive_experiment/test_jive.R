jive.list = myjive(X.train.list)
L = jive.list$L
K = jive.list$K
F.train.list = jive.list$F.list
U.train.list = jive.list$U.list
R.train.list = jive.list$R.list
P.list = jive.list$P.list

FnU.test.list = lapply(1:n_label, function(ix) FnU.svd2(X.test.list[[ix]], L)) 
F.test.list = lapply(FnU.test.list, function(list) list$F) 
U.test.list = lapply(FnU.test.list, function(list) list$U)







Y.train.mean = lapply(Y.train.list, mean)
gamma_new = matrix(0, ncol = n_label, nrow = K)

mse.new.list = list()

Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]]- Fgamma.train.list[[l]] - Y.train.mean[[l]]))
data.train.WLS = data.frame(Y = Y.train.WLS, X.train)
ml.lm.WLS = lm(Y~., data = data.train.WLS)
Yhat.lm.WLS.train = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.train.list[[ix]])))

res.Y.train.list = lapply(1:n_label, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]] - Yhat.lm.WLS.train[[ix]])
mse.lm.WLS.train.vec = sapply(1:n_label, function(ix) mean((res.Y.train.list[[ix]]-Fgamma.train.list[[ix]])^2))
mse.new.list = list.append(mse.new.list, mse.lm.WLS.train.vec)

#gamma = return_gamma_mtx_lasso(res.Y.train.list, F.train.list, .001)$gamma
gamma = gamma_new
gamma_new = return_gamma_mtx_ridge(res.Y.train.list, F.train.list, .1)$gamma
print(gamma_new)

Fgamma.train.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%gamma_new[,ix])
mse.lm.WLS.train.vec = sapply(1:n_label, function(ix) mean((res.Y.train.list[[ix]]-Fgamma.train.list[[ix]])^2))
mse.new.list = list.append(mse.new.list, mse.lm.WLS.train.vec)





Yhat.lm.WLS.test = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.test.list[[ix]])))
mse.lm.WLS.vec11 = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS11 = sum(mse.lm.WLS.vec11*n.test.vec)/sum(n.test.vec)
mse.new.list = list.append(mse.new.list, mse.lm.WLS.vec11)

Fgamma_hat.test = lapply(1:n_label, function(ix) (F.test.list[[ix]])%*%gamma[,ix])
mse.lm.WLS.vec1 = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]+Fgamma_hat.test[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS1 = sum(mse.lm.WLS.vec1*n.test.vec)/sum(n.test.vec)
mse.new.list = list.append(mse.new.list, mse.lm.WLS.vec1)
mse.new.list

