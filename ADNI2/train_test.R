PYhat.train.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
PYhat.ridge.train.list = lapply(1:n_label, function(ix) predict(ml.ridge.F.list[[ix]], s = ml.ridge.F.list[[ix]]$lambda.min, newx = F.train.list[[ix]]))
HYhat.train.list = lapply(1:n_label, function(ix) cbind(1,U.train.list[[ix]])%*%ml.lm.U$coefficients)
Yhat.train.list = lapply(1:n_label, function(ix) PYhat.train.list[[ix]] + HYhat.train.list[[ix]])
Yhat.ridge.train.list = lapply(1:n_label, function(ix) PYhat.ridge.train.list[[ix]] + HYhat.train.list[[ix]])
mse.train.vec = sapply(1:n_label, function(ix) mean((Yhat.train.list[[ix]]-Y.train.list[[ix]])^2))
mse.ridge.train.vec = sapply(1:n_label, function(ix) mean((Yhat.ridge.train.list[[ix]]-Y.train.list[[ix]])^2))

Yhat.lm.WLS.train = lapply(1:n_label, function(ix) ml.lm.WLS$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]]+Y.train.mean[[ix]])

diff.PYhat.train.vec = sapply(1:n_label, function(ix) mean((PYhat.WLSbeta.train.list[[ix]] - PYhat.train.list[[ix]])^2))
diff.HYhat.train.vec = sapply(1:n_label, function(ix) mean((HYhat.WLSbeta.train.list[[ix]] - HYhat.train.list[[ix]])^2))

Yhat.WLSbeta.train.list = lapply(1:n_label, function(ix) PYhat.WLSbeta.train.list[[ix]] + HYhat.WLSbeta.train.list[[ix]])
mse.WLSbeta.train.vec = sapply(1:n_label, function(ix) mean((Yhat.WLSbeta.train.list[[ix]]-Y.train.list[[ix]])^2))
