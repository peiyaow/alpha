# input:
# Y.train.list, X.train.list, X.test.list, threshold and original X.train.list and X.test.list
# output: 
# Yhat.test.list
lm.U.threshold.method2 = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list1 = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list1[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list1[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list1[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  # new prediction
  X.combine.list = lapply(1:n_label, function(ix) rbind(X.train.list[[ix]], X.test.list[[ix]]))
  X.combine.mean = lapply(X.combine.list, colMeans)
  X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
  
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.combine.list[[ix]], K = K.list[[ix]], plot = F))
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
  F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  F.train.new.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][1:n.train.vec[ix],]))
  U.train.new.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
  
  for (ix in 1:n_label){
    F.list[[ix]][,-1] = t(t(F.list[[ix]][,-1])*sign(diag(as.matrix(cor(F.train.new.list[[ix]][,-1], F.train.list[[ix]][,-1])))))
  }
  U.test.list = lapply(1:n_label, function(ix) U.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],])
  F.test.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],]))
  
  ml.lm.U = lm(Y~., data = data.U.train)  
  ix.vec = c(0,cumsum(n.train.vec))
  sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
  w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
  
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  
  return(list(Yhat.test.list = Yhat.test.list))
}

ridge.U.threshold.method2 = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list1 = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list1[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list1[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list1[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  # new prediction
  X.combine.list = lapply(1:n_label, function(ix) rbind(X.train.list[[ix]], X.test.list[[ix]]))
  X.combine.mean = lapply(X.combine.list, colMeans)
  X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
  
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.combine.list[[ix]], K = K.list[[ix]], plot = F))
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
  F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  F.train.new.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][1:n.train.vec[ix],]))
  U.train.new.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
  
  for (ix in 1:n_label){
    F.list[[ix]][,-1] = t(t(F.list[[ix]][,-1])*sign(diag(as.matrix(cor(F.train.new.list[[ix]][,-1], F.train.list[[ix]][,-1])))))
  }
  U.test.list = lapply(1:n_label, function(ix) U.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],])
  F.test.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],]))
  
  ml.lm.U = lm(Y~., data = data.U.train)  
  ix.vec = c(0,cumsum(n.train.vec))
  sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
  w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
  
  # weighted ridge: weights from OLS.U
  ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
  HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  return(list(Yhat.test.list = Yhat.test.list))
}

lm.U2.threshold.method2 = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list1 = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list1[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list1[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list1[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  # new prediction
  X.combine.list = lapply(1:n_label, function(ix) rbind(X.train.list[[ix]], X.test.list[[ix]]))
  X.combine.mean = lapply(X.combine.list, colMeans)
  X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
  
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.combine.list[[ix]], K = K.list[[ix]], plot = F))
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
  F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  F.train.new.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][1:n.train.vec[ix],]))
  U.train.new.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
  
  for (ix in 1:n_label){
    F.list[[ix]][,-1] = t(t(F.list[[ix]][,-1])*sign(diag(as.matrix(cor(F.train.new.list[[ix]][,-1], F.train.list[[ix]][,-1])))))
  }
  U.test.list = lapply(1:n_label, function(ix) U.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],])
  F.test.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],]))
  
  
  # compute weights
  ml.lm.U = lm(Y~., data = data.U.train)  
  HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
  sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
  w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
  
  # WLS.U: weights from OLS.U and OLS.F
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
  res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) cbind(1, U.test.list[[ix]])%*%ml.lm.U$coefficients)
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  
  return(list(Yhat.test.list = Yhat.test.list))
}

ridge.U2.threshold.method2 = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list1 = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list1[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list1[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list1[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  # new prediction
  X.combine.list = lapply(1:n_label, function(ix) rbind(X.train.list[[ix]], X.test.list[[ix]]))
  X.combine.mean = lapply(X.combine.list, colMeans)
  X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
  
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.combine.list[[ix]], K = K.list[[ix]], plot = F))
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
  F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  F.train.new.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][1:n.train.vec[ix],]))
  U.train.new.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
  
  for (ix in 1:n_label){
    F.list[[ix]][,-1] = t(t(F.list[[ix]][,-1])*sign(diag(as.matrix(cor(F.train.new.list[[ix]][,-1], F.train.list[[ix]][,-1])))))
  }
  U.test.list = lapply(1:n_label, function(ix) U.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],])
  F.test.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],]))
  
  # compute weights
  ml.lm.U = lm(Y~., data = data.U.train)  
  HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
  sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
  w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
  
  # weighted ridge: weights from OLS.U
  ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
  HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  return(list(Yhat.test.list = Yhat.test.list))
}




threshold = 0.1
hehe = lm.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)


cv.select_threshold.method2 = function(threshold.vec, X.train.list, Y.train.list, nfolds){
  n_label = length(Y.train.list)
  flds.list = lapply(1:n_label, function(l) createFolds(Y.train.list[[l]], k = nfolds, list = TRUE, returnTrain = FALSE))
  n_thres = length(threshold.vec)
  mse.list = list()
  for (k in 1:nfolds){
    mse.list[[k]] = list()
    for (t in 1:n_thres){
      X.train.list1 = lapply(1:n_label, function(l) X.train.list[[l]][unlist(flds.list[[l]][-k]),])
      X.train.list2 = lapply(1:n_label, function(l) X.train.list[[l]][unlist(flds.list[[l]][k]),])
      Y.train.list1 = lapply(1:n_label, function(l) Y.train.list[[l]][unlist(flds.list[[l]][-k])])
      Y.train.list2 = lapply(1:n_label, function(l) Y.train.list[[l]][unlist(flds.list[[l]][k])])
      
      ml.res = lm.U.threshold.method2(X.train.list1, Y.train.list1, X.train.list2, threshold.vec[t])
      
      n.vec2 = sapply(1:n_label, function(l) length(Y.train.list2[[l]]))
      mse.vec = sapply(1:n_label, function(l)  mean((ml.res$Yhat.test.list[[l]] - Y.train.list2[[l]])^2))
      mse.list[[k]][[t]] = sum(mse.vec*n.vec2)/sum(n.vec2)
    }
    mse.list[[k]] = do.call(c, mse.list[[k]])
  }
  mse.vec = apply(do.call(rbind, mse.list), 2, mean)
  threshold = threshold.vec[which.min(mse.vec)]
  return(list(threshold = threshold, mse.vec = mse.vec))
}

haha = cv.select_threshold.method2(threshold.vec, X.train.list, Y.train.list, nfolds = 10)


Yhat.test.list = lm.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.lm.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U1 = sum(mse.lm.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = lm.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.lm.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U2 = sum(mse.lm.U.vec2*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.ridge.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U1 = sum(mse.ridge.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.ridge.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U2 = sum(mse.ridge.U.vec2*n.test.vec)/sum(n.test.vec)





lm.U.threshold.method2.single = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)

  # F.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Fx)))
  F.test.list = lapply(1:n_label, function(ix) matrix(t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Fx)), nrow = n.test.vec[ix]))
  U.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Ux)))
  
  ml.lm.U = lm(Y~., data = data.U.train)  
  ix.vec = c(0,cumsum(n.train.vec))
  sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
  w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
  
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  
  return(list(Yhat.test.list = Yhat.test.list))
}

lm.U2.threshold.method2.single = function(X.train.list, Y.train.list, X.test.list, threshold){
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  n.vec = n.train.vec + n.test.vec
  n_label = length(X.train.list)
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  
  K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list[[l]], threshold)$K)
  X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = K.list[[ix]], plot = F))
  
  H.list = lapply(X2U.list, function(list) list$H)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  F.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Fx)))
  U.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Ux)))
  
  # compute weights
  ml.lm.U = lm(Y~., data = data.U.train)  
  HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
  sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
  w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
  
  # WLS.U: weights from OLS.U and OLS.F
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
  res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
  HYhat.test.list = lapply(1:n_label, function(ix) cbind(1, U.test.list[[ix]])%*%ml.lm.U$coefficients)
  Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
  
  return(list(Yhat.test.list = Yhat.test.list))
}

cv.select_threshold.method2.single = function(threshold.vec, X.train.list, Y.train.list, nfolds){
  n_label = length(Y.train.list)
  flds.list = lapply(1:n_label, function(l) createFolds(Y.train.list[[l]], k = nfolds, list = TRUE, returnTrain = FALSE))
  n_thres = length(threshold.vec)
  mse.list = list()
  for (k in 1:nfolds){
    mse.list[[k]] = list()
    for (t in 1:n_thres){
      X.train.list1 = lapply(1:n_label, function(l) X.train.list[[l]][unlist(flds.list[[l]][-k]),])
      X.train.list2 = lapply(1:n_label, function(l) X.train.list[[l]][unlist(flds.list[[l]][k]),])
      Y.train.list1 = lapply(1:n_label, function(l) Y.train.list[[l]][unlist(flds.list[[l]][-k])])
      Y.train.list2 = lapply(1:n_label, function(l) Y.train.list[[l]][unlist(flds.list[[l]][k])])
      
      ml.res = lm.U.threshold.method2.single(X.train.list1, Y.train.list1, X.train.list2, threshold.vec[t])
      
      n.vec2 = sapply(1:n_label, function(l) length(Y.train.list2[[l]]))
      mse.vec = sapply(1:n_label, function(l)  mean((ml.res$Yhat.test.list[[l]] - Y.train.list2[[l]])^2))
      mse.list[[k]][[t]] = sum(mse.vec*n.vec2)/sum(n.vec2)
    }
    mse.list[[k]] = do.call(c, mse.list[[k]])
  }
  mse.vec = apply(do.call(rbind, mse.list), 2, mean)
  threshold = threshold.vec[which.min(mse.vec)]
  return(list(threshold = threshold, mse.vec = mse.vec))
}

hehe = cv.select_threshold.method2.single(threshold.vec, X.train.list, Y.train.list, nfolds = 10)
haha = lm.U.threshold.method2.single(X.train.list, Y.train.list, X.test.list, 0.2)







