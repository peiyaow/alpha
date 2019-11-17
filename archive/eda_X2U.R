library(caret)
library(glmnet)

# bc = BoxCoxTrans(Y)
# Y = predict(bc, Y)
# boxplot(Y[Y<7.98])
# bcb = function(Y, bc) {return((Y*bc$lambda+1)^{1/bc$lambda})}

source('~/Documents/Research/coding/R/alpha/functions.R')
load("~/Documents/Research/coding/R/realdata/ADNI1.RData")
load("~/Documents/Research/coding/R/realdata/myseeds.RData")

n = dim(X)[1]
p = dim(X)[2]

mse.U.vec.list = list()
mse.U.list = list()

# myseeds = floor(1e4*runif(50))
for (i in 1:50){
  set.seed(myseeds[i])
  # set.seed(76)
  ix.train = unlist(createDataPartition(label, times = 1, p = 3/4))
  ix.test = (1:n)[-ix.train]
  ix.list = list(ix.train, ix.test)
  Y.list = lapply(1:2, function(x) Y[ix.list[[x]]]) # train,test
  X.list = lapply(1:2, function(x) X[ix.list[[x]],]) 
  label.list = lapply(1:2, function(x) label[ix.list[[x]]])
  
  X.train = X.list[[1]]
  X.test = X.list[[2]]
  Y.train = Y.list[[1]]
  Y.test = Y.list[[2]]
  label.train = label.list[[1]]
  label.test = label.list[[2]]
  
  # scale the data overall
  X.train.mean = apply(X.train, 2, mean)
  X.train.sd =  apply(X.train, 2, sd)
  X.train = sweep(sweep(X.train, 2, X.train.mean), 2, X.train.sd, "/")
  X.test = sweep(sweep(X.test, 2, X.train.mean), 2, X.train.sd, "/")
  
  # groupwise design matrix
  X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
  X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
  Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
  Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  
  # standardize Y (subtract mean)
  Y.train.mean = lapply(Y.train.list, mean)
  mse.mean.vec = sapply(1:3, function(ix) mean((Y.test.list[[ix]] - Y.train.mean[[ix]])^2))
  Y.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]])
  Y.test.list = lapply(1:3, function(ix) Y.test.list[[ix]] - Y.train.mean[[ix]])
  
  # standardize X (subtract mean)
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  
  # data.frame format
  data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
  data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
  
  # ------------------------------- ALPHA -----------------------------------------
  # par(mfcol = c(2,3))
  X2U.list = lapply(X.train.list, function(X)  X2U(X, plot = T))
  
  # mycut = X2U3(X.train.list)
  # X2U.list = lapply(X.train.list, function(X)  X2U.cut(X, mycut))
  
  H.list = lapply(X2U.list, function(list) list$H)
  F.list = lapply(X2U.list, function(list) list$F_)
  
  # ------------- plot y versus factor --------------
  # par(mfrow=c(3,2))
  # for (ix in 1:3){
  #   plot(F.list[[ix]], Y.train.list[[ix]])
  #   boxplot(F.list[[ix]], horizontal = T)
  # }
  
  # ----------- regress Y on F groupwisely --------------
  # ml.lm.F.class = lapply(1:3, function(ix) lm(Y.train.list[[ix]]~F.list[[ix]]))
  # for (ix in 1:3){
  #   print(summary(ml.lm.F.class[[ix]]))
  # }
  
  U.list = lapply(1:3, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  # PY.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y_.list[[ix]])
  
  # ------------------------------ groupwise regress Y_ on U ------------------------------------------------
  
  # lm (low-rank)
  data.U.list = lapply(1:3, function(ix) data.frame(Y = Y_.list[[ix]], U.list[[ix]]))
  ml.lm.U.class = lapply(1:3, function(ix) lm(Y~., data = data.U.list[[ix]]))
  Yhat.lm.U.class = lapply(1:3, function(ix) predict(ml.lm.U.class[[ix]], newdata = data.test.list[[ix]][,-1]))
  mse.lm.U.class.vec = sapply(1:3, function(ix) mean((Yhat.lm.U.class[[ix]] - Y.test[[ix]])^2))
  mse.lm.U.class = sum(mse.lm.U.class.vec*n.test.vec)/sum(n.test.vec)
  
  # lasso
  ml.lasso.U.class = lapply(1:3, function(ix) cv.glmnet(x=U.list[[ix]], y = Y_.list[[ix]]))
  Yhat.lasso.U.class = lapply(1:3, function(ix) predict(ml.lasso.U.class[[ix]], s=ml.lasso.U.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.lasso.U.class.vec = sapply(1:3, function(ix) mean((Y.test.list[[ix]] - Yhat.lasso.U.class[[ix]])^2))
  mse.lasso.U.class = sum(mse.lasso.U.class.vec*n.test.vec)/sum(n.test.vec)
  
  # ridge
  ml.ridge.U.class = lapply(1:3, function(ix) cv.glmnet(x=U.list[[ix]], y = Y_.list[[ix]], alpha = 0))
  Yhat.ridge.U.class = lapply(1:3, function(ix) predict(ml.ridge.U.class[[ix]], s=ml.ridge.U.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.ridge.U.class.vec = sapply(1:3, function(ix) mean((Y.test.list[[ix]] - Yhat.ridge.U.class[[ix]])^2))
  mse.ridge.U.class = sum(mse.ridge.U.class.vec*n.test.vec)/sum(n.test.vec)
  
  # ------------------------------- global regress Y_ on U ---------------------------------------------------
  U.train = do.call(rbind, U.list)
  Y_.train = do.call(c, Y_.list)
  # PY.train = do.call(c, PY.list)
  
  # lm
  data.U.train = data.frame(Y = Y_.train, U.train)
  ml.lm.U = lm(Y~., data = data.U.train)  
  Yhat.lm.U.test = lapply(1:3, function(ix) predict(ml.lm.U, newdata = data.test.list[[ix]][, -1]))
  mse.lm.U.vec = sapply(1:3, function(ix) mean((Yhat.lm.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.lm.U = sum(mse.lm.U.vec*n.test.vec)/sum(n.test.vec)
  
  # lasso
  ml.lasso.U = cv.glmnet(x=U.train, y=Y_.train)
  Yhat.lasso.U.test = lapply(1:3, function(ix) predict(ml.lasso.U, s=ml.lasso.U$lambda.min, newx = X.test.list[[ix]]))
  mse.lasso.U.vec = sapply(1:3, function(ix) mean((Yhat.lasso.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.lasso.U = sum(mse.lasso.U.vec*n.test.vec)/sum(n.test.vec)
  
  # ridge
  ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, alpha = 0)
  Yhat.ridge.U.test = lapply(1:3, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
  mse.ridge.U.vec = sapply(1:3, function(ix) mean((Yhat.ridge.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.ridge.U = sum(mse.ridge.U.vec*n.test.vec)/sum(n.test.vec)
  
  mse.U.vec = rbind(mse.lm.U.class.vec, mse.lasso.U.class.vec, mse.ridge.U.class.vec, mse.lm.U.vec, mse.lasso.U.vec, mse.ridge.U.vec)
  mse.U = c(mse.lm.U.class, mse.lasso.U.class, mse.ridge.U.class, mse.lm.U, mse.lasso.U, mse.ridge.U)
  names(mse.U) = c("lm.U.class", "lasso.U.class", "ridge.U.class", "lm.U", "lasso.U", "ridge.U")
  mse.U.vec.list[[i]] = mse.U.vec
  mse.U.list[[i]] = mse.U
}

mse.U.result = do.call(rbind, mse.U.list)
mse.U.vec.result = array(as.numeric(unlist(mse.U.vec.list)), dim=c(6, 3, 50))
mse.U.vec.result = aperm(mse.U.vec.result, c(3, 1, 2))

save(mse.U.result, mse.U.vec.result, mse.mean.vec.result, file = "X2U.RData")


#######
mse.lasso.mtx = do.call(rbind, mse.lasso.list)
mse.lm.mtx = do.call(rbind, mse.lm.list)
mse.U.class.mtx = do.call(rbind, mse.U.class.list)
mse.lm.X.class.mtx = do.call(rbind, mse.lm.X.class.list)
mse.lasso.U.class.mtx = do.call(rbind, mse.lasso.U.class.list)
mse.lasso.X.class.mtx = do.call(rbind, mse.lasso.X.class.list)

lasso.result = cbind(apply(mse.lasso.mtx, 2, mean), apply(mse.lasso.mtx, 2, sd)/sqrt(50))
OLS.result = cbind(apply(mse.lm.mtx, 2, mean), apply(mse.lm.mtx, 2, sd)/sqrt(50))
row.names(lasso.result) = c("global", "groupwise", "proposed")
colnames(lasso.result) = c("MSE", "SE")
row.names(OLS.result) = c("global", "groupwise", "proposed")
colnames(OLS.result) = c("MSE", "SE")

OLS.U.class.result = rbind(apply(mse.U.class.mtx, 2, mean), apply(mse.U.class.mtx, 2, sd)/sqrt(50))
OLS.X.class.result  = rbind(apply(mse.lm.X.class.mtx, 2, mean), apply(mse.lm.X.class.mtx, 2, sd)/sqrt(50))
row.names(OLS.U.class.result) = c("MSE", "SE")
colnames(OLS.U.class.result) = c("AD", "MCI", "NC")
row.names(OLS.X.class.result) = c("MSE", "SE")
colnames(OLS.X.class.result) = c("AD", "MCI", "NC")

lasso.U.class.result = rbind(apply(mse.lasso.U.class.mtx, 2, mean), apply(mse.lasso.U.class.mtx, 2, sd))
lasso.X.class.result = rbind(apply(mse.lasso.X.class.mtx, 2, mean), apply(mse.lasso.X.class.mtx, 2, sd))

row.names(lasso.U.class.result) = c("MSE", "SE")
colnames(lasso.U.class.result) = c("AD", "MCI", "NC")
row.names(lasso.X.class.result) = c("MSE", "SE")
colnames(lasso.X.class.result) = c("AD", "MCI", "NC")

library(xtable)
xtable(lasso.result)
xtable(OLS.result)

xtable(lasso.U.class.result)
xtable(lasso.X.class.result)
xtable(OLS.U.class.result)
xtable(OLS.X.class.result)

