library(readr)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)

# t = 1
# X0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/X1a.txt", col_names = FALSE))
# Y0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/Y5T.txt", col_names = FALSE))[, 4+5*(t-1)]
# label0 = as.ordered(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/label1.txt", col_names = FALSE)$X1)
# source('~/Documents/Research/coding/R/alpha/functions.R')
# 
# # remove NAs in Y
# id.cut1 = !is.na(Y0)
# X0 = X0[id.cut1, ]
# Y0 = Y0[id.cut1]
# label0 = label0[id.cut1]
# 
# # remove negatives in Y
# id.cut2 = (Y0 > 0)
# X0 = X0[id.cut2,]
# Y0 = Y0[id.cut2]
# label0 = label0[id.cut2]
# 
# # only select MRI+PET
# X0 = X0[, 1:93]
# X = X0
# Y = Y0
# label = label0

# bc = BoxCoxTrans(Y)
# Y = predict(bc, Y)
# boxplot(Y[Y<7.98])
# bcb = function(Y, bc) {return((Y*bc$lambda+1)^{1/bc$lambda})}
source('~/Documents/Research/coding/R/alpha/functions.R')
load("~/Documents/Research/coding/R/realdata/ADNI1.RData")

n = dim(X)[1]
p = dim(X)[2]

mse.lasso.list = list()
mse.lm.list = list()
mse.U.class.list = list()
mse.lm.X.class.list = list()
mse.lasso.U.class.list = list()
mse.lasso.X.class.list = list()

myseeds = floor(1e4*runif(50))
for (i in 1:50){
  set.seed(myseeds[i])
  set.seed(76)
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
  
  X.train.mean = apply(X.train, 2, mean)
  X.train.sd =  apply(X.train, 2, sd)
  X.train = sweep(sweep(X.train, 2, X.train.mean), 2, X.train.sd, "/")
  X.test = sweep(sweep(X.test, 2, X.train.mean), 2, X.train.sd, "/")
  
  # ------------------------------ overall model -----------------------------------
  # lm
  data.train = data.frame(Y = Y.train, X.train)
  data.test = data.frame(Y = Y.test, X.test)
  ml.lm.X = lm(Y~., data = data.train)
  Yhat.lm.X.test = predict(ml.lm.X, newdata = data.test[, -1])
  mse.lm.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.lm.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.lm.X = mean((Yhat.lm.X.test - Y.test)^2)
  
  # lasso
  ml.lasso.X = cv.glmnet(x=X.train, y=Y.train)
  Yhat.lasso.X.test = predict(ml.lasso.X, s = ml.lasso.X$lambda.min, newx = X.test)
  mse.lasso.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.lasso.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.lasso.X = mean((Yhat.lasso.X.test - Y.test)^2)
  
  # ridge
  ml.ridge.X = cv.glmnet(x=X.train, y=Y.train, alpha = 0)
  Yhat.ridge.X.test = predict(ml.ridge.X, s = ml.ridge.X$lambda.min, newx = X.test)
  mse.ridge.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.ridge.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.ridge.X = mean((Yhat.ridge.X.test - Y.test)^2)
  
  # rf
  ml.rf.X = randomForest(x = X.train, y = Y.train)
  Yhat.rf.X.test = predict(ml.rf.X, newdata = X.test)
  mse.rf.X = mean((Yhat.rf.X.test - Y.test)^2)
  
  # svr
  ml.svr.X = svm(Y~., data.train)
  Yhat.svr.X.test = predict(ml.svr.X, newdata = data.test[, -1])
  mse.svr.X = mean((Yhat.svr.X.test - Y.test)^2)  
    
  # ------------------------------- groupwise model ---------------------------------
  X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
  X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
  Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
  Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  
  # standardize Y
  Y.train.mean = lapply(Y.train.list, mean)
  Y.train.sd = lapply(Y.train.list, sd)
  mse.mean.vec = sapply(1:3, function(ix) mean((Y.test.list[[ix]] - Y.train.mean[[ix]])^2))
  
  # do not divide sd
  Y.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]])
  Y.test.list = lapply(1:3, function(ix) Y.test.list[[ix]] - Y.train.mean[[ix]])
  
  # standardize X
  X.train.mean = lapply(X.train.list, colMeans)
  # X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
  X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  # X.train.list = lapply(1:3, function(ix) sweep(sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
  # X.test.list = lapply(1:3, function(ix) sweep(sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
  
  data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
  data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
  
  # lm
  ml.lm.X.class = lapply(1:3, function(ix) lm(Y~., data=data.train.list[[ix]]))
  Yhat.lm.X.class.test = lapply(1:3, function(ix) predict(ml.lm.X.class[[ix]], newdata = data.test.list[[ix]][, -1]))
  # do not divide sd
  mse.lm.X.class.vec = sapply(1:3, function(ix) mean((Yhat.lm.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.lm.X.class = sum(mse.lm.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # lasso
  ml.lasso.X.class = lapply(1:3, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]]))
  Yhat.lasso.X.class.test = lapply(1:3, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.lasso.X.class.vec = sapply(1:3, function(ix) mean((Yhat.lasso.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # ridge
  ml.ridge.X.class = lapply(1:3, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
  Yhat.ridge.X.class.test = lapply(1:3, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
  mse.ridge.X.class.vec = sapply(1:3, function(ix) mean((Yhat.ridge.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # rf
  ml.rf.X.class = lapply(1:3, function(ix) randomForest(x = X.train.list[[ix]], y = Y.train.list[[ix]]))
  Yhat.rf.X.class.test = lapply(1:3, function(ix) predict(ml.rf.X.class[[ix]], newdata = X.test.list[[ix]]))
  mse.rf.X.class.vec = sapply(1:3, function(ix) mean((Yhat.rf.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.rf.X.class = sum(mse.rf.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # svr
  ml.svr.X.class = lapply(1:3, function(ix) svm(Y~., data=data.train.list[[ix]]))
  Yhat.svr.X.class.test = lapply(1:3, function(ix) predict(ml.svr.X.class[[ix]], newdata = data.test.list[[ix]][, -1]))
  # do not divide sd
  mse.svr.X.class.vec = sapply(1:3, function(ix) mean((Yhat.svr.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.X.class = sum(mse.svr.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  # ------------------------------- ALPHA -----------------------------------------
  X2U.list = lapply(X.train.list, function(X)  X2U(X, plot = T))
  
  mycut = X2U3(X.train.list)
  X2U.list = lapply(X.train.list, function(X)  X2U.cut(X, mycut))
  
  H.list = lapply(X2U.list, function(list) list$H)
  F.list = lapply(X2U.list, function(list) list$F_)
  
  par(mfrow=c(3,2))
  for (ix in 1:3){
    plot(F.list[[ix]], Y.train.list[[ix]])
    boxplot(F.list[[ix]], horizontal = T)
  }
  
  # regress Y on F groupwisely
  ml.lm.F.class = lapply(1:3, function(ix) lm(Y.train.list[[ix]]~F.list[[ix]]))
  for (ix in 1:3){
    print(summary(ml.lm.F.class[[ix]]))
  }
  
  U.list = lapply(1:3, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  PY.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y_.list[[ix]])
  
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
  PY.train = do.call(c, PY.list)
  
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
  
  #WLS
  # n.train.vec = sapply(X.train.list, nrow)
  # w = 1/c(rep(Y.train.sd[[1]], n.train.vec[[1]]), rep(Y.train.sd[[2]], n.train.vec[[2]]), rep(Y.train.sd[[3]], n.train.vec[[3]]))
  # ml.U = lm(Y~., data = data.U.train, weights = w)  
  # Yhat.U.class.test = lapply(1:3, function(ix) predict(ml.U, newdata = data.test.list[[ix]][, -1]))
  # mse.U.class.vec = sapply(1:3, function(ix) mean((Yhat.U.class.test[[ix]]-Y.test.list[[ix]])^2))
  # mse.U.class = sum(mse.U.class.vec*n.test.vec)/sum(n.test.vec)
  # 
  # 
  # ml.lasso.U = cv.glmnet(x=U.train, y=Y_.train, weights = w)
  # Yhat.lasso.U.class.test = lapply(1:3, function(ix) predict(ml.lasso.U, s=ml.lasso.U$lambda.min, newx = X.test.list[[ix]]))
  # mse.lasso.U.class.vec = sapply(1:3, function(ix) mean((Yhat.lasso.U.class.test[[ix]]-Y.test.list[[ix]])^2))
  # mse.lasso.U.class = sum(mse.lasso.U.class.vec*n.test.vec)/sum(n.test.vec)
  
  
  mse.lasso = c(mse.lasso.X, mse.lasso.U, mse.lasso.X.class, mse.lasso.U.class)
  mse.ridge = c(mse.ridge.X, mse.ridge.U, mse.ridge.X.class, mse.ridge.U.class)
  mse.lm = c(mse.lm.X, mse.lm.U, mse.lm.X.class, mse.lm.U.class)
  
  mse.lm.vec = rbind(mse.lm.X.vec, mse.lm.X.class.vec, mse.lm.U.vec, mse.lm.U.class.vec, mse.mean.vec)
  mse.lasso.vec = rbind(mse.lasso.X.vec, mse.lasso.X.class.vec, mse.lasso.U.vec, mse.lasso.U.class.vec, mse.mean.vec)
  mse.ridge.vec = rbind(mse.ridge.X.vec, mse.ridge.X.class.vec, mse.ridge.U.vec, mse.ridge.U.class.vec, mse.mean.vec)
  
  # mse.lasso.list[[i]] = mse.lasso
  # mse.lm.list[[i]] = mse.lm
  # mse.U.class.list[[i]] = mse.U.class.vec
  # mse.lm.X.class.list[[i]] = mse.lm.X.class.vec
  # mse.lasso.U.class.list[[i]] = mse.lasso.U.class.vec
  # mse.lasso.X.class.list[[i]] = mse.lasso.X.class.vec
}

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

