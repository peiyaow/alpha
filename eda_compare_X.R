library(caret)
library(glmnet)
library(randomForest)
library(e1071)

# bc = BoxCoxTrans(Y)
# Y = predict(bc, Y)
# boxplot(Y[Y<7.98])
# bcb = function(Y, bc) {return((Y*bc$lambda+1)^{1/bc$lambda})}

source('~/Documents/Research/coding/R/alpha/functions.R')
load("~/Documents/Research/coding/R/realdata/ADNI1.RData")
load("~/Documents/Research/coding/R/realdata/myseeds.RData")

X = X[-c(167,770),]
label = label[-c(167,770)]
Y = Y[-c(167,770)]

X = X[-151,]
label = label[-151]
Y = Y[-151]

n = dim(X)[1]
p = dim(X)[2]

mse.X.vec.list = list()
mse.X.list = list()

for (i in 1:50){
  set.seed(myseeds[i])
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
  mse.rf.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.rf.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.rf.X = mean((Yhat.rf.X.test - Y.test)^2)
  
  # svr
  ml.svr.radial.X = tune.svm(Y~., data = data.train, cost = 10^(-4:0), kernel = "radial")
  Yhat.svr.radial.X.test = predict(ml.svr.radial.X$best.model, newdata = data.test[, -1])
  mse.svr.radial.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.svr.radial.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.svr.radial.X = mean((Yhat.svr.radial.X.test - Y.test)^2)  
  
  ml.svr.linear.X = tune.svm(Y~., data = data.train, cost = 10^(-4:0), kernel = "linear")
  Yhat.svr.linear.X.test = predict(ml.svr.linear.X$best.model, newdata = data.test[, -1])
  mse.svr.linear.X.vec = sapply(c(4,3,0), function(l) mean((Yhat.svr.linear.X.test[label.test==l] - Y.test[label.test==l])^2))
  mse.svr.linear.X = mean((Yhat.svr.linear.X.test - Y.test)^2) 
  
  # ------------------------------- groupwise model ---------------------------------
  X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
  X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
  Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
  Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
  n.train.vec = sapply(X.train.list, nrow)
  n.test.vec = sapply(X.test.list, nrow)
  
  # standardize Y
  Y.train.mean = lapply(Y.train.list, mean)
  mse.mean.vec = sapply(1:3, function(ix) mean((Y.test.list[[ix]] - Y.train.mean[[ix]])^2))
  mse.mean = sum(mse.mean.vec*n.test.vec)/sum(n.test.vec)
  
  # do not divide sd
  Y.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]])
  Y.test.list = lapply(1:3, function(ix) Y.test.list[[ix]] - Y.train.mean[[ix]])
  
  # standardize X
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  
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
  # ml.svr.X.class = lapply(1:3, function(ix) svm(Y~., data=data.train.list[[ix]]))
  # Yhat.svr.X.class.test = lapply(1:3, function(ix) predict(ml.svr.X.class[[ix]], newdata = data.test.list[[ix]][, -1]))
  # # do not divide sd
  # mse.svr.X.class.vec = sapply(1:3, function(ix) mean((Yhat.svr.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  # mse.svr.X.class = sum(mse.svr.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  ml.svr.radial.X.class = lapply(1:3, function(ix) tune.svm(Y~., data=data.train.list[[ix]], cost = 10^(-4:0), kernel = "radial"))
  Yhat.svr.radial.X.class.test = lapply(1:3, function(ix) predict(ml.svr.radial.X.class[[ix]]$best.model, newdata = data.test.list[[ix]][, -1]))
  # do not divide sd
  mse.svr.radial.X.class.vec = sapply(1:3, function(ix) mean((Yhat.svr.radial.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.radial.X.class = sum(mse.svr.radial.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  ml.svr.linear.X.class = lapply(1:3, function(ix) tune.svm(Y~., data=data.train.list[[ix]], cost = 10^(-4:0), kernel = "linear"))
  Yhat.svr.linear.X.class.test = lapply(1:3, function(ix) predict(ml.svr.linear.X.class[[ix]]$best.model, newdata = data.test.list[[ix]][, -1]))
  # do not divide sd
  mse.svr.linear.X.class.vec = sapply(1:3, function(ix) mean((Yhat.svr.linear.X.class.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.linear.X.class = sum(mse.svr.linear.X.class.vec*n.test.vec)/sum(n.test.vec)
  
  mse.X.vec = rbind(mse.lm.X.class.vec, mse.lasso.X.class.vec, mse.ridge.X.class.vec, mse.rf.X.class.vec, mse.svr.radial.X.class.vec, mse.svr.linear.X.class.vec, mse.lm.X.vec, mse.lasso.X.vec, mse.ridge.X.vec, mse.rf.X.vec, mse.svr.radial.X.vec, mse.svr.linear.X.vec, mse.mean.vec)
  mse.X = c(mse.lm.X.class, mse.lasso.X.class, mse.ridge.X.class, mse.rf.X.class, mse.svr.radial.X.class, mse.svr.linear.X.class, mse.lm.X, mse.lasso.X, mse.ridge.X, mse.rf.X, mse.svr.radial.X, mse.svr.linear.X, mse.mean)
  names(mse.X) = c("lm.X.class", "lasso.X.class", "ridge.X.class", "rf.X.class", "svr.radial.X.class", "svr.linear.X.class", "lm.X", "lasso.X", "ridge.X", "rf.X", "svr.radial.X", "svr.linear.X", "group_mean")
  mse.X.vec.list[[i]] = mse.X.vec
  mse.X.list[[i]] = mse.X
}

mse.X.result = do.call(rbind, mse.X.list)
mse.X.vec.result = array(as.numeric(unlist(mse.X.vec.list)), dim=c(11, 3, 50))
mse.X.vec.result = aperm(mse.X.vec.result, c(3, 1, 2))

save(mse.X.result, mse.X.vec.result, file = "compare_X_0.RData")

