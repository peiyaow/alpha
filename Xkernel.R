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

mse.U.vec.list = list()
mse.U.list = list()

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
  n.train = sum(n.train.vec)
  n.test = sum(n.test.vec)
  
  
  # standardize Y
  Y.train.mean = lapply(Y.train.list, mean)
  Y.train.sd = lapply(Y.train.list, sd)
  
  # standardize Y (subtract mean and divide sd)
  Y.train.list = lapply(1:3, function(ix) (Y.train.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])
  Y.test.list = lapply(1:3, function(ix) (Y.test.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])
  
  # standardize X (subtract mean)
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
  X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))
  
  # data.frame format
  data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
  data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
  
  # ------------------------------- svr linear PCA --------------------------------------
  mycut = X2U4(X.train.list, plot = T)
  X2U.list = lapply(X.train.list, function(X)  X2U.cut(X, mycut))
  H.list = lapply(X2U.list, function(list) list$H)
  F.list = lapply(X2U.list, function(list) list$F_)
  U.list = lapply(1:3, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])

  U.train = do.call(rbind, U.list)
  Y_.train = do.call(c, Y_.list)
  data.U.train = data.frame(Y = Y_.train, U.train)

  # svr
  ml.svr.U = svm(Y~., data.U.train,  kernel = "radial")
  Yhat.svr.U.test = lapply(1:3, function(ix) predict(ml.svr.U, newdata = data.test.list[[ix]]))
  mse.svr.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.svr.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.U = sum(mse.svr.U.vec*n.test.vec)/sum(n.test.vec)

  # ------------------------------- svr gaussian PCA --------------------------------------
  mycut = X2U4.kernel(X.train.list, plot = T)
  X2U.list = lapply(X.train.list, function(X)  X2U.kernel.cut(X, mycut))
  
  H.list = lapply(X2U.list, function(list) list$H)
  F.list = lapply(X2U.list, function(list) list$F_)
  Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
  Y_.train = do.call(c, Y_.list)
  # ---------------------------- directly multiply X by H ---------------------------------
  U.list = lapply(1:3, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
  U.train = do.call(rbind, U.list)
  data.U.train = data.frame(Y = Y_.train, U.train)
  ml.svr.HX.U = svm(Y~., data.U.train,  kernel = "radial")
  Yhat.svr.HX.U.test = lapply(1:3, function(ix) predict(ml.svr.HX.U, newdata = data.test.list[[ix]]))
  mse.svr.HX.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.svr.HX.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.HX.U = sum(mse.svr.HX.U.vec*n.test.vec)/sum(n.test.vec)
  
  # ------------------------------- gaussian PCA + Kernel trick --------------------------------------
  gamma = 1/p
  rbf = rbfdot(sigma = gamma)
  
  # demeaned kernel matrix for training
  K.mat = matrix(0, nrow = n.train, ncol = n.train)
  ix.vec = c(0,cumsum(n.train.vec))
  for (i in 1:2){
    for (j in (i+1):3){
      K.mat[(ix.vec[i]+1):ix.vec[i+1], (ix.vec[j]+1):ix.vec[j+1]] = X2K(X.train.list[[i]], X.train.list[[j]])
    }
  }
  K.mat = K.mat + t(K.mat)
  for (i in 1:3){
    K.mat[(ix.vec[i]+1):ix.vec[i+1], (ix.vec[i]+1):ix.vec[i+1]] = X2K(X.train.list[[i]], X.train.list[[i]])
  }
  kermat = kernelMatrix(kernel = rbf, x = X.train)
  H_diag = as.matrix(bdiag(H.list))
  kermat@.Data = H_diag%*%K.mat%*%H_diag
  
  ml.ksvm = ksvm(x = kermat, y = Y_.train)
  
  # demeaned kernel matrix for testing
  K.mat.test = matrix(0, nrow = n.test, ncol = n.train)
  ix.test.vec = c(0,cumsum(n.test.vec))
  ix.train.vec = c(0,cumsum(n.train.vec))
  for (i in 1:3){
    for (j in 1:3){
      K.mat.test[(ix.test.vec[i]+1):ix.test.vec[i+1], (ix.train.vec[j]+1):ix.train.vec[j+1]] = X2K.test(X.test.list[[i]], X.train.list[[j]])
    }
  }
  kermat.test = kernelMatrix(kernel = rbf, x = X.test, y = X.train)
  kermat.test@.Data = K.mat.test
  kermat.test@.Data = kermat.test[,ml.ksvm@alphaindex]
  Yhat.svr.ker.U.test = lapply(c(4,3,0), function(ix) predict(ml.ksvm, kermat.test)[label.test==ix])
  mse.svr.ker.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.svr.ker.U.test[[ix]]-Y.test.list[[ix]])^2))
  mse.svr.ker.U = sum(mse.svr.ker.U.vec*n.test.vec)/sum(n.test.vec)
  
  
  mse.U.vec = rbind(mse.svr.U.vec, mse.svr.ker.U.vec)
  mse.U = c(mse.svr.U, mse.svr.ker.U)
  
  mse.U.vec.list[[i]] = mse.U.vec
  mse.U.list[[i]] = mse.U
}

mse.U.result = do.call(rbind, mse.U.list)
mse.U.vec.result = array(as.numeric(unlist(mse.U.vec.list)), dim=c(2, 3, 50))
mse.U.vec.result = aperm(mse.U.vec.result, c(3, 1, 2))
