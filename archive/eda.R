library(readr)
library(caret)
library(glmnet)
library(corpcor)


X2U = function(X, K = 93){
  n = nrow(X)
  p = ncol(X)
  PCA.res = eigen(X%*%t(X)/n/p)
  
  eigen_vals = PCA.res$values
  eigen_vals = c(eigen_vals, Inf)
  K = min(which.max(eigen_vals[1:p]/eigen_vals[2:(p+1)]), K)
  
  F_ = PCA.res$vectors[,1:K]*sqrt(n)
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  U = H%*%X
  return(list(U, H, P, K))
}

t = 1
# mac: /Users/MonicaW/Documents/GitHub/LWPR
X0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/X1a.txt", col_names = FALSE))
Y0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/Y5T.txt", col_names = FALSE))[, 4+5*(t-1)]
label0 = as.ordered(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/label1.txt", col_names = FALSE)$X1)

# remove NAs in Y
id.cut1 = !is.na(Y0)
X0 = X0[id.cut1, ]
Y0 = Y0[id.cut1]
label0 = label0[id.cut1]

# remove negatives in Y
id.cut2 = (Y0 > 0)
X0 = X0[id.cut2,]
Y0 = Y0[id.cut2]
label0 = label0[id.cut2]

# only select MRI+PET
X0 = X0[, 1:93]

X = X0
Y = Y0
#Y = log(Y+1)
label = label0
n = dim(X)[1]
p = dim(X)[2]

result.list = list()
myseeds = floor(1e4*runif(50))
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
  
  X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
  X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
  Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
  Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
  
  X.train.mean = lapply(X.train.list, colMeans)
  X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
  X.train.list = lapply(1:3, function(ix) sweep(sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
  X.test.list = lapply(1:3, function(ix) sweep(sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
  n.test.vec = sapply(X.test.list, nrow)
  
  # train lasso within each class
  ml.X.class = lapply(1:3, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], standardize = F))
  Y.hat.train.class.list = lapply(1:3, function(ix) predict(ml.X.class[[ix]], s = ml.X.class[[ix]]$lambda.min, newx = X.train.list[[ix]]))

  mse.X.class = sapply(1:3, function(ix) mean((predict(ml.X.class[[ix]], s = ml.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]) - Y.test.list[[ix]])^2))
  mse.X.class = sum(mse.X.class*n.test.vec)/sum(n.test.vec)
  
  # overall model
  ml.X = cv.glmnet(x=X.train, y=Y.train, standardize = F)
  Y.hat.train.X.list = lapply(c(4,3,0), function(l) predict(ml.X, s = ml.X$lambda.min, newx = X.train)[label.train == l])
  mse.train.X.list = lapply(1:3, function(ix) mean((Y.hat.train.X.list[[ix]] - Y.train.list[[ix]])^2))
  mse.X = mean((predict(ml.X, s = ml.X$lambda.min, newx = X.test) - Y.test)^2)
  
  mse.list = list()
  for (K in seq(10,90,10)){
    U.train.list = lapply(X.train.list, function(X) X2U(X, 50)[[1]])
    H.train.list = lapply(X.train.list, function(X) X2U(X, 50)[[2]])
    P.train.list = lapply(X.train.list, function(X) X2U(X, 50)[[3]])
    K.train.list = lapply(X.train.list, function(X) X2U(X, 50)[[4]])
    Y_.train.list = lapply(1:3, function(ix) H.train.list[[ix]]%*%Y.train.list[[ix]])
    #PY.train.list = lapply(1:3, function(ix) P.train.list[[ix]]%*%Y.train.list[[ix]])
    #PY.train = do.call(c, PY.train.list)
    
    U_.train.list = lapply(1:3, function(ix) H.train.list[[ix]]%*%U.train.list[[ix]])
    Y_.train = do.call(c, Y_.train.list)
    U_.train = do.call(rbind, U_.train.list)
    ml.U_ = cv.glmnet(x=U_.train, y=Y_.train, standardize = F)
    
    
    coef(ml.U_, s=ml.U_$lambda.min)
    
    
    
    
    
    
    
    
    
    Y_.hat.train.U_.list = lapply(c(4,3,0), function(l) predict(ml.U_, s = ml.U_$lambda.min, newx = U_.train)[label.train == l])
    mse.Y_.train.U_.list = lapply(1:3, function(ix) mean((Y_.hat.train.U_.list[[ix]] - Y_.train.list[[ix]])^2))
    Y_.hat.train.class.list = lapply(1:3, function(ix) H.train.list[[ix]]%*%Y.hat.train.class.list[[ix]])
    Y_.hat.train.X.list = lapply(1:3, function(ix) H.train.list[[ix]]%*%Y.hat.train.X.list[[ix]])
    mse.Y_.train.class.list = lapply(1:3, function(ix) mean((Y_.hat.train.class.list[[ix]] - Y_.train.list[[ix]])^2))
    mse.Y_.train.X.list = lapply(1:3, function(ix) mean((Y_.hat.train.X.list[[ix]] - Y_.train.list[[ix]])^2))
    
    Y.hat.train.U_.list = lapply(1:3, function(ix) pseudoinverse(H.train.list[[ix]])%*%Y_.hat.train.U_.list[[ix]])
    mse.train.U_.list = lapply(1:3, function(ix) mean((Y.hat.train.U_.list[[ix]] - Y.train.list[[ix]])^2))
    
    
    
    U.test.list = lapply(X.test.list, function(X) X%*%X2U(X, K)[[2]])
    F.train.list = lapply(1:3, function(ix) X.train.list[[ix]] - U.train.list[[ix]])
    F.test.list = lapply(1:3, function(ix) X.test.list[[ix]] - U.test.list[[ix]])
    
    ml.F.class = lapply(1:3, function(ix) cv.glmnet(x=F.train.list[[ix]], y=Y.train.list[[ix]], standardize = T))
    mse.F.class = sapply(1:3, function(ix) mean((predict(ml.F.class[[ix]], s = ml.F.class[[ix]]$lambda.min, newx = F.test.list[[ix]]) - Y.test.list[[ix]])^2))
    
    res.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - predict(ml.F.class[[ix]], s = ml.F.class[[ix]]$lambda.min, newx = F.train.list[[ix]]))
    res.train = do.call(c, res.train.list)
    plot(res.train)
    U.train = do.call(rbind, U.train.list)
    U.test = do.call(rbind, U.test.list)
    ml.U = cv.glmnet(x=U.train, y=res.train, standardize = T)
    res.hat.test = predict(ml.U, s = ml.U$lambda.min, newx = U.test)
    Y.hat.test = do.call(c, lapply(1:3, function(ix) predict(ml.F.class[[ix]], s = ml.F.class[[ix]]$lambda.min, newx = F.test.list[[ix]]))) + res.hat.test
    mse.my = mean((do.call(c, Y.test.list) - Y.hat.test)^2)
    
    
    F.train = X.train - U.train
    F.test = X.test - U.test
    ml.F = cv.glmnet(x=F.train, y=Y.train, standardize = T)
    mse.F = mean((exp(predict(ml.F, s = ml.F$lambda.min, newx = F.test)) - exp(Y.test))^2)
    mse.list[[as.character(K)]] = mse.F
  }
  result.list[[i]] = c(mse.X.class, mse.X, do.call(c, mse.list))
}

result = do.call(rbind, result.list)
colnames(result)[1:2] = c("class", "global")

PCA.min = apply(result[,3:11], 1, min)

rbind(apply(cbind(result, PCA.min), 2, mean), apply(cbind(result, PCA.min), 2, sd))




