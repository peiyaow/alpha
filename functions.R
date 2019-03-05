library(kernlab)
library(CVST)
X2U = function(X, K=50, plot = T){
  n = nrow(X)
#  p = ncol(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  if (plot){
    par(mfrow=c(3,1))
    plot(eigen_vals[1:K])
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
    plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  }
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = PCA.res$vectors[,1:K]*sqrt(n)
  #F1 = PCA.res$vectors[,1]*sqrt(n)
  #F2 = PCA.res$vectors[,2]*sqrt(n)
  F2 = PCA.res$vectors[,1:2]*sqrt(n) # first 2 factors
  F_res = PCA.res$vectors[,3:n]*sqrt(n) # other factors
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  L2 = solve(t(F2)%*%F2)%*%t(F2)%*%X
  L_res = solve(t(F_res)%*%F_res)%*%t(F_res)%*%X
  #return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1, L = L, L2 = L2))
  return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, L = L, L2 = L2, L_res = L_res))
}

X2U1 = function(X, K=50, plot = F){
  n = nrow(X)
  #  p = ncol(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  if (plot){
    par(mfrow=c(3,1))
    plot(eigen_vals[1:K])
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
    plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  }
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = cbind(1,PCA.res$vectors[,1:K]*sqrt(n))
  F1 = F_[,1]
  #F1 = PCA.res$vectors[,1]*sqrt(n)
  #F2 = PCA.res$vectors[,2]*sqrt(n)
  F2 = PCA.res$vectors[,1:2]*sqrt(n) # first 2 factors
  F_res = PCA.res$vectors[,3:n]*sqrt(n) # other factors
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  P1 = F1%*%solve(t(F1)%*%F1)%*%t(F1)
  H = diag(n) - P
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  L2 = solve(t(F2)%*%F2)%*%t(F2)%*%X
  L_res = solve(t(F_res)%*%F_res)%*%t(F_res)%*%X
  #return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1, L = L, L2 = L2))
  return(list(H = H, P = P, P1 = P1, K = K, F_ = F_, F2 = F2, L = L, L2 = L2, L_res = L_res))
}

X2U2 = function(X, K=NULL, plot = F){ # estimate the factors give K
  n = nrow(X)
  p = floor(ncol(X)/2)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  if (plot){
    par(mfrow=c(3,1))
    plot(eigen_vals[1:(p-1)])
    plot(eigen_vals[1:(p-1)]/eigen_vals[2:(p)])
    plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  }
  if (is.null(K)){
    K = which.max(eigen_vals[1:(p-1)]/eigen_vals[2:(p)])
  }
  F_ = cbind(1,PCA.res$vectors[,1:K]*sqrt(n))
  F1 = F_[,1]
  #F1 = PCA.res$vectors[,1]*sqrt(n)
  #F2 = PCA.res$vectors[,2]*sqrt(n)
  F2 = PCA.res$vectors[,1:2]*sqrt(n) # first 2 factors
  F_res = PCA.res$vectors[,3:n]*sqrt(n) # other factors
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  P1 = F1%*%solve(t(F1)%*%F1)%*%t(F1)
  H = diag(n) - P
  U = H%*%X
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  L2 = solve(t(F2)%*%F2)%*%t(F2)%*%X
  L_res = solve(t(F_res)%*%F_res)%*%t(F_res)%*%X
  #return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1, L = L, L2 = L2))
  return(list(H = H, P = P, P1 = P1, K = K, F_ = F_, F2 = F2, L = L, L2 = L2, L_res = L_res, U = U))
}


FnU = function(X0, x, F0, K){
  # give single x, get F and U from x
  X = rbind(X0,x)
  n = nrow(X)
  X.mean = apply(X, 2, mean)
  X = sweep(X, 2, X.mean)
  FF = X2U2(X, K = K)$F_
  U = X2U2(X, K = K)$U
  F0.new = FF[1:(n-1), ]
  Fx = FF[n,]
  Fx[-1] = t(t(Fx[-1])*sign(diag(as.matrix(cor(F0[,-1], F0.new[,-1])))))
  Ux = U[n,]
  return(list(F0.new = F0.new, Fx = Fx, Ux = Ux))
}

predict_method1_single = function(X.train.list, Y.train.list, x, l){
  n.train.vec = sapply(X.train.list, nrow)
  n_label = length(X.train.list)
  X.combine.list = X.train.list
  X.combine.list[[l]] = rbind(X.combine.list[[l]],x)
  X.combine.mean = lapply(X.combine.list, colMeans)
  X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
  
  X2U.list = lapply(X.combine.list, function(X)  X2U2(X))
  H.list = lapply(X2U.list, function(list) list$H)
  P1.list = lapply(X2U.list, function(list) list$P1)
  K.list = lapply(X2U.list, function(list) list$K)
  P.list = lapply(X2U.list, function(list) list$P)
  
  U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
  F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
  
  HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]][1:n.train.vec[ix],1:n.train.vec[ix]]%*%Y.train.list[[ix]])
  PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]][1:n.train.vec[ix],1:n.train.vec[ix]]%*%Y.train.list[[ix]])
  U.train.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
  F.train.list = lapply(1:n_label, function(ix) F.list[[ix]][1:n.train.vec[ix],])
  
  u.test = U.list[[l]][n.train.vec[l]+1,]
  f.test = F.list[[l]][n.train.vec[l]+1,]
  
  U.train = do.call(rbind, U.train.list)
  HY.train = do.call(c, HY.train.list)
  data.U.train = data.frame(Y = HY.train, U.train)
  
  #---------------- get weights from OLS.U ----------------
  ml.lm.U = lm(Y~., data = data.U.train)  
  ix.vec = c(0,cumsum(n.train.vec))
  sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
  w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
  
  #---------------- OLS.U with weights --------------------
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
  
  # OLS.F on resid.OLS.U
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.OLS.U.test = f.test%*%ml.lm.F.list[[l]]$coefficients
  HYhat.OLS.U.test = c(1,u.test)%*%ml.lm.U$coefficients
  Yhat.OLS.U.test = PYhat.OLS.U.test + HYhat.OLS.U.test
  
  #---------------- ridge.U with weights ------------------
  ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
  HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
  res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
  
  # OLS.F on resid.ridge.U
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.ridge.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.ridge.U.test = f.test%*%ml.lm.F.list[[l]]$coefficients
  HYhat.ridge.U.test = predict(ml.ridge.U, s=ml.ridge.U$lambda.min, t(u.test))
  Yhat.ridge.U.test = PYhat.ridge.U.test + HYhat.ridge.U.test
  
  #-------------- get weights from OLS.U and OLS.F -----------------
  ml.lm.U = lm(Y~., data = data.U.train)  
  HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
  sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
  w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
  
  #---------------- OLS.U with weights --------------------
  ml.lm.U = lm(Y~., data = data.U.train, weights = w)
  HYhat.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
  res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
  
  # OLS.F on resid.OLS.U
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.OLS.U.test = f.test%*%ml.lm.F.list[[l]]$coefficients
  HYhat.OLS.U.test = c(1,u.test)%*%ml.lm.U$coefficients
  Yhat.OLS.U.test1 = PYhat.OLS.U.test + HYhat.OLS.U.test
  
  #---------------- ridge.U with weights ------------------
  ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
  HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
  res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
  
  # OLS.F on resid.ridge.U
  data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.ridge.train.list[[ix]], F.train.list[[ix]][,-1]))
  ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
  
  PYhat.ridge.U.test = f.test%*%ml.lm.F.list[[l]]$coefficients
  HYhat.ridge.U.test = predict(ml.ridge.U, s=ml.ridge.U$lambda.min, t(u.test))
  Yhat.ridge.U.test1 = PYhat.ridge.U.test + HYhat.ridge.U.test
  
  
  return(list(Yhat = c(Yhat.OLS.U.test, Yhat.ridge.U.test, Yhat.OLS.U.test1, Yhat.ridge.U.test1)))
}

X2U.K = function(X, K = 50, plot = T){
  n = nrow(X)
  #  p = ncol(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  F_ = PCA.res$vectors[,1:K]*sqrt(n)
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  return(list(H = H, P = P, K = K, F_ = F_))
}

X2U3 = function(X.list, K = 50, plot = T){
  eigen.vec = c()
  for (X in X.train.list){
    n = nrow(X)
    PCA.res = eigen(X%*%t(X)/n)
    eigen_vals = PCA.res$values
    eigen.vec = c(eigen.vec, eigen_vals)
  }
  sort_eigen.vec = sort(eigen.vec, decreasing = T)
  if (plot){
    par(mfrow=c(2,1))
    plot(sort_eigen.vec[1:K])
    plot(sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)])
  }
  K = which.max(sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)])
  cut.eigen = sort_eigen.vec[K]
  
  return(cut.eigen)
}

X2U4 = function(X.list, K = 50, plot = T){
  eigen.vec = c()
  for (X in X.list){
    n = nrow(X)
    PCA.res = eigen(X%*%t(X)/n)
    eigen_vals = PCA.res$values
    eigen.vec = c(eigen.vec, eigen_vals)
  }
  sort_eigen.vec = sort(eigen.vec, decreasing = T)
  if (plot){
    par(mfrow=c(2,1))
    plot(sort_eigen.vec[1:K])
    plot(sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)])
    
  }
  ratio.eigen = sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)]
  K = 1
  while (ratio.eigen[K+1] > ratio.eigen[K]){
    K = K+1
  }
  cut.eigen = sort_eigen.vec[K]
  
  return(cut.eigen)
}

X2U4.kernel = function(X.list, K = 50, gamma = NULL, plot = T){
  eigen.vec = c()
  for (X in X.list){
    n = nrow(X)
    p = ncol(X)
    if(is.null(gamma)){
      gamma = 1/p
    }
    rbf = rbfdot(sigma = gamma)
    kermat = kernelMatrix(kernel = rbf, x = X)
    I_n = matrix(1/n, nrow = n, ncol = n)
    kermat = kermat - I_n%*%kermat - kermat%*%I_n + I_n%*%kermat%*%I_n
    kermat = kermat/n
    PCA.res = eigen(kermat)
    eigen_vals = PCA.res$values
    
    eigen.vec = c(eigen.vec, eigen_vals)
  }
  sort_eigen.vec = sort(eigen.vec, decreasing = T)
  if (plot){
    par(mfrow=c(2,1))
    plot(sort_eigen.vec[1:K])
    plot(sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)])
  }
  ratio.eigen = sort_eigen.vec[1:K]/sort_eigen.vec[2:(K+1)]
  K = 1
  while (ratio.eigen[K+1] > ratio.eigen[K]){
    K = K+1
  }
  cut.eigen = sort_eigen.vec[K]
  #print(K)
  
  return(cut.eigen)
}

# X2U.cut = function(X, cut){
#   n = nrow(X)
#   #  p = ncol(X)
#   PCA.res = eigen(X%*%t(X)/n)
#   eigen_vals = PCA.res$values
#   K = n - length(eigen_vals[eigen_vals<cut])
#   if (K == 0){
#     F_ = 0
#     P = matrix(0, nrow = n, ncol = n)
#     H = diag(n) - P
#   }else{
#     F_ = PCA.res$vectors[,1:K]*sqrt(n)
#     P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
#     H = diag(n) - P
#   }
#   return(list(H = H, P = P, K = K, F_ = F_))
# }

X2U.cut = function(X, cut){
  n = nrow(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  K = n - length(eigen_vals[eigen_vals<cut])
  if (K == 0){
    F_ = as.matrix(rep(1, n))
    P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
    H = diag(n) - P
  }else{
    F_ = cbind(1,PCA.res$vectors[,1:K]*sqrt(n))
    P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
    H = diag(n) - P
  }
  F1 = F_[,1]
  P1 = F1%*%solve(t(F1)%*%F1)%*%t(F1)
  return(list(H = H, P = P, K = K, F_ = F_, P1 = P1))
}

X2U.kernel.cut = function(X, cut, gamma = NULL){
  n = nrow(X)
  p = ncol(X)
  if(is.null(gamma)){
    gamma = 1/p
  }
  rbf = rbfdot(sigma = gamma)
  kermat = kernelMatrix(kernel = rbf, x = X)
  I_n = matrix(1/n, nrow = n, ncol = n)
  kermat = kermat - I_n%*%kermat - kermat%*%I_n + I_n%*%kermat%*%I_n
  kermat = kermat/n
  PCA.res = eigen(kermat)
  eigen_vals = PCA.res$values
  print(eigen_vals)
  K = n - length(eigen_vals[eigen_vals<cut])
  print(K)
  if (K == 0){
    F_ = 0
    P = matrix(0, nrow = n, ncol = n)
    H = diag(n) - P
  }else{
    F_ = PCA.res$vectors[,1:K]*sqrt(n)
    P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
    H = diag(n) - P
  }
  return(list(H = H, P = P, K = K, F_ = F_))
}

  
# X2U2 = function(X, K=50, plot = T){
#   n = nrow(X)
#   PCA.res = eigen(X%*%t(X))
#   
#   eigen_vals = PCA.res$values
#   if (plot){
#     par(mfrow=c(2,1))
#     plot(eigen_vals[1:K])
#     plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
#   }
#   K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
#   F_ = PCA.res$vectors[,1:K]*sqrt(n)
#   P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
#   H = diag(n) - P
#   return(list(H = H, P = P, K = K))
# }


# for training data X1 and X2
# input X1, X2
# output standardized kernel K(X1, X2) n1 by n2
X2K = function(X1, X2, gamma = NULL){
  n1 = nrow(X1)
  n2 = nrow(X2)
  p = ncol(X1) # X2 should be the same dim as X1
  if(is.null(gamma)){
    gamma = 1/p
  }
  rbf = rbfdot(sigma = gamma)
  kermat = kernelMatrix(kernel = rbf, x = X1, y = X2)
  I_n1 = matrix(1/n1, nrow = n1, ncol = n1)
  I_n2 = matrix(1/n2, nrow = n2, ncol = n2)
  
  kermat = kermat - I_n1%*%kermat - kermat%*%I_n2 + I_n1%*%kermat%*%I_n2
  return(kermat)
}

X2K.test = function(X.test, X, gamma = NULL){
  n1 = nrow(X)
  n2 = nrow(X.test)
  p = ncol(X) # X2 should be the same dim as X1
  if(is.null(gamma)){
    gamma = 1/p
  }
  rbf = rbfdot(sigma = gamma)
  kermat = kernelMatrix(kernel = rbf, x = X)
  kermat.U.test = kernelMatrix(kernel = rbf, x = X.test, y = X) # n2 by n1
  I_n1 = matrix(1/n1, nrow = n1, ncol = n1)
  I_n2n1 = matrix(1/n1, nrow = n2, ncol = n1)
  kermat.U.test = kermat.U.test - I_n2n1%*%kermat - kermat.U.test%*%I_n1 + I_n2n1%*%kermat%*%I_n1
  return(kermat.U.test)
}

krr = function(K.train, Y.train, lambda, K.test = NULL){
  n = dim(K.train)[1]
  alpha_hat = solve(K.train + lambda*diag(n)) %*% Y.train
  Y.test = NULL
  if(!is.null(K.test)){
    Y.test = K.test%*%alpha_hat
  }
  return(list(alpha_hat = alpha_hat, Y.test = Y.test))
}

predict.krr = function(ml.krr, K.test){
  # K.test: n.test by n.train
  Y.test = K.test%*%ml.krr$alpha_hat
  return(Y.test)
}

cv.krr = function(K, Y, lambda.vec, nfolds){
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  llam = length(lambda.vec)
  
  mse.list = list()
  for (k in 1:nfolds){
    K.train = K[unlist(flds[-k]), unlist(flds[-k])]
    K.val = K[unlist(flds[k]), unlist(flds[-k])] # n.test by n.train
    Y.train = Y[unlist(flds[-k])]
    Y.val = Y[unlist(flds[k])]
    Yhat.mtx = sapply(lambda.vec, function(lam) krr(K.train, Y.train, lam, K.val)$Y.test) # n.val by llam
    mse.list[[k]] = apply(sweep(Yhat.mtx, 1, Y.val), 2, function(col) mean(col^2))
  }
  mse.vec = apply(do.call(rbind,mse.list), 2, mean)
  lambda.min = lambda.vec[which.min(mse.vec)]
  best.ml = krr(K, Y, lambda.min)
  return(list(best.ml = best.ml, mse.vec = mse.vec, lambda.min = lambda.min))
}


regkrr = constructKRRLearner()
cv.regkrr = function(X, Y, gamma = NULL, nfolds = 10, lambda.vec){
  p = dim(X)[2]
  if (is.null(gamma)){
    gamma = 1/p
  }
  data.train = constructData(X, Y)
  params = constructParams(kernel="rbfdot", sigma=gamma, lambda=lambda.vec)
  best.para = CV(data.train, regkrr, params, fold = nfolds, verbose = F)
  m = regkrr$learn(data.train, best.para[[1]])
  best.lam = best.para[[1]]$lambda
  return(list(best.ml = m, lambda.min = best.lam))
}

X2U.kernel = function(X, K=50, gamma = NULL, plot = T){
  n = nrow(X)
  p = ncol(X)
  if(is.null(gamma)){
    gamma = 1/p
  }
  rbf = rbfdot(sigma = gamma)
  kermat = kernelMatrix(kernel = rbf, x = X)
  I_n = matrix(1/n, nrow = n, ncol = n)
  kermat = kermat - I_n%*%kermat - kermat%*%I_n + I_n%*%kermat%*%I_n
  kermat = kermat/n
  PCA.res = eigen(kermat)
  eigen_vals = PCA.res$values
  
  if (plot){
    par(mfrow=c(3,1))
    plot(eigen_vals[1:K])
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
    plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  }
  
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = PCA.res$vectors[,1:K]*sqrt(n)
  F1 = PCA.res$vectors[,1]*sqrt(n)
  F2 = PCA.res$vectors[,2]*sqrt(n)
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1))
}

X2U1.kernel = function(X, K=50, gamma = NULL, plot = F){
  n = nrow(X)
  p = ncol(X)
  if(is.null(gamma)){
    gamma = 1/p
  }
  rbf = rbfdot(sigma = gamma)
  kermat = kernelMatrix(kernel = rbf, x = X)
  I_n = matrix(1/n, nrow = n, ncol = n)
  kermat = kermat - I_n%*%kermat - kermat%*%I_n + I_n%*%kermat%*%I_n
  kermat = kermat/n
  PCA.res = eigen(kermat)
  eigen_vals = PCA.res$values
  
  if (plot){
    par(mfrow=c(3,1))
    plot(eigen_vals[1:K])
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
    plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))
  }
  
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  # F_ = PCA.res$vectors[,1:K]*sqrt(n)
  # F1 = PCA.res$vectors[,1]*sqrt(n)
  # F2 = PCA.res$vectors[,2]*sqrt(n)
  # P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  # H = diag(n) - P
  # return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1))
  
  F_ = cbind(1,PCA.res$vectors[,1:K]*sqrt(n))
  F1 = F_[,1]
  #F1 = PCA.res$vectors[,1]*sqrt(n)
  #F2 = PCA.res$vectors[,2]*sqrt(n)
  F2 = PCA.res$vectors[,1:2]*sqrt(n) # first 2 factors
  F_res = PCA.res$vectors[,3:n]*sqrt(n) # other factors
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  P1 = F1%*%solve(t(F1)%*%F1)%*%t(F1)
  H = diag(n) - P
  #L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  #L2 = solve(t(F2)%*%F2)%*%t(F2)%*%X
  #L_res = solve(t(F_res)%*%F_res)%*%t(F_res)%*%X
  #return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1, L = L, L2 = L2))
  return(list(H = H, P = P, P1 = P1, K = K, F_ = F_, F2 = F2))
}

cv.glmnet_ = function(X, Y, U, Y_, label, w, alpha = 0, lambda.vec, nfolds){
  flds = createFolds(label, k = nfolds, list = TRUE, returnTrain = FALSE)
  llam = length(lambda.vec)
  
  mse.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]),]
    X.val = X[unlist(flds[k]), ] 
    Y.train = Y[unlist(flds[-k])]
    Y.val = Y[unlist(flds[k])]
    U.train = U[unlist(flds[-k]),]
    U.val = U[unlist(flds[k]), ] 
    Y_.train = Y_[unlist(flds[-k])]
    Y_.val = Y_[unlist(flds[k])]
    w.train = w[unlist(flds[-k])]
    ml.ridge = glmnet(x = U.train, y = Y_.train, weights = w.train, lambda = lambda.vec, alpha = alpha)
    Yhat.mtx = predict(ml.ridge, newx = X.val) # n.train by llam
    mse.list[[k]] = apply(sweep(Yhat.mtx, 1, Y.val)^2, 2, mean)
  }
  mse.vec = apply(do.call(rbind, mse.list), 2, mean)
  lambda.min = lambda.vec[which.min(mse.vec)]
  best.ml = glmnet(x = U, y = Y_, weights = w, lambda = lambda.min, alpha = alpha)
  return(list(best.ml = best.ml, mse.vec = mse.vec, lambda.min = lambda.min))
}







