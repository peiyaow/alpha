library(kernlab)
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
  F1 = PCA.res$vectors[,1]*sqrt(n)
  F2 = PCA.res$vectors[,2]*sqrt(n)
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  return(list(H = H, P = P, K = K, F_ = F_, F2 = F2, F1 = F1))
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
  
  return(cut.eigen)
}

X2U.cut = function(X, cut){
  n = nrow(X)
  #  p = ncol(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  K = n - length(eigen_vals[eigen_vals<cut])
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
  
  K = n - length(eigen_vals[eigen_vals<cut])
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

  
X2U2 = function(X, K=50, plot = T){
  n = nrow(X)
  PCA.res = eigen(X%*%t(X))
  
  eigen_vals = PCA.res$values
  if (plot){
    par(mfrow=c(2,1))
    plot(eigen_vals[1:K])
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  }
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = PCA.res$vectors[,1:K]*sqrt(n)
  P = F_%*%solve(t(F_)%*%F_)%*%t(F_)
  H = diag(n) - P
  return(list(H = H, P = P, K = K))
}


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
