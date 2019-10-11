commonFL = function(X, K=35, plot = FALSE){
  # given concatXf output L
  n = nrow(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  if (plot){
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  }
  K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = PCA.res$vectors[, 1:K]*sqrt(n)
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  return(list(L = L, K = K, F = F_))
}

commonFL_given_K = function(X, K=10, plot = FALSE){
  # given concatXf output L
  n = nrow(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  if (plot){
    plot(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  }
  # K = which.max(eigen_vals[1:K]/eigen_vals[2:(K+1)])
  F_ = PCA.res$vectors[, 1:K]*sqrt(n)
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  return(list(L = L, K = K, F = F_))
}

FnU.svd2 = function(X, L){
  # X: standardized testing data (subtracting mean from X.train) 
  # L: loading matrix estimated from training data
  # Return F and U for testing data using svd
  n = nrow(X)
  res.svd = svd(L%*%t(X))
  F = sqrt(n)*res.svd$v%*%t(res.svd$u) # not including all 1 column
  F1 = cbind(1,F)
  P = F1%*%solve(t(F1)%*%F1)%*%t(F1)
  Xu = X - F%*%L
  U = (diag(n) - F%*%t(F)/n)%*%Xu
  R = Xu - U
  Xf = X - U
  return(list(F = F, F1 = F1, U = U, R = R, Xf = Xf, Xu = Xu, P = P))
}

myjive = function(X.list){
  n.vec = sapply(X.list, function(X) nrow(X))
  n_label = length(n.vec)
  a_prev = 4
  a = 5
  Xf.list = X.list
  while (abs(a-a_prev) > 1e-10){
    Xf = do.call(rbind, Xf.list)
    res = commonFL(Xf)
    L = res$L
    K = res$K
    FUP.list = lapply(1:n_label, function(l) FnU.svd2(X.list[[l]], L))
    Xf.list = lapply(1:n_label, function(l) FUP.list[[l]]$Xf)
    a_prev = a
    a = sum(sapply(1:n_label, function(l) sum(R.list[[l]]^2)))
  }
  F.list = lapply(1:n_label, function(l) FUP.list[[l]]$F)
  U.list = lapply(1:n_label, function(l) FUP.list[[l]]$U)
  P.list = lapply(1:n_label, function(l) FUP.list[[l]]$P)
  return(list(F.list = F.list, U.list = U.list, L = L, P.list = P.list, K = K))
}

allFL = function(X){
  n = nrow(X)
  PCA.res = eigen(X%*%t(X)/n)
  eigen_vals = PCA.res$values
  F_ = PCA.res$vectors*sqrt(n)
  L = solve(t(F_)%*%F_)%*%t(F_)%*%X
  return(list(L = L, F = F_))
}

FnU.svd1 = function(X, L){
  # X: standardized testing data (subtracting mean from X.train) 
  # L: loading matrix estimated from training data
  # Return F and U for testing data using svd
  n = nrow(X)
  if (nrow(L) == 0){ # no factors identified
    F1.test = as.matrix(rep(1, n))
    P = F1.test%*%solve(t(F1.test)%*%F1.test)%*%t(F1.test)
    U.test = X
  }else{
    res.svd = svd(L%*%t(X))
    F.test = sqrt(n)*res.svd$v%*%t(res.svd$u) # not including all 1 column
    F1.test = cbind(1, F.test) # including all 1 column
    P = F1.test%*%solve(t(F1.test)%*%F1.test)%*%t(F1.test)
    U.test = X - F.test%*%L
  }
  return(list(F_ = F1.test, U = U.test, P = P))
}

X_list = function(F, L){
  K = nrow(L)
  X.list = lapply(1:K, function(k) as.matrix(F[,k])%*%t(as.matrix(L[k,])))
  X_cum.list = list()
  for (k in 1:K){
    X_cum.list[[k]] = Reduce('+', X.list[1:k])
  }
  return(list(X = X.list, cumX = X_cum.list))
}
