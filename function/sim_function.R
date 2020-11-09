library(MASS)
library(POET)
library(pracma)

SOLVE = function(x){ # input must be matrix
  if (sum(dim(x))){
    return(solve(x))
  }else{
    return(x)
  }
}

DIAG = function(e){
  if (length(e) > 1){
    return(diag(e))
  }else{
    return(matrix(e))
  }
}

FactorModelPara = function(n = 100, p = 200, K = 3, spike = c(40, 20, 4)/2, d = 0.1, du = 0.3, rrho = 0.1){
  # K: number of factors
  # spike: K spiked eigenvalues
  # d: dilute parameter for L2
  # du: dilute parameter for Sigma_U
  
  R = matrix(, nrow = K, ncol = K)
  for (i in 1:K){
    for (j in 1:K){
      R[i,j] = rrho^(abs(i-j))
    }
  }
  LtL = R*sqrt(spike)%*%t(sqrt(spike))
  eigen.res = eigen(LtL)
  V = eigen.res$vectors
  Q = randortho(K, type = "orthonormal")
  L1 = Q%*%sqrt(DIAG(eigen.res$values))%*%t(V)
  
  TT = matrix(runif((p-K)*K, min = -d, max = d), nrow = K)
  L2 = Q%*%TT
  
  L = cbind(L1, L2) # K by p
  SigmaX = t(L)%*%L
  
  Sigma_U = du*diag(p)
  
  return(list(L = L, SigmaU = Sigma_U, SigmaX = SigmaX))
}

generateX = function(n = 100, p = 200, K = 3, spike = c(40, 20, 4)/2, d = 0.1, du = 0.3, rrho = 0.1){
  # K: number of factors
  # spike: K spiked eigenvalues
  # d: dilute parameter for L2
  # du: dilute parameter for Sigma_U
  
  R = matrix(, nrow = K, ncol = K)
  for (i in 1:K){
    for (j in 1:K){
      R[i,j] = rrho^(abs(i-j))
    }
  }
  LtL = R*sqrt(spike)%*%t(sqrt(spike))
  svd_res = svd(LtL)
  D = sqrt(diag(svd_res$d))%*%t(svd_res$v)
  Q = svd_res$u #Q
  if (rrho == 0){
    Q = randortho(K, type = "orthonormal")
  }
  
  TT = d*matrix(runif((p-K)*K), nrow = K)
  L2 = Q%*%TT
  
  # explode L1 by sqrt(spike)
  L1 = Q%*%D
  
  L = cbind(L1, L2) # K by p
  SigmaX = t(L)%*%L
  
  Sigma_U = du*diag(p)
  
  F1 = mvrnorm(n, rep(0, K), diag(K))
  U1 = mvrnorm(n, rep(0, p), Sigma_U)
  
  X1 = F1%*%L + U1
  return(list(X = X1, F = F1, U = U1, L = L, SigmaU = Sigma_U, SigmaX = SigmaX))
}

# generateX = function(K = 3, p = 200, spike = c(40, 20, 4)/2, d = 0.1, du = 0.3){
#   # K: number of factors
#   # spike: K spiked eigenvalues
#   # d: dilute parameter for L2
#   # du: dilute parameter for Sigma_U
#   
#   L1 = randortho(K, type = "orthonormal")
#   TT = d*matrix(runif((p-K)*K), nrow = K)
#   L2 = L1%*%TT
#   
# # explode L1 by sqrt(spike)
#   L1 = t(t(L1)*sqrt(spike))
#   
#   L = cbind(L1, L2) # K by p
#   SigmaX = t(L)%*%L
#   
#   Sigma_U = du*diag(p)
#   
#   F1 = mvrnorm(n, rep(0, K), diag(K))
#   U1 = mvrnorm(n, rep(0, p), Sigma_U)
#   
#   X1 = F1%*%L + U1
#   return(list(X = X1, F = F1, U = U1, L = L, SigmaU = Sigma_U, SigmaX = SigmaX))
# }



simulateFnU.rd0 = function(n, p, K, mu_B, Sigma_B, Sigma_U){
  L = t(mvrnorm(p, mu_B, Sigma_B))
  F_ = mvrnorm(n, mu = rep(0, K), diag(K))
  U = mvrnorm(n, mu = rep(0, p), Sigma_U)
  F_.test = mvrnorm(n, mu = rep(0, K), diag(K))
  U.test = mvrnorm(n, mu = rep(0, p), Sigma_U)
  X = F_%*%L + U
  X.test = F_.test%*%L + U.test
  return(list(X = X, X.test = X.test, L = L, F_ = F_, U = U, F.test = F_.test, U.test = U.test))
}

simulateFnU.hd0 = function(n, p, K, mu_B, Sigma_B, Sigma_U){
  L.list = list()
  Sigma_U.list = list()
  L0 = t(mvrnorm(80, mu_B, Sigma_B))
  for (i in 1:(p/80)){
    L.list[[i]] = L0
    Sigma_U.list[[i]] = Sigma_U
  }
  L = do.call(cbind, L.list)
  bSigma_U = bdiag(Sigma_U.list)
  
  F_ = mvrnorm(n, mu = rep(0, K), diag(K))  
  U = mvrnorm(n, mu = rep(0, p), bSigma_U)
  F_.test = mvrnorm(n, mu = rep(0, K), diag(K))
  U.test = mvrnorm(n, mu = rep(0, p), bSigma_U)
  X = F_%*%L + U
  X.test = F_.test%*%L + U.test
  
  return(list(X = X, X.test = X.test, L = L, L0 = L0, F_ = F_, U = U, F.test = F_.test, U.test = U.test))
}

simulateFnU0 = function(n, p, K, mu_B, Sigma_B, Sigma_U){
  if (p > 80){
    return(simulateFnU.hd(n, p, K, mu_B, Sigma_B, Sigma_U))
  }else{
    return(simulateFnU.rd(n, p, K, mu_B, Sigma_B, Sigma_U))
  }
}

simulateFnU.rd = function(n, p, K, L, Sigma_U){
  F_ = mvrnorm(n, mu = rep(0, K), diag(K))
  U = mvrnorm(n, mu = rep(0, p), Sigma_U)
  F_.test = mvrnorm(n, mu = rep(0, K), diag(K))
  U.test = mvrnorm(n, mu = rep(0, p), Sigma_U)
  X = F_%*%L + U
  X.test = F_.test%*%L + U.test
  return(list(X = X, X.test = X.test, L = L, F_ = F_, U = U, F.test = F_.test, U.test = U.test))
}

simulateFnU.hd = function(n, p, K, L0, Sigma_U){
  L.list = list()
  Sigma_U.list = list()
  for (i in 1:(p/80)){
    L.list[[i]] = L0
    Sigma_U.list[[i]] = Sigma_U
  }
  L = do.call(cbind, L.list)
  bSigma_U = bdiag(Sigma_U.list)
  
  F_ = mvrnorm(n, mu = rep(0, K), diag(K))  
  U = mvrnorm(n, mu = rep(0, p), bSigma_U)
  F_.test = mvrnorm(n, mu = rep(0, K), diag(K))
  U.test = mvrnorm(n, mu = rep(0, p), bSigma_U)
  X = F_%*%L + U
  X.test = F_.test%*%L + U.test
  return(list(X = X, X.test = X.test, L = L, L0 = L0, F_ = F_, U = U, F.test = F_.test, U.test = U.test))
}

simulateFnU = function(n, p, K, L0, Sigma_U){
  if (p > 80){
    return(simulateFnU.hd(n, p, K, L0, Sigma_U))
  }else{
    return(simulateFnU.rd(n, p, K, L0, Sigma_U))
  }
}

square_root = function(a){
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
  a.sqrtinv <- a.eig$vectors %*% diag(1/sqrt(a.eig$values)) %*% solve(a.eig$vectors)
  return(list(sqrt=a.sqrt,sqrt.inv=a.sqrtinv))
}

denman.beavers <- function(mat,maxit=50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1,nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5*(y+solve(z))
    z <- 0.5*(z+solve(y))
    y <- y.temp
  }
  return(list(sqrt=y,sqrt.inv=z))
}

est_sigma2 = function(Y, X, Sigma){
  p = dim(X)[2]
  n = dim(X)[1]
  Y = as.matrix(Y)
  A = (p + n + 1)/(n*(n+1))*t(Y)%*%Y
  B = 1/(n*(n+1))*(t(Y)%*%X%*%solve(Sigma)%*%t(X)%*%Y)
  #    t(Sigma.sqrtinv%*%t(X)%*%Y)%*%(Sigma.sqrtinv%*%t(X)%*%Y)
  return(list(A = A, B = B, sigma2 = A - B))
}

var_moment_est = function(Y, X, Sigma){
  p = dim(X)[2]
  n = dim(X)[1]
  Y = as.matrix(Y)
  A = (p)/(n*(n+1))*t(Y)%*%Y
  B = 1/(n*(n+1))*(t(Y)%*%X%*%solve(Sigma)%*%t(X)%*%Y)
  return(list(A = A, B = B))
}

POET = function (Y, K = -Inf, C = -Inf, thres = "soft", matrix = "cor") 
{
  p = nrow(Y)
  n = ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean))) %*% matrix(1, 1, n)
  if (K == -Inf) {
    K1 = 0.25 * (POETKhat(Y)$K1HL + POETKhat(Y)$K2HL + POETKhat(Y)$K1BN + 
                   POETKhat(Y)$K2BN)
    K = floor(K1) + 1
  }
  if (K > 0) {
    V <- eigen(t(Y) %*% Y)$vectors
    V = as.matrix(V)
    Dd <- eigen(t(Y) %*% Y)$values
    Dd = as.vector(Dd)
    W <- sort(diag(Dd), index.return = TRUE)$x
    W = as.matrix(W)
    Id <- sort(diag(Dd), index.return = TRUE)$ix
    Id = as.matrix(Id)
    F <- sqrt(n) * V[, 1:K]
    LamPCA = Y %*% F/n
    uhat = Y - LamPCA %*% t(F)
    Lowrank = LamPCA %*% t(LamPCA)
    rate = 1/sqrt(p) + sqrt((log(p))/n)
  }
  else {
    uhat = Y
#    rate = sqrt((log(p))/n)
    rate = 1/sqrt(p) + sqrt((log(p))/n)
    Lowrank = matrix(0, p, p)
    LamPCA = matrix(NA, nrow = p, ncol = n)
  }
  SuPCA = uhat %*% t(uhat)/n
  SuDiag = diag(diag(SuPCA))
  if (matrix == "cor") {
    R = solve(SuDiag^(1/2)) %*% SuPCA %*% solve(SuDiag^(1/2))
  }
  if (matrix == "vad") {
    R = SuPCA
  }
  if (C == -Inf) {
    C1 = POETCmin(Y, K, thres, matrix)
    C = C1 + 0.1
  }
  uu = array(0, dim = c(p, p, n))
  roottheta = array(0, dim = c(p, p))
  lambda = array(0, dim = c(p, p))
  for (i in 1:p) {
    for (j in 1:i) {
      uu[i, j, ] = uhat[i, ] * uhat[j, ]
      roottheta[i, j] = sd(uu[i, j, ])
      lambda[i, j] = roottheta[i, j] * rate * C
      lambda[j, i] = lambda[i, j]
    }
  }
  Rthresh = matrix(0, p, p)
  if (thres == "soft") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (abs(R[i, j]) < lambda[i, j] && j < i) {
          Rthresh[i, j] = 0
        }
        else {
          if (j == i) {
            Rthresh[i, j] = R[i, j]
          }
          else {
            Rthresh[i, j] = sign(R[i, j]) * (abs(R[i, 
                                                   j]) - lambda[i, j])
          }
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  if (thres == "hard") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (abs(R[i, j]) < lambda[i, j] && j < i) {
          Rthresh[i, j] = 0
        }
        else {
          Rthresh[i, j] = R[i, j]
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  if (thres == "scad") {
    for (i in 1:p) {
      for (j in 1:i) {
        if (j == i) {
          Rthresh[i, j] = R[i, j]
        }
        else {
          if (abs(R[i, j]) < lambda[i, j]) {
            Rthresh[i, j] = 0
          }
          else {
            if (abs(R[i, j]) < 2 * lambda[i, j]) {
              Rthresh[i, j] = sign(R[i, j]) * (abs(R[i, 
                                                     j]) - lambda[i, j])
            }
            else {
              if (abs(R[i, j]) < 3.7 * lambda[i, j]) {
                Rthresh[i, j] = ((3.7 - 1) * R[i, j] - 
                                   sign(R[i, j]) * 3.7 * lambda[i, j])/(3.7 - 
                                                                          2)
              }
              else {
                Rthresh[i, j] = R[i, j]
              }
            }
          }
        }
        Rthresh[j, i] = Rthresh[i, j]
      }
    }
  }
  SigmaU = matrix(0, p, p)
  if (matrix == "cor") {
    SigmaU = SuDiag^(1/2) %*% Rthresh * SuDiag^(1/2)
  }
  if (matrix == "vad") {
    SigmaU = Rthresh
  }
  SigmaY = SigmaU + Lowrank
  result <- list(SigmaU = SigmaU, SigmaY = SigmaY, factors = t(F), 
                 loadings = LamPCA)
  return(result)
}
