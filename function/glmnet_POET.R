library(CVXR)
library(caret)
library(scout)

## using scout
scout_POET = function(X, Y, Sigma, lambda = 0, alpha = 0){
  n = nrow(X)
  p = ncol(X)

  # xtx,xty,rho,thr=1e-4,maxit=100,beta.init=NULL
  if (alpha == 0){
    xtx = Sigma + 2*lambda*diag(p)
    rho = 0
  }else if (alpha == 1){
    xtx = Sigma
    rho = lambda
  }else{
    print('alpha can only take values 0 or 1')
  }
  xty = t(X)%*%Y/n
  scout.beta = crossProdLasso(xtx = xtx, xty = xty, rho = rho)$beta
  return(scout.beta)
}

glmnet.lambda_max= function(X, Y, alpha){
  n = nrow(X)
  max(abs(t(Y - mean(Y)*(1-mean(Y))) %*% X ))/(alpha * n) # largest lambda value
}

## using CVXR package
elastic_reg <- function(beta, lambda = 0, alpha = 0) {
  ridge <- (1 - alpha)*sum(beta^2)
  lasso <- alpha*p_norm(beta, 1)
  lambda*(lasso + ridge)
}

glmnet_POET <- function(X, Y, Sigma, lambda = 0, alpha = 0){
  n = nrow(X)
  p = ncol(X)
  beta = Variable(p)
  loss = quad_form(beta, Sigma)/2 - t(Y)%*%X%*%beta/n
  obj = loss + elastic_reg(beta, lambda, alpha)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  result$getValue(beta)
}

cv.glmnet_POET = function(X, Y, lambda, alpha = 0, C, nfolds = 10){
  flds = createFolds(Y, k = nfolds, list = TRUE, returnTrain = FALSE)
  MSE.list = list()
  for (k in 1:nfolds){
    X.train = X[unlist(flds[-k]), ]
    X.val = X[unlist(flds[k]), ]
    Y.train = matrix(Y[unlist(flds[-k])])
    Y.val = matrix(Y[unlist(flds[k])])
    
    SigmaU = POET(t(X.train), K = 0, C = C, thres = "soft", matrix = "vad")$SigmaU
    beta.list = lapply(lambda, function(lam) glmnet_POET(X.train, Y.train, SigmaU, lambda = lam, alpha = alpha))
    Yhat.list = lapply(beta.list, function(beta) X.val%*%beta)
    MSE.list[[k]] = sapply(Yhat.list, function(Yhat) mean((Y.val - Yhat)^2))
  }
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  lambda[which.min(rMSE)]
}
  