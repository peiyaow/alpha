suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(Matrix)

#K 
#n_label
return_gamma_mtx = function(res.Y.train.list, F.train.list){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)
  
  
  
  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}


return_gamma_mtx_lasso = function(res.Y.train.list, F.train.list, lambda){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)

  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda * p_norm(gamma, 1))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}

return_gamma_mtx_lasso_no_constraint = function(res.Y.train.list, F.train.list, lambda){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)

  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda * p_norm(gamma, 1))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  # constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}

return_gamma_mtx_ridge = function(res.Y.train.list, F.train.list, weight, lambda){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) weight[ix]*diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)

  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda*sum(gamma^2))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}

return_gamma_mtx_EN = function(res.Y.train.list, F.train.list, lambda, alpha){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)
  
  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda*(1-alpha)* sum(gamma^2) + lambda * alpha * p_norm(gamma, 1))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}

return_gamma_mtx_EN = function(res.Y.train.list, F.train.list, weight, lambda, alpha){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)
  
  # weight_matrix for gamma for the penalty
  weight_matrix.list = lapply(1:n_label, function(ix) diag(K)*weight[ix])
  weight_matrix = bdiag(weight_matrix.list)
  
  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda*(1-alpha)* sum((sqrt(weight_matrix)%*%gamma)^2) + lambda * alpha * p_norm(weight_matrix%*%gamma, 1))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}

return_gamma_mtx_EN = function(res.Y.train.list, F.train.list, weight, lambda, alpha){
  n_label = length(Y.train.list)
  K = ncol(F.train.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)
  
  # weight_matrix for gamma for the penalty
  weight_matrix.list = lapply(1:n_label, function(ix) diag(K)*weight[ix])
  weight_matrix = bdiag(weight_matrix.list)
  
  gamma = Variable(n_label*K)
  X = bdiag(F.train.list)
  Y = do.call(c, res.Y.train.list)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% gamma))^2) + lambda*(1-alpha)* sum((gamma)^2) + lambda * alpha * p_norm(weight_matrix%*%gamma, 1))
  # objective = Minimize(sum((Y - X %*% gamma)^2))
  constraint = list(matrix(rep(diag(K), n_label), ncol = K*n_label)%*%gamma == 0)
  problem <- Problem(objective, constraints = constraint)
  solution = solve(problem)
  gamma.mtx = matrix(solution$getValue(gamma), ncol = n_label)
  #gamma.list = apply(gamma.mtx, 2, function(x) x)
  return(list(gamma = gamma.mtx))
}


