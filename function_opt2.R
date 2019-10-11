suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(Matrix)

return_beta_mtx = function(Y.train.WLS, F.train.list, L.list){
  n_label = length(Y.train.WLS)
  p = ncol(L.list[[1]])
  n.train.vec = sapply(1:n_label, function(ix) nrow(F.train.list[[ix]]))
  scale_matrix.list = lapply(1:n_label, function(ix) diag(n.train.vec[ix])/n.train.vec[ix])
  scale_matrix = bdiag(scale_matrix.list)
  
  beta = Variable(n_label*p)
  FL.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%L.list[[ix]])
  X = bdiag(FL.list)
  Y = do.call(c, Y.train.WLS)
  objective = Minimize(sum((scale_matrix%*%(Y - X %*% beta))^2))
  constraint = list(matrix(rep(diag(p), n_label), ncol = p*n_label)%*%beta == 0)
  problem <- Problem(objective, constraints = constraint)
  
  solution = solve(problem)
  beta.mtx = matrix(solution$getValue(beta), ncol = n_label)
  gamma.list = lapply(1:n_label, function(ix) L.list[[ix]]%*%beta.mtx[,ix])
  return(list(beta = beta.mtx, gamma.list = gamma.list))
}
