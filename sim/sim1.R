library(MASS)
n = 100
p = 200
K = 3
mu = rep(0, p)
Sigma = diag(p)

F_all = mvrnorm(n, mu, Sigma)
F1 = F_all[,1:K]
F2 = F_all[,(K+1):p] 
# L1 = rbind(rep(1,p), rep(0.8,p), rep(0.6,p)) # K by p
L1 = rbind(c(1, 1,rep(0, p-2)), c(0, 1, 1, rep(0, p-3)), c(0, 0, 1,1, rep(0, p-4))) # K by p
L2 = (diag(p)*sqrt(2))[(K+1):p,]
#  matrix(rep(rep(0.2,p), p-K), nrow = p-K) # (p-K) by p

F_ = F1
U =  F2%*%L2
X = F_%*%L1 + U
beta_f = rep(2, p)
beta_u = rep(0.005, p)
gamma = L1%*%beta_f
L1%*%beta_f
L2%*%beta_u
Y = F_%*%gamma + U%*%beta_u
plot(abs(cor(X, Y)))

X.mean = apply(X, 2, mean)
X = sweep(X, 2, X.mean)

Fhat = X2F(X)$F_
values = X2F(X)$eigen_vals
plot(abs(cor(Fhat[,-1], Y)))
