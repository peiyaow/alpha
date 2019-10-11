library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(CVST)
require(methods)
library(huge)
library(igraph)
library(gridExtra)

library(igraph)
library(dplyr)
library(ggplot2)
library(clime)
library(readr)
library(xtable)
library(POET)

source("~/Documents/GitHub/alpha/functions.R")
source("~/Documents/GitHub/alpha/sim_func.R")

load("~/Documents/GitHub/alpha/sim/parameter/sim_seeds.RData")
load("~/Documents/GitHub/alpha/sim/parameter/para_grid.RData")
load("~/Documents/GitHub/alpha/sim/parameter/para_123.RData")
setwd("~/Documents/GitHub/alpha/sim/123_hd/")

set.seed(111)
para = para_grid.list[[17]]

n = 100
p = para$p
n_label = 3
beta_F.list = lapply(1:n_label, function(ix) c(rep(ix, 5), rep(0, p-5)))

sigma.vec = para$sigma.vec
beta_U = c(c(rep(para$beta, 5), rep(0, p-5)))

n.train.vec = c(100, 100, 100)
n.test.vec = c(100, 100, 100)
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

sim_FnU.list = lapply(1:n_label, function(ix) simulateFnU(n, p, K.list[[ix]], mu_B.list[[ix]], Sigma_B.list[[ix]], Sigma_U))

# true loading matrix
L.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$L)
gamma.list = lapply(1:n_label, function(ix) L.list[[ix]]%*%beta_F.list[[ix]])

F.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$F_)
U.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$U)
F.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$F.test)
U.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$U.test)
X.train.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$X)
X.test.list = lapply(1:n_label, function(ix) sim_FnU.list[[ix]]$X.test)
tY.train.list = lapply(1:n_label, function(ix) F.list[[ix]]%*%gamma.list[[ix]] + U.list[[ix]]%*%beta_U)
tY.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%gamma.list[[ix]] + U.test.list[[ix]]%*%beta_U)
Y.train.list = lapply(1:n_label, function(ix) F.list[[ix]]%*%gamma.list[[ix]] + U.list[[ix]]%*%beta_U + rnorm(n)*sigma.vec[ix]^2)
Y.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%gamma.list[[ix]] + U.test.list[[ix]]%*%beta_U + rnorm(n)*sigma.vec[ix]^2)

# plot true Y and noisy Y
par(mfrow = c(3,3))
for (i in 1:3){
  plot(Y.test.list[[i]])
  plot(tY.test.list[[i]])
  plot(tY.test.list[[i]], Y.test.list[[i]])
}

X.train = do.call(rbind, X.train.list)
X.test = do.call(rbind, X.test.list)
Y.train = do.call(c, Y.train.list)
Y.test = do.call(c, Y.test.list)

# standardize X (subtract mean)
X.train.mean = lapply(X.train.list, colMeans)
X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))

X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))
H.list = lapply(X2U.list, function(list) list$H)
P.list = lapply(X2U.list, function(list) list$P)
L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 
K.list = lapply(X2U.list, function(list) list$K)

F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
U.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)
XF.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_[,-1]%*%L.list[[ix]])

U.train = do.call(rbind, U.train.list)
Sigma_u_hat = t(U.train)%*%U.train/sum(n.train.vec)
PCA.res = eigen(Sigma_u_hat)
A = PCA.res$vectors%*%diag(1/sqrt(PCA.res$values))

XA.train.list = lapply(1:n_label, function(ix) X.train.list[[ix]]%*%A)
X2U_new.list = lapply(1:n_label, function(ix) X2U1(XA.train.list[[ix]], plot = F))

H_new.list = lapply(X2U_new.list, function(list) list$H)
P_new.list = lapply(X2U_new.list, function(list) list$P)
L_new.list = lapply(X2U_new.list, function(list) matrix(list$L[-1,], ncol = p)) 
K_new.list = lapply(X2U_new.list, function(list) list$K)

F_new.train.list = lapply(1:n_label, function(ix) X2U_new.list[[ix]]$F_)
U_new.train.list = lapply(1:n_label, function(ix) X2U_new.list[[ix]]$U)
XF_new.train.list = lapply(1:n_label, function(ix) X2U_new.list[[ix]]$F_[,-1]%*%L_new.list[[ix]])

t(U_new.train.list[[2]])%*%U_new.train.list[[2]]
t(U.train.list[[2]])%*%U.train.list[[2]]

beta_tilde = t(A)%*%(Sigma_u_hat)%*%(A)%*%t(U.train%*%A)%*%Y.train
XFA.train.list = lapply(1:n_label, function(ix) XF.train.list[[ix]]%*%A)

XFA.train.list[[1]]%*%t(U.train%*%A)




