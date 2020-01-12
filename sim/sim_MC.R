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

setwd("~/Documents/GitHub/alpha/")
load("./data/ADNI2_clean3.RData")
source("functions.R")
ROI_names <- read_csv("ROI_names.csv", col_names = FALSE)$X1

n = dim(X)[1]
p = dim(X)[2]
X.mean = apply(X, 2, mean)
X.sd =  apply(X, 2, sd)
X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
label.level = levels(label)
n_label = length(label.level)
X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) Y[label == l])
label.list = lapply(label.level, function(l) label[label == l])
n.vec = as.vector(table(label))
X.mean.list = lapply(X.list, colMeans)
X.list = lapply(1:n_label, function(ix) sweep(X.list[[ix]], 2, X.mean.list[[ix]]))

# choose a group
# para_1: X.list[[1]], K = 1
# para_2: X.list[[3]], K = 2
# para_3: X.list[[5]], K = 3

# myX = do.call(rbind, X.list[1:3])
myX = X.list[[4]]
n = 100
# ix = sample(nrow(myX), 100)
ix = 1:100
X = myX[ix,]
p = ncol(X)
p_trunc = 80
X2U.list = X2U1(X[,1:p_trunc], K = 10, plot = F)
K = X2U.list$K
L = matrix(X2U.list$L[-1,], ncol = p_trunc)
mu_B = apply(L, 1, mean)
Sigma_B = L%*%t(L)/p_trunc
U = X2U.list$H%*%X[,1:p_trunc]
#Sigma_U = t(U)%*%U/n
#Cor_U = cor(U)
#sd_U = apply(U, 2, sd)*sqrt((n-1)/n)
#diag(sd_U)%*%Cor_U%*%diag(sd_U) - Sigma_U
POET.list = POET(t(X[, 1:p_trunc]), K = X2U.list$K, C = 0.5, thres = "soft", matrix = "vad")
Sigma_U = POET.list$SigmaU
save(K, mu_B, Sigma_B, Sigma_U, file = "para_3_1.RData")



# generate list of data for the setting K1, K2, K3 = 3
mu_B.list = list()
Sigma_B.list = list()
K.list = list()
for (i in 1:3){
  load(paste0("para_3_",as.character(i),".RData"))
  mu_B.list[[i]]= mu_B
  Sigma_B.list[[i]] = Sigma_B
#  Sigma_U.list[[i]] = Sigma_U
  K.list[[i]] = K
}
#mu_B.list = mu_B.list[c(2,3,1)]
#Sigma_B.list = Sigma_B.list[c(2,3,1)]
#Sigma_U.list = Sigma_U.list[c(2,3,1)]
save(K.list, mu_B.list, Sigma_B.list, Sigma_U, file = "para_333.RData")

# generate list of data for the setting K1 = 1, K2 = 2, K3 = 3
mu_B.list = list()
Sigma_B.list = list()
K.list = list()
for (i in 1:3){
  load(paste0("para_", as.character(i), "_1.RData"))
  mu_B.list[[i]]= mu_B
  Sigma_B.list[[i]] = Sigma_B
  #  Sigma_U.list[[i]] = Sigma_U
  K.list[[i]] = K
}
#mu_B.list = mu_B.list[c(2,3,1)]
#Sigma_B.list = Sigma_B.list[c(2,3,1)]
#Sigma_U.list = Sigma_U.list[c(2,3,1)]
save(K.list, mu_B.list, Sigma_B.list, Sigma_U, file = "para_123.RData")

# save seeds
seeds = floor(runif(20)*10000)
save(seeds, file = "sim_seeds.RData")




