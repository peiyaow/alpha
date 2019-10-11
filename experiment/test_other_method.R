library(rlist)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(CVST)
require(methods)

setwd("~/Documents/GitHub/alpha/")
load("./data/ADNI2_clean3.RData")
source("functions.R")
source("jive_function.R")

n = dim(X)[1]
p = dim(X)[2]

ix.train = unlist(createDataPartition(label, times = 1, p = 3/4))
ix.test = (1:n)[-ix.train]
ix.list = list(ix.train, ix.test)
Y.list = lapply(1:2, function(x) Y[ix.list[[x]]]) # train,test
X.list = lapply(1:2, function(x) X[ix.list[[x]],])
label.list = lapply(1:2, function(x) label[ix.list[[x]]])

X.train = X.list[[1]]
X.test = X.list[[2]]
Y.train = Y.list[[1]]
Y.test = Y.list[[2]]
label.train = label.list[[1]]
label.test = label.list[[2]]

X.train = X.train[label.train < 4,]
Y.train = Y.train[label.train < 4]
X.test = X.test[label.test < 4,]
Y.test = Y.test[label.test < 4]
label.train = droplevels(label.train[label.train < 4])
label.test = droplevels(label.test[label.test < 4])

# scale the data overall
X.train.mean = apply(X.train, 2, mean)
X.train.sd =  apply(X.train, 2, sd)
X.train = sweep(sweep(X.train, 2, X.train.mean), 2, X.train.sd, "/")
X.test = sweep(sweep(X.test, 2, X.train.mean), 2, X.train.sd, "/")


label.level = levels(label.train)
n_label = length(label.level)

# groupwise design matrix
X.train.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])
Y.train.list = lapply(label.level, function(l) Y.train[label.train == l])
Y.test.list = lapply(label.level, function(l) Y.test[label.test == l])
label.train.list = lapply(label.level, function(l) label.train[label.train == l])




n.train.vec = sapply(X.train.list, nrow)
n.test.vec = sapply(X.test.list, nrow)
n.vec = as.vector(table(label))
n.train = sum(n.train.vec)
n.test = sum(n.test.vec)



# standardize X (subtract mean)
X.train.mean = lapply(X.train.list, colMeans)
X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))

Y.train.mean = lapply(Y.train.list, mean)
Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]] - Y.train.mean[[l]]))
X.train = do.call(rbind, X.train.list)
data.train.WLS = data.frame(Y = Y.train.WLS, X.train)

# WLS
ml.lm.WLS = lm(Y~., data = data.train.WLS)
Yhat.lm.WLS.test = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.test.list[[ix]])))
mse.lm.WLS.vec = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS = sum(mse.lm.WLS.vec*n.test.vec)/sum(n.test.vec)


X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))
H.list = lapply(X2U.list, function(list) list$H)
P.list = lapply(X2U.list, function(list) list$P)
K.list = lapply(X2U.list, function(list) list$K)
L.list =  lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) # loading matrix
F.train.list = lapply(1:n_label, function(ix) as.matrix(X2U.list[[ix]]$F_[,-1]))

Y.train.list
Y.train.WLS = lapply(1:n_label, function(l) Y.train.list[[l]]- Y.train.mean[[l]])

data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = Y.train.WLS[[ix]], 
                                                              F.train.list[[ix]]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))

fitted = do.call(c, lapply(1:n_label, function(ix) ml.lm.F.list[[ix]]$fitted.values))

plot(fitted)
plot(do.call(c, Y.train.WLS))
plot(do.call(c, Y.train.WLS) - fitted)
gamma.list = return_beta_mtx(Y.train.WLS, F.train.list, L.list)$gamma.list
t(gamma.list[[1]])%*%gamma.list[[1]]

U.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)
U.train = do.call(rbind, U.train.list)
sigma_4 = t(U.train.list[[4]])%*%U.train.list[[4]]/n.train.vec[4]
sigma_u = t(U.train)%*%U.train/n.train


HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
HY.train = do.call(c, HY.train.list)
data.U.train = data.frame(Y = HY.train, U.train)
ml.lm.U = lm(Y~., data = data.U.train)
ml.lm.U$coefficients[-1]%*%sigma_4%*%as.matrix(ml.lm.U$coefficients[-1])

