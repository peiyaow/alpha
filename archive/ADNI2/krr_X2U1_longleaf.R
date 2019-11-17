# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
# cat(args, "\n")
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(CVST)
require(methods)

source("/nas/longleaf/home/peiyao/alpha/functions.R")
load("/nas/longleaf/home/peiyao/alpha/data/ADNI.RData")

# X = X[label!=3,]
# Y = Y[label!=3]
# label = label[label!=3]
# label = droplevels(label)

# X = X[label!=4,]
# Y = Y[label!=4]
# label = label[label!=4]
# label = droplevels(label)

# X = X[label!=3,]
# Y = Y[label!=3]
# label = label[label!=3]
# label = droplevels(label)

Y = log(Y+1)

n = dim(X)[1]
p = dim(X)[2]

set.seed(myseed)

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

# scale the data overall
X.train.mean = apply(X.train, 2, mean)
X.train.sd =  apply(X.train, 2, sd)
X.train = sweep(sweep(X.train, 2, X.train.mean), 2, X.train.sd, "/")
X.test = sweep(sweep(X.test, 2, X.train.mean), 2, X.train.sd, "/")

label.level = levels(label)
n_label = length(label.level)

# groupwise design matrix
X.train.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])
Y.train.list = lapply(label.level, function(l) Y.train[label.train == l])
Y.test.list = lapply(label.level, function(l) Y.test[label.test == l])
n.train.vec = sapply(X.train.list, nrow)
n.test.vec = sapply(X.test.list, nrow)
n.train = sum(n.train.vec)
n.test = sum(n.test.vec)

# # standardize Y
# Y.train.mean = lapply(Y.train.list, mean)
# Y.train.sd = lapply(Y.train.list, sd)
# 
# # standardize Y (subtract mean and divide sd)
# Y.train.list = lapply(1:n_label, function(ix) (Y.train.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])
# Y.test.list = lapply(1:n_label, function(ix) (Y.test.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])

# standardize X (subtract mean)
X.train.mean = lapply(X.train.list, colMeans)
X.train.list = lapply(1:n_label, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:n_label, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))

# data.frame format
data.train.list = lapply(1:n_label, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
data.test.list = lapply(1:n_label, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
data.krr.test.list = lapply(1:n_label, function(ix) constructData(y = Y.test.list[[ix]], x=X.test.list[[ix]]))

Y.train.mean = lapply(Y.train.list, mean)
Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]] - Y.train.mean[[l]]))
X.train.WLS = do.call(rbind, X.train.list)
data.train.WLS = data.frame(Y = Y.train.WLS, X.train.WLS)
ml.lm.WLS = lm(Y~., data = data.train.WLS)
ix.vec = c(0,cumsum(n.train.vec))
sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.WLS$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
ml.ridge.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 0, weights = w)
Yhat.ridge.WLS.test = lapply(1:n_label, function(ix) predict(ml.ridge.WLS, s=ml.ridge.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.WLS.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.WLS.test[[ix]]+Y.train.mean[[ix]])-exp(Y.test.list[[ix]]))^2))
# mse.ridge.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))

mse.ridge.WLS = sum(mse.ridge.WLS.vec*n.test.vec)/sum(n.test.vec)

# ridge
ml.ridge.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
Yhat.ridge.X.class.test = lapply(1:n_label, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
#mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.X.class.test[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)

# kernel ridge
lambda.vec = exp(1)^seq(log(10^-4), log(10^1), length.out = 100)
gamma = gamma
ml.krr.X.class = lapply(1:n_label, function(ix) cv.regkrr(X.train.list[[ix]], Y.train.list[[ix]], gamma = gamma, lambda.vec = lambda.vec))
Yhat.krr.X.class.test = lapply(1:n_label, function(ix) regkrr$predict(ml.krr.X.class[[ix]]$best.ml, data.krr.test.list[[ix]]))
mse.krr.X.class.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.krr.X.class.test[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.krr.X.class = sum(mse.krr.X.class.vec*n.test.vec)/sum(n.test.vec)

# ------------------------------- ALPHA -----------------------------------------
X2U.list = lapply(X.train.list, function(X)  X2U1(X, plot = F))

H.list = lapply(X2U.list, function(list) list$H)
P1.list = lapply(X2U.list, function(list) list$P1)
K.list = lapply(X2U.list, function(list) list$K)

U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
Y_.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
Y_mean.list = lapply(1:n_label, function(ix) (P1.list[[ix]]%*%Y.train.list[[ix]])[1])

U.train = do.call(rbind, U.list)
Y_.train = do.call(c, Y_.list)
data.U.train = data.frame(Y = Y_.train, U.train)

# lm and wls ridge
ml.lm.U = lm(Y~., data = data.U.train)  
Yhat.lm.U.train = lapply(1:n_label, function(ix) predict(ml.lm.U, newdata = data.train.list[[ix]][, -1])+Y_mean.list[[ix]])

# weight computed from sx
sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train[[l]])^2))
w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, weights = w, alpha = 0)
Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
# mse.ridge.sx.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.sx.U.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.ridge.sx.U = sum(mse.ridge.sx.U.vec*n.test.vec)/sum(n.test.vec)

# weight computed from su
ix.vec = c(0,cumsum(n.train.vec))
sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, weights = w, alpha = 0)
Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
# mse.ridge.su.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.su.U.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.ridge.su.U = sum(mse.ridge.su.U.vec*n.test.vec)/sum(n.test.vec)

# regular ridge no weight
ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, alpha = 0)
Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
#mse.ridge.U.vec = sapply(1:n_label, function(ix) mean((Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.ridge.U = sum(mse.ridge.U.vec*n.test.vec)/sum(n.test.vec)

# weight computed from diag(H) individual
w = 1/do.call(c, sapply(H.list, function(H) diag(H)))
ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, weights = w, alpha = 0)
Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
#mse.ridge.wHi.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.wHi.U.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.ridge.U.test[[ix]]+Y_mean.list[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.ridge.wHi.U = sum(mse.ridge.wHi.U.vec*n.test.vec)/sum(n.test.vec)
# ------------------------------ kernel PCA -----------------------------------------
gamma = gamma
X2U.list = lapply(X.train.list, function(X)  X2U1.kernel(X, gamma = gamma))

H.list = lapply(X2U.list, function(list) list$H)
P1.list = lapply(X2U.list, function(list) list$P1)
Y_.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
Y_mean.list = lapply(1:n_label, function(ix) (P1.list[[ix]]%*%Y.train.list[[ix]])[1])

Y_.train = do.call(c, Y_.list)

# construct kernel matrix
# demeaned kernel matrix for training
K.mat = matrix(0, nrow = n.train, ncol = n.train)
ix.vec = c(0,cumsum(n.train.vec))

for (i in 1:(n_label-1)){
  for (j in (i+1):n_label){
    K.mat[(ix.vec[i]+1):ix.vec[i+1], (ix.vec[j]+1):ix.vec[j+1]] = X2K(X.train.list[[i]], X.train.list[[j]], gamma = gamma)
  }
}

K.mat = K.mat + t(K.mat)

for (i in 1:n_label){
  K.mat[(ix.vec[i]+1):ix.vec[i+1], (ix.vec[i]+1):ix.vec[i+1]] = X2K(X.train.list[[i]], X.train.list[[i]], gamma = gamma)
}

H_diag = as.matrix(bdiag(H.list))
K.mat = H_diag%*%K.mat%*%H_diag

# kernel PCA and corresponding kernel ridge
lambda.vec = ml.ridge.U$lambda
ml.krr = cv.krr(K.mat, Y_.train, lambda.vec, nfolds = 10)

# construct kernel test matrix
K.mat.test = matrix(0, nrow = n.test, ncol = n.train)
ix.test.vec = c(0,cumsum(n.test.vec))
ix.train.vec = c(0,cumsum(n.train.vec))
for (i in 1:n_label){
  for (j in 1:n_label){
    K.mat.test[(ix.test.vec[i]+1):ix.test.vec[i+1], (ix.train.vec[j]+1):ix.train.vec[j+1]] = X2K.test(X.test.list[[i]], X.train.list[[j]], gamma = gamma)
  }
}
K.mat.test = K.mat.test%*%H_diag # n.test by n.train

Yhat.krr.U.test = predict.krr(ml.krr$best.ml, K.mat.test)
Yhat.krr.U.test = lapply(1:n_label, function(ix) Yhat.krr.U.test[(ix.test.vec[ix]+1):ix.test.vec[ix+1]])
mse.krr.U.vec = sapply(1:n_label, function(ix) mean((exp(Yhat.krr.U.test[[ix]]+Y_mean.list[[ix]])-exp(Y.test.list[[ix]]))^2))
mse.krr.U = sum(mse.krr.U.vec*n.test.vec)/sum(n.test.vec)

file.name = c("gamma_", as.character(-log10(gamma)),".csv")
file.name = paste(file.name, collapse ="")
write.table(t(c(mse.ridge.WLS, mse.ridge.X.class, mse.krr.X.class, mse.ridge.U, mse.ridge.sx.U, mse.ridge.su.U, mse.ridge.wHi.U, mse.krr.U, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = c("class_gamma_", as.character(-log10(gamma)),".csv")
file.name = paste(file.name, collapse ="")
write.table(t(c(mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.krr.X.class.vec, mse.ridge.U.vec, mse.ridge.sx.U.vec, mse.ridge.su.U.vec, mse.ridge.wHi.U.vec, mse.krr.U.vec, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)




