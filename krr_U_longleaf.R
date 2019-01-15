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
load("/nas/longleaf/home/peiyao/alpha/data/ADNI1.RData")

X = X[-c(167,770),]
label = label[-c(167,770)]
Y = Y[-c(167,770)]

X = X[-151,]
label = label[-151]
Y = Y[-151]

X = X[-c(112,169),]
label = label[-c(112,169)]
Y = Y[-c(112,169)]

X = X[-c(70,126),]
label = label[-c(70,126)]
Y = Y[-c(70,126)]

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

# groupwise design matrix
X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
n.train.vec = sapply(X.train.list, nrow)
n.test.vec = sapply(X.test.list, nrow)
n.train = sum(n.train.vec)
n.test = sum(n.test.vec)

# standardize Y
Y.train.mean = lapply(Y.train.list, mean)
Y.train.sd = lapply(Y.train.list, sd)

# standardize Y (subtract mean and divide sd)
Y.train.list = lapply(1:3, function(ix) (Y.train.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])
Y.test.list = lapply(1:3, function(ix) (Y.test.list[[ix]] - Y.train.mean[[ix]])/Y.train.sd[[ix]])

# standardize X (subtract mean)
X.train.mean = lapply(X.train.list, colMeans)
X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))

# data.frame format
data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))
data.krr.test.list = lapply(1:3, function(ix) constructData(y = Y.test.list[[ix]], x=X.test.list[[ix]]))


# ridge
ml.ridge.X.class = lapply(1:3, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
Yhat.ridge.X.class.test = lapply(1:3, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.X.class.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.ridge.X.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)

print(paste0("linear ridge within class performance:", as.character(mse.ridge.X.class)))

# kernel ridge
lambda.vec = exp(1)^seq(log(10^-4), log(10^1), length.out = 100)
gamma = gamma
ml.krr.X.class = lapply(1:3, function(ix) cv.regkrr(X.train.list[[ix]], Y.train.list[[ix]], gamma = gamma, lambda.vec = lambda.vec))
Yhat.krr.X.class.test = lapply(1:3, function(ix) regkrr$predict(ml.krr.X.class[[ix]]$best.ml, data.krr.test.list[[ix]]))
mse.krr.X.class.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.krr.X.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.krr.X.class = sum(mse.krr.X.class.vec*n.test.vec)/sum(n.test.vec)

print(paste0("kernel ridge within class performance:", as.character(mse.krr.X.class)))

# ------------------------------- ALPHA -----------------------------------------
mycut = X2U4(X.train.list, plot = F)
X2U.list = lapply(X.train.list, function(X)  X2U.cut(X, mycut))

H.list = lapply(X2U.list, function(list) list$H)
#F.list = lapply(X2U.list, function(list) list$F_)

U.list = lapply(1:3, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])

U.train = do.call(rbind, U.list)
Y_.train = do.call(c, Y_.list)
# PY.train = do.call(c, PY.list)
# data.U.train = data.frame(Y = Y_.train, U.train)

# ridge
ml.ridge.U = cv.glmnet(x=U.train, y=Y_.train, alpha = 0)
Yhat.ridge.U.test = lapply(1:3, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.ridge.U.test[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U = sum(mse.ridge.U.vec*n.test.vec)/sum(n.test.vec)

print(paste0("linear PCA linear ridge performance:", as.character(mse.ridge.U)))

# kernel ridge
gamma = gamma
data.U.train = constructData(U.train, Y_.train)
krr.U = constructKRRLearner()
lambda.vec = exp(1)^seq(log(10^-4), log(10^0), length.out = 100)
params = constructParams(kernel="rbfdot", sigma=gamma, lambda=lambda.vec)
best.para = CV(data.U.train, krr.U, params, fold = 10, verbose = T)
m = krr.U$learn(data.U.train, best.para[[1]])
data.test = constructData(X.test, Y.test)
Yhat.kridge.U.test = krr.U$predict(m, data.test)
Yhat.kridge.U.test = lapply(c(4,3,0), function(ix) Yhat.kridge.U.test[label.test==ix])
mse.kridge.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.kridge.U.test[[ix]]-Y.test.list[[ix]])^2))
mse.kridge.U = sum(mse.kridge.U.vec*n.test.vec)/sum(n.test.vec)

print(paste0("linear PCA kernel ridge performance:", as.character(mse.kridge.U)))
# ------------------------------ kernel PCA -----------------------------------------
gamma = gamma
mycut = X2U4.kernel(X.train.list, gamma = gamma, plot = F)
X2U.list = lapply(X.train.list, function(X)  X2U.kernel.cut(X, mycut))

H.list = lapply(X2U.list, function(list) list$H)
F.list = lapply(X2U.list, function(list) list$F_)
Y_.list = lapply(1:3, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
Y_.train = do.call(c, Y_.list)

# construct kernel matrix

# demeaned kernel matrix for training
K.mat = matrix(0, nrow = n.train, ncol = n.train)
ix.vec = c(0,cumsum(n.train.vec))
for (i in 1:2){
  for (j in (i+1):3){
    K.mat[(ix.vec[i]+1):ix.vec[i+1], (ix.vec[j]+1):ix.vec[j+1]] = X2K(X.train.list[[i]], X.train.list[[j]], gamma = gamma)
  }
}
K.mat = K.mat + t(K.mat)
for (i in 1:3){
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
for (i in 1:3){
  for (j in 1:3){
    K.mat.test[(ix.test.vec[i]+1):ix.test.vec[i+1], (ix.train.vec[j]+1):ix.train.vec[j+1]] = X2K.test(X.test.list[[i]], X.train.list[[j]], gamma = gamma)
  }
}
K.mat.test = K.mat.test%*%H_diag # n.test by n.train
Yhat.krr.U.test = lapply(c(4,3,0), function(ix) predict.krr(ml.krr$best.ml, K.mat.test[label.test==ix,]))
mse.krr.U.vec = sapply(1:3, function(ix) Y.train.sd[[ix]]^2*mean((Yhat.krr.U.test[[ix]]-Y.test.list[[ix]])^2))
mse.krr.U = sum(mse.krr.U.vec*n.test.vec)/sum(n.test.vec)

print(paste0("kernel PCA kernel ridge performance:", as.character(mse.krr.U)))

print(c(mse.ridge.X.class, mse.krr.X.class, mse.ridge.U, mse.kridge.U, mse.krr.U, myseed))

file.name = c("gamma_", as.character(-log10(gamma)),".csv")
file.name = paste(file.name, collapse ="")
write.table(t(c(mse.ridge.X.class, mse.krr.X.class, mse.ridge.U, mse.kridge.U, mse.krr.U, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)


  


