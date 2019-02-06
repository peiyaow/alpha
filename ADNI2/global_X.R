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
load("/nas/longleaf/home/peiyao/alpha/data/ADNI2_clean.RData")

X = X[label!=4,]
Y = Y[label!=4]
label = label[label!=4]
label = droplevels(label)

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

# ridge
ml.ridge.X = cv.glmnet(x=X.train, y=Y.train, alpha = 0)
Yhat.ridge.X.test = predict(ml.ridge.X, s = ml.ridge.X$lambda.min, newx = X.test)
mse.ridge.X.vec = sapply(label.level, function(l) mean((Yhat.ridge.X.test[label.test==l] - Y.test[label.test==l])^2))
mse.ridge.X = mean((Yhat.ridge.X.test - Y.test)^2)

# kernel ridge
lambda.vec = exp(1)^seq(log(10^-4), log(10^1), length.out = 100)
gamma = 0.0001
ml.krr.X = cv.regkrr(X.train, Y.train, gamma = gamma, lambda.vec = lambda.vec)
Yhat.krr.X.test = regkrr$predict(ml.krr.X$best.ml, constructData(y=Y.test, x=X.test))
mse.krr.X.vec = sapply(label.level, function(l) mean((Yhat.krr.X.test[label.test==l]-Y.test[label.test==l])^2))
mse.krr.X = sum(mse.krr.X.class.vec*n.test.vec)/sum(n.test.vec)


file.name = c("gamma_", as.character(-log10(gamma)),".csv")
file.name = paste(file.name, collapse ="")
write.table(t(c(mse.ridge.X, mse.krr.X, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = c("class_gamma_", as.character(-log10(gamma)),".csv")
file.name = paste(file.name, collapse ="")
write.table(t(c(mse.ridge.X.vec, mse.krr.X.vec, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)




