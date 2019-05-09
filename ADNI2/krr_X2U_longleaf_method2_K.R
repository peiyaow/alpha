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
load("/nas/longleaf/home/peiyao/alpha/data/ADNI2_clean3.RData")

# X = X[-c(770),]
# label = label[-c(770)]
# Y = Y[-c(770)]
# 
# X = X[label!=4,]
# Y = Y[label!=4]
# label = label[label!=4]
# label = droplevels(label)
# 
# X = X[label!=3,]
# Y = Y[label!=3]
# label = label[label!=3]
# label = droplevels(label)
# 
# X = X[Y < 37,]
# label = label[Y<37]
# Y = Y[Y<37]

n = dim(X)[1]
p = dim(X)[2]

set.seed(myseed)

ix.train = unlist(createDataPartition(label, times = 1, p = 3/4))
# ix.train = sample(n, floor(n*3/4))
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

# data.frame format
data.train.list = lapply(1:n_label, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
data.test.list = lapply(1:n_label, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))

# global lm
data.train = data.frame(Y=Y.train, X.train)
ml.lm.global = lm(Y~., data = data.train)
Yhat.lm.global.test = predict(ml.lm.global, new = data.frame(X.test))
mse.lm.global.vec = sapply(label.level, function(l) mean((Yhat.lm.global.test[label.test==l] - Y.test[label.test==l])^2))
mse.lm.global = sum(mse.lm.global.vec*n.test.vec)/sum(n.test.vec)

# global ridge
ml.ridge.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0)
Yhat.ridge.global.test = predict(ml.ridge.global, s=ml.ridge.global$lambda.min, newx = X.test)
mse.ridge.global.vec = sapply(label.level, function(l) mean((Yhat.ridge.global.test[label.test==l] - Y.test[label.test==l])^2))
mse.ridge.global = sum(mse.ridge.global.vec*n.test.vec)/sum(n.test.vec)

# WLS
Y.train.mean = lapply(Y.train.list, mean)
Y.train.WLS = do.call(c, lapply(1:n_label, function(l) Y.train.list[[l]] - Y.train.mean[[l]]))
X.train.WLS = do.call(rbind, X.train.list)
data.train.WLS = data.frame(Y = Y.train.WLS, X.train.WLS)
ml.lm.WLS = lm(Y~., data = data.train.WLS)
ix.vec = c(0,cumsum(n.train.vec))
sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.WLS$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
# w = rep(1, n.train)
# WLS lm
ml.lm.WLS = lm(Y~., data = data.train.WLS, weights = w)
Yhat.lm.WLS.test = lapply(1:n_label, function(ix) predict(ml.lm.WLS, new = data.frame(X.test.list[[ix]])))
mse.lm.WLS.vec = sapply(1:n_label, function(ix) mean((Yhat.lm.WLS.test[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS = sum(mse.lm.WLS.vec*n.test.vec)/sum(n.test.vec)

# WLS ridge
ml.ridge.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 0, weights = w)
Yhat.ridge.WLS.test = lapply(1:n_label, function(ix) predict(ml.ridge.WLS, s=ml.ridge.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.WLS = sum(mse.ridge.WLS.vec*n.test.vec)/sum(n.test.vec)

# class ridge
ml.ridge.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
Yhat.ridge.X.class.test = lapply(1:n_label, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)

# ------------------------------- ALPHA -----------------------------------------
# calculate K
X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))

H.list = lapply(X2U.list, function(list) list$H)
P.list = lapply(X2U.list, function(list) list$P)
K.list = lapply(X2U.list, function(list) list$K)

F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])

data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
coef.list = lapply(ml.lm.F.list, function(list) as.matrix(list$coefficients[-1]))

L.list =  lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p))
#L.list2 =  lapply(X2U.list, function(list) list$L[-1,])
coef.list2 = lapply(L.list, function(L) as.matrix(L%*%as.matrix(ml.lm.WLS$coefficients[-1])))
coef.list.all = lapply(1:n_label, function(ix) cbind(coef.list[[ix]], coef.list2[[ix]]))

p_value.list = lapply(1:n_label, function(ix) pt(-abs(coef.list.all[[ix]][,2] - coef.list.all[[ix]][,1] )/(sqrt(diag(L.list[[ix]]%*%solve(t(X.train.WLS)%*%diag(w)%*%X.train.WLS)%*%t(L.list[[ix]])))*sigma(ml.lm.WLS)), df = ml.lm.WLS$df.residual))

K.list = lapply(1:n_label, function(ix) screenK(p_value.list[[ix]], forward = T))
# --------------

## calculate threshold
# Y_concat = do.call(c, Y.train.list)
# X_concat = do.call(rbind, X.train.list)
# threshold = sort(abs(cor(Y_concat, X_concat)), decreasing = T)[1]

X.train.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])

Yhat.test.list = lm.U.K.method2(X.train.list, Y.train.list, X.test.list, K.list)$Yhat.test.list
mse.lm.U.vec10 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U10 = sum(mse.lm.U.vec10*n.test.vec)/sum(n.test.vec)

# Yhat.test.list = lm.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
# mse.lm.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.lm.U1 = sum(mse.lm.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = lm.U2.K.method2(X.train.list, Y.train.list, X.test.list, K.list)$Yhat.test.list
mse.lm.U.vec20 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U20 = sum(mse.lm.U.vec20*n.test.vec)/sum(n.test.vec)

# Yhat.test.list = lm.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
# mse.lm.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.lm.U2 = sum(mse.lm.U.vec2*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U.K.method2(X.train.list, Y.train.list, X.test.list, K.list)$Yhat.test.list
mse.ridge.U.vec10 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U10 = sum(mse.ridge.U.vec10*n.test.vec)/sum(n.test.vec)

# Yhat.test.list = ridge.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
# mse.ridge.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.ridge.U1 = sum(mse.ridge.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U2.K.method2(X.train.list, Y.train.list, X.test.list, K.list)$Yhat.test.list
mse.ridge.U.vec20 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U20 = sum(mse.ridge.U.vec20*n.test.vec)/sum(n.test.vec)

# Yhat.test.list = ridge.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
# mse.ridge.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.ridge.U2 = sum(mse.ridge.U.vec2*n.test.vec)/sum(n.test.vec)


file.name = "result_overall.csv"
write.table(t(c(mse.lm.global, mse.ridge.global, mse.lm.WLS, mse.ridge.WLS, mse.ridge.X.class, mse.lm.U10, mse.ridge.U10, mse.lm.U20, mse.ridge.U20, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

# file.name = paste0("result_class_", threshold*100/5, ".csv")
file.name = "result_class.csv"
write.table(t(c(mse.lm.global.vec, mse.ridge.global.vec, mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.lm.U.vec10, mse.ridge.U.vec10, mse.lm.U.vec20, mse.ridge.U.vec20, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)





