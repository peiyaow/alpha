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

# X = X[label!=4,]
# Y = Y[label!=4]
# label = label[label!=4]
# label = droplevels(label)

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
X2U.list = lapply(X.train.list, function(X)  X2U2(X, plot = F))

H.list = lapply(X2U.list, function(list) list$H)
P1.list = lapply(X2U.list, function(list) list$P1)
K.list = lapply(X2U.list, function(list) list$K)
P.list = lapply(X2U.list, function(list) list$P)
U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
Y_mean.train.list = lapply(1:n_label, function(ix) (P1.list[[ix]]%*%Y.train.list[[ix]])[1])

U.train = do.call(rbind, U.train.list)
HY.train = do.call(c, HY.train.list)
data.U.train = data.frame(Y = HY.train, U.train)

F.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Fx)))
U.test.list = lapply(1:n_label, function(ix) t(apply(X.test.list[[ix]], 1, function(row) FnU(X.train.list[[ix]], row, F.train.list[[ix]], K.list[[ix]])$Ux)))

# compute weights from OLS.U
ml.lm.U = lm(Y~., data = data.U.train)  
ix.vec = c(0,cumsum(n.train.vec))
sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))

# WLS.U: weights from OLS.U
ml.lm.U = lm(Y~., data = data.U.train, weights = w)
HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))

PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
mse.lm.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U1 = sum(mse.lm.U.vec1*n.test.vec)/sum(n.test.vec)

# weighted ridge: weights from OLS.U
ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.ridge.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
mse.ridge.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U1 = sum(mse.ridge.U.vec1*n.test.vec)/sum(n.test.vec)

# PYhat.WLSbeta.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%X2U.list[[ix]]$L%*%ml.lm.WLS$coefficients[-1]+Y_mean.train.list[[ix]])
# HYhat.WLSbeta.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.WLS$coefficients)
# 
# diff.PYhat.vec = sapply(1:n_label, function(ix) mean((PYhat.WLSbeta.test.list[[ix]] - PYhat.test.list[[ix]])^2))
# diff.PYhat.ridge.vec = sapply(1:n_label, function(ix) mean((PYhat.WLSbeta.test.list[[ix]] - PYhat.ridge.test.list[[ix]])^2))
# diff.HYhat.vec = sapply(1:n_label, function(ix) mean((HYhat.WLSbeta.test.list[[ix]] - HYhat.test.list[[ix]])^2))

# --------------------------------------- new w lm and wls ridge -----------------------------------------------------
# compute weights from OLS.U and OLS.F
ml.lm.U = lm(Y~., data = data.U.train)  
HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))

# WLS.U: weights from OLS.U and OLS.F
ml.lm.U = lm(Y~., data = data.U.train, weights = w)
HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))

PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
mse.lm.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U2 = sum(mse.lm.U.vec2*n.test.vec)/sum(n.test.vec)

# weighted ridge: weights from OLS.U and OLS.F
ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.ridge.train.list[[ix]], F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
mse.ridge.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U2 = sum(mse.ridge.U.vec2*n.test.vec)/sum(n.test.vec)
# -------------------------------------------------------------------------------------------------------------------------------------

file.name = c("result_overall.csv")
write.table(t(c(mse.lm.global, mse.ridge.global, mse.lm.WLS, mse.ridge.WLS, mse.ridge.X.class, mse.lm.U1, mse.ridge.U1, mse.lm.U2, mse.ridge.U2,myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = c("result_class.csv")
write.table(t(c(mse.lm.global.vec, mse.ridge.global.vec, mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.lm.U.vec1, mse.ridge.U.vec1, mse.lm.U.vec2, mse.ridge.U.vec2, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
rbind(mse.lm.global.vec, mse.ridge.global.vec, mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.lm.U.vec1, mse.ridge.U.vec1, mse.lm.U.vec2, mse.ridge.U.vec2, myseed)



