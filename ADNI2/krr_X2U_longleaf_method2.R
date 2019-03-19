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

X = X[-c(770),]
label = label[-c(770)]
Y = Y[-c(770)]
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

threshold.vec = floor(sort(sapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list[[l]], threshold)$cor.vec), decreasing = T)*1000)/1000
threshold.vec = threshold.vec[threshold.vec>0.1]
res = cv.select_threshold.method2(threshold.vec, X.train.list, Y.train.list, nfolds = 10)
threshold = res$threshold

# K.list = lapply(1:n_label, function(l) getK(Y.train.list[[l]], X.train.list[[l]], threshold)$K)
# # X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]],plot = T))
# X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = K.list[[ix]], plot = F))
# 
# H.list = lapply(X2U.list, function(list) list$H)
# # P1.list = lapply(X2U.list, function(list) list$P1)
# P.list = lapply(X2U.list, function(list) list$P)
# # K.list = lapply(X2U.list, function(list) list$K)
# 
# U.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.train.list[[ix]])
# F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
# HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
# PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
# # Y_mean.train.list = lapply(1:n_label, function(ix) (P1.list[[ix]]%*%Y.train.list[[ix]])[1])
# 
# U.train = do.call(rbind, U.train.list)
# HY.train = do.call(c, HY.train.list)
# data.U.train = data.frame(Y = HY.train, U.train)

# ----------------------------- get F.test.list and U.test.list ------------------------------------
X.train.list = lapply(label.level, function(l) X.train[label.train == l,])
X.test.list = lapply(label.level, function(l) X.test[label.test == l,])
# X.combine.list = lapply(1:n_label, function(ix) rbind(X.train.list[[ix]], X.test.list[[ix]]))
# 
# X.combine.mean = lapply(X.combine.list, colMeans)
# X.combine.list = lapply(1:n_label, function(ix) sweep(X.combine.list[[ix]], 2, X.combine.mean[[ix]]))
# 
# X2U.list = lapply(1:n_label, function(ix) X2U2(X.combine.list[[ix]], K = K.list[[ix]], plot = F))
# H.list = lapply(X2U.list, function(list) list$H)
# # P1.list = lapply(X2U.list, function(list) list$P1)
# # K.list = lapply(X2U.list, function(list) list$K)
# P.list = lapply(X2U.list, function(list) list$P)
# 
# U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.combine.list[[ix]])
# F.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
# F.train.new.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][1:n.train.vec[ix],]))
# U.train.new.list = lapply(1:n_label, function(ix) U.list[[ix]][1:n.train.vec[ix],])
# 
# for (ix in 1:n_label){
#   F.list[[ix]][,-1] = t(t(F.list[[ix]][,-1])*sign(diag(as.matrix(cor(F.train.new.list[[ix]][,-1], F.train.list[[ix]][,-1])))))
# }
# U.test.list = lapply(1:n_label, function(ix) U.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],])
# F.test.list = lapply(1:n_label, function(ix) as.matrix(F.list[[ix]][(n.train.vec[ix]+1):n.vec[ix],]))
# # XF.test.list = lapply(1:n_label, function(ix) X.combine.list[[ix]][(n.train.vec[ix]+1): n.vec[ix],] - U.test.list[[ix]])
# # PYhat.ridge.XF.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.XF.list[[ix]], s=ml.ridge.XF.list[[ix]]$lambda.min, newx = XF.test.list[[ix]]))
# # -------------------------------------------------------------------------------------------------
# 
# # ------------------------------- compute weights from OLS.U --------------------------------------
# ml.lm.U = lm(Y~., data = data.U.train)  
# ix.vec = c(0,cumsum(n.train.vec))
# sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
# w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
# # w = rep(1, n.train)
# 
# # WLS.U: weights from OLS.U
# ml.lm.U = lm(Y~., data = data.U.train, weights = w)
# HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
# # res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# 
# PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
# HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
# Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
# mse.lm.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.lm.U1 = sum(mse.lm.U.vec1*n.test.vec)/sum(n.test.vec)
# 
# # weighted ridge: weights from OLS.U
# ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
# HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
# # res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = PY.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
# HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
# Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
# # Yhat.test.list = lapply(1:n_label, function(ix) PYhat.ridge.XF.test.list[[ix]] + HYhat.test.list[[ix]])
# mse.ridge.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.ridge.U1 = sum(mse.ridge.U.vec1*n.test.vec)/sum(n.test.vec)
# # -------------------------------------------------------------------------------------------------
# 
# # ----------------------- compute weights from OLS.U and OLS.F ------------------------------------
# ml.lm.U = lm(Y~., data = data.U.train)  
# HYhat.lm.U.train.list = lapply(1:n_label, function(l) ml.lm.U$fitted.values[(ix.vec[l]+1):ix.vec[l+1]])
# res.lm.U.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.lm.U.train.list[[l]])
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.lm.U.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] + ml.lm.F.list[[ix]]$fitted.values)
# sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
# w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
# # w = rep(1, n.train)
# # WLS.U: weights from OLS.U and OLS.F
# ml.lm.U = lm(Y~., data = data.U.train, weights = w)
# HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
# res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
# HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
# Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
# mse.lm.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.lm.U2 = sum(mse.lm.U.vec2*n.test.vec)/sum(n.test.vec)
# 
# # weighted ridge: weights from OLS.U and OLS.F
# ml.ridge.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
# HYhat.ridge.train.list = lapply(1:n_label, function(l) predict(ml.ridge.U, s = ml.ridge.U$lambda.min, newx = U.train)[(ix.vec[l]+1):ix.vec[l+1]])
# res.ridge.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.ridge.train.list[[l]])
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.ridge.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
# HYhat.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, U.test.list[[ix]]))
# Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
# # Yhat.test.list = lapply(1:n_label, function(ix) PYhat.ridge.XF.test.list[[ix]] + HYhat.test.list[[ix]])
# mse.ridge.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.ridge.U2 = sum(mse.ridge.U.vec2*n.test.vec)/sum(n.test.vec)
# 

Yhat.test.list = lm.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.lm.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U1 = sum(mse.lm.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = lm.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.lm.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.lm.U2 = sum(mse.lm.U.vec2*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.ridge.U.vec1 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U1 = sum(mse.ridge.U.vec1*n.test.vec)/sum(n.test.vec)

Yhat.test.list = ridge.U2.threshold.method2(X.train.list, Y.train.list, X.test.list, threshold)$Yhat.test.list
mse.ridge.U.vec2 = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
mse.ridge.U2 = sum(mse.ridge.U.vec2*n.test.vec)/sum(n.test.vec)

# -------------------------------------------------------------------------------------------------
# c(mse.lm.global, mse.ridge.global, mse.lm.WLS, mse.ridge.WLS, mse.ridge.X.class, mse.lm.U1, mse.ridge.U1, mse.lm.U2, mse.ridge.U2)
# rbind(mse.lm.global.vec, mse.ridge.global.vec, mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.lm.U.vec1, mse.ridge.U.vec1, mse.lm.U.vec2, mse.ridge.U.vec2)


# # lm and wls ridge
# ml.lm.U = lm(Y~., data = data.U.train)  
# Yhat.lm.U.train = lapply(1:n_label, function(ix) predict(ml.lm.U, newdata = data.train.list[[ix]][, -1])+Y_mean.train.list[[ix]])
# Yhat.lm.U.test = lapply(1:n_label, function(ix) predict(ml.lm.U, newdata = data.test.list[[ix]][, -1])+Y_mean.train.list[[ix]])
# mse.lm.U.vec = sapply(1:n_label, function(ix) mean((Yhat.lm.U.test[[ix]]-(Y.test.list[[ix]]))^2))
# mse.lm.U = sum(mse.lm.U.vec*n.test.vec)/sum(n.test.vec)
# 
# # weight computed from sx, ridge
# sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train[[l]])^2))
# w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l])))
# ml.ridge.U = cv.glmnet(x=U.train, y=HY.train, weights = w, alpha = 0)
# Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
# mse.ridge.sx.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.U.test[[ix]]+Y_mean.train.list[[ix]])-(Y.test.list[[ix]]))^2))
# mse.ridge.sx.U = sum(mse.ridge.sx.U.vec*n.test.vec)/sum(n.test.vec)
# 
# # weight computed from su, ridge
# ix.vec = c(0,cumsum(n.train.vec))
# sigma2 = sapply(1:n_label, function(ix) sum((ml.lm.U$residuals[(ix.vec[ix]+1):ix.vec[ix+1]])^2)/n.train.vec[ix])
# w = do.call(c, lapply(1:n_label, function(ix) rep(1/sigma2[ix], n.train.vec[ix])))
# ml.ridge.U = cv.glmnet(x=U.train, y=HY.train, weights = w, alpha = 0)
# Yhat.ridge.U.test = lapply(1:n_label, function(ix) predict(ml.ridge.U, s=ml.ridge.U$lambda.min, newx = X.test.list[[ix]]))
# mse.ridge.su.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.U.test[[ix]]+Y_mean.train.list[[ix]])-(Y.test.list[[ix]]))^2))
# mse.ridge.su.U = sum(mse.ridge.su.U.vec*n.test.vec)/sum(n.test.vec)
# 
# # WLS U lm
# ml.lm.U = lm(Y~., data = data.U.train, weights = w)
# Yhat.lm.sx.U.test = lapply(1:n_label, function(ix) predict(ml.lm.U, new = data.frame(X.test.list[[ix]])))
# mse.lm.su.U.vec = sapply(1:n_label, function(ix) mean(((Yhat.lm.sx.U.test[[ix]]+Y_mean.train.list[[ix]])-(Y.test.list[[ix]]))^2))
# mse.lm.su.U = sum(mse.lm.su.U.vec*n.test.vec)/sum(n.test.vec)
# 
# # 2 test set zheng ti
# HYhat.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]])
# res.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - HYhat.train.list[[l]])
# 
# data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = res.train.list[[ix]], F.train.list[[ix]][,-1]))
# ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))
# ml.ridge.F.list = lapply(1:n_label, function(ix) cv.glmnet(x = F.train.list[[ix]], y = res.train.list[[ix]], alpha = 0))
# # PYhat.train.list = lapply(1:n_label, function(ix) predict(ml.lm.F.list[[ix]], new = data.frame(F.train.list[[ix]][,-1])))
# 
# #######
# PYhat.WLSbeta.train.list = lapply(1:n_label, function(ix) F.train.list[[ix]]%*%X2U.list[[ix]]$L%*%ml.lm.WLS$coefficients[-1]+Y_mean.train.list[[ix]])
# HYhat.WLSbeta.train.list = lapply(1:n_label, function(ix) cbind(1,U.train.list[[ix]])%*%ml.lm.WLS$coefficients)
# #######
# 
# # X.train = X.list[[1]]
# # X.test = X.list[[2]]
# 
# 
# PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)
# PYhat.ridge.test.list = lapply(1:n_label, function(ix) predict(ml.ridge.F.list[[ix]], s = ml.ridge.F.list[[ix]]$lambda.min, newx = F.test.list[[ix]]))
# HYhat.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.U$coefficients)
# Yhat.test.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.list[[ix]])
# Yhat.ridgeF.test.list = lapply(1:n_label, function(ix) PYhat.ridge.test.list[[ix]] + HYhat.test.list[[ix]])
# mse.vec = sapply(1:n_label, function(ix) mean((Yhat.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.ridgeF.vec = sapply(1:n_label, function(ix) mean((Yhat.ridgeF.test.list[[ix]]-Y.test.list[[ix]])^2))
# mse.lm.WLS.vec
# 
# PYhat.WLSbeta.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%X2U.list[[ix]]$L%*%ml.lm.WLS$coefficients[-1]+Y_mean.train.list[[ix]])
# HYhat.WLSbeta.test.list = lapply(1:n_label, function(ix) cbind(1,U.test.list[[ix]])%*%ml.lm.WLS$coefficients)
# 
# diff.PYhat.vec = sapply(1:n_label, function(ix) mean((PYhat.WLSbeta.test.list[[ix]] - PYhat.test.list[[ix]])^2))
# diff.PYhat.ridge.vec = sapply(1:n_label, function(ix) mean((PYhat.WLSbeta.test.list[[ix]] - PYhat.ridge.test.list[[ix]])^2))
# diff.HYhat.vec = sapply(1:n_label, function(ix) mean((HYhat.WLSbeta.test.list[[ix]] - HYhat.test.list[[ix]])^2))
# 
# 
# 
# ##########

# PYhat.ridge.X.class.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%X2U.list[[ix]]$L%*%coef(ml.ridge.X.class[[ix]], s = ml.ridge.X.class[[ix]]$lambda.min)[-1]+Y_mean.train.list[[ix]])
# HYhat.ridge.X.class.test.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%coef(ml.ridge.X.class[[ix]], s = ml.ridge.X.class[[ix]]$lambda.min)[-1])
# 
# Yhat.ridge.X.class.test.list = lapply(1:n_label, function(ix) PYhat.ridge.X.class.test.list[[ix]] + HYhat.ridge.X.class.test.list[[ix]])
# mse.ridge.X.class.vec2 = sapply(1:n_label, function(ix) mean((Yhat.ridge.X.class.test.list[[ix]]-Y.test.list[[ix]])^2))
# 
# 
# PYhat.ridge.X.class.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%X2U.list[[ix]]$L%*%coef(ml.ridge.X.class[[ix]], s = ml.ridge.X.class[[ix]]$lambda.min)[-1]+Y_mean.train.list[[ix]])
# HYhat.ridge.WLS.test.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%coef(ml.ridge.WLS, s = ml.ridge.WLS$lambda.min)[-1])
# 
# Yhat.ridge.X.class.test.list = lapply(1:n_label, function(ix) PYhat.ridge.X.class.test.list[[ix]] + HYhat.ridge.WLS.test.list[[ix]])
# mse.ridge.X.class.vec2 = sapply(1:n_label, function(ix) mean((Yhat.ridge.X.class.test.list[[ix]]-Y.test.list[[ix]])^2))
# 
# 
# diff.PYhat.vec = sapply(1:n_label, function(ix) mean((PYhat.ridge.X.class.test.list[[ix]] - PYhat.test.list[[ix]])^2))
# diff.HYhat.vec = sapply(1:n_label, function(ix) mean((HYhat.ridge.X.class.test.list[[ix]] - HYhat.test.list[[ix]])^2))
# 
# cross.HYhatPYhat.vec = sapply(1:n_label, function(ix) mean(2*(HYhat.ridge.X.class.test.list[[ix]] - HYhat.test.list[[ix]])*(PYhat.ridge.X.class.test.list[[ix]] - PYhat.test.list[[ix]])))
# 
# diff.PYhat.vec + diff.HYhat.vec + cross.HYhatPYhat.vec

# file.name = paste0("result_overall_", threshold*100/5, ".csv")
file.name = "result_overall.csv"
write.table(t(c(mse.lm.global, mse.ridge.global, mse.lm.WLS, mse.ridge.WLS, mse.ridge.X.class, mse.lm.U1, mse.ridge.U1, mse.lm.U2, mse.ridge.U2,threshold, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

# file.name = paste0("result_class_", threshold*100/5, ".csv")
file.name = "result_class.csv"
write.table(t(c(mse.lm.global.vec, mse.ridge.global.vec, mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.ridge.X.class.vec, mse.lm.U.vec1, mse.ridge.U.vec1, mse.lm.U.vec2, mse.ridge.U.vec2, threshold, myseed)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)





