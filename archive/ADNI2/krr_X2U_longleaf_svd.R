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

# global EN
ml.EN.global = cv.glmnet(x=X.train, y= Y.train, alpha = 0.5)
Yhat.EN.global.test = predict(ml.EN.global, s=ml.EN.global$lambda.min, newx = X.test)
mse.EN.global.vec = sapply(label.level, function(l) mean((Yhat.EN.global.test[label.test==l] - Y.test[label.test==l])^2))
mse.EN.global = sum(mse.EN.global.vec*n.test.vec)/sum(n.test.vec)

# global lasso
ml.lasso.global = cv.glmnet(x=X.train, y= Y.train, alpha = 1)
Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
mse.lasso.global.vec = sapply(label.level, function(l) mean((Yhat.lasso.global.test[label.test==l] - Y.test[label.test==l])^2))
mse.lasso.global = sum(mse.lasso.global.vec*n.test.vec)/sum(n.test.vec)

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

# WLS EN
ml.EN.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 0.5, weights = w)
Yhat.EN.WLS.test = lapply(1:n_label, function(ix) predict(ml.EN.WLS, s=ml.EN.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.EN.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
mse.EN.WLS = sum(mse.EN.WLS.vec*n.test.vec)/sum(n.test.vec)

# WLS lasso
ml.lasso.WLS = cv.glmnet(x=X.train.WLS, y=Y.train.WLS, alpha = 1, weights = w)
Yhat.lasso.WLS.test = lapply(1:n_label, function(ix) predict(ml.lasso.WLS, s=ml.lasso.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.lasso.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
mse.lasso.WLS = sum(mse.lasso.WLS.vec*n.test.vec)/sum(n.test.vec)

# class ridge
ml.ridge.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0))
Yhat.ridge.X.class.test = lapply(1:n_label, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)

# class EN
ml.EN.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0.5))
Yhat.EN.X.class.test = lapply(1:n_label, function(ix) predict(ml.EN.X.class[[ix]], s=ml.EN.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.EN.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.EN.X.class = sum(mse.EN.X.class.vec*n.test.vec)/sum(n.test.vec)

# class lasso
ml.lasso.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 1))
Yhat.lasso.X.class.test = lapply(1:n_label, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.lasso.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)

# ------------------------------- ALPHA -----------------------------------------
# ----------------- perform the two-step screening to calculate K in each group ---------------------
X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))

# H.list = lapply(X2U.list, function(list) list$H)
# P.list = lapply(X2U.list, function(list) list$P)
K.list = lapply(X2U.list, function(list) list$K)
L.list =  lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) # loading matrix

F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
# PY.train.list = lapply(1:n_label, function(ix) P.list[[ix]]%*%Y.train.list[[ix]])
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                              F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))

# gamma under H1
coef.list = lapply(ml.lm.F.list, function(list) as.matrix(list$coefficients[-1])) 
# gamma under H0
coef.list2 = lapply(L.list, function(L) as.matrix(L%*%as.matrix(ml.lm.WLS$coefficients[-1]))) 
# concat two gammas
coef.list.all = lapply(1:n_label, function(ix) cbind(coef.list[[ix]], coef.list2[[ix]])) 

LXTXL.list = lapply(1:n_label, function(ix) L.list[[ix]]%*%solve(t(X.train.WLS)%*%diag(w)%*%X.train.WLS)%*%t(L.list[[ix]]))
coef.diff.list = lapply(1:n_label, function(ix) coef.list.all[[ix]][,2] - coef.list.all[[ix]][,1])

t.stat.list = lapply(1:n_label, function(ix) -abs(coef.diff.list[[ix]])/(sqrt(diag(LXTXL.list[[ix]]))*sigma(ml.lm.WLS)))
p_value.t.list = lapply(1:n_label, function(ix) 2*pt(t.stat.list[[ix]], df = ml.lm.WLS$df.residual))

f.stat.list = lapply(1:n_label, function(ix) sapply(1:K.list[[ix]], function(K) t(coef.diff.list[[ix]][1:K])%*%solve(LXTXL.list[[ix]])[1:K, 1:K]%*%(coef.diff.list[[ix]])[1:K]/(K*sigma(ml.lm.WLS)^2)))
p_value.f.list = lapply(1:n_label, function(ix) sapply(1:K.list[[ix]], function(K) pf(f.stat.list[[ix]][K], df1 = K, df2 = ml.lm.WLS$df.residual, lower.tail = F)))

# bonferroni K
K.list = lapply(1:n_label, function(ix) screenK.bonferroni(p_value.t.list[[ix]])) # K after screening

# Scheffe K forward
K.list = lapply(1:n_label, function(ix) screenK(p_value.f.list[[ix]])) # K after screening

# Scheffe K backward
K.list = lapply(1:n_label, function(ix) screenK(p_value.f.list[[ix]], forward = F)) # K after screening

# --------------
X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = K.list[[ix]], plot = F))
H.list = lapply(X2U.list, function(list) list$H)
P.list = lapply(X2U.list, function(list) list$P)
L.list = lapply(X2U.list, function(list) matrix(list$L[-1,], ncol = p)) 

F.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$F_)
U.train.list = lapply(1:n_label, function(ix) X2U.list[[ix]]$U)

FnU.test.list = lapply(1:n_label, function(ix) FnU.svd(X.test.list[[ix]], L.list[[ix]])) 
F.test.list = lapply(FnU.test.list, function(list) list$F_) 
U.test.list = lapply(FnU.test.list, function(list) list$U) 

# OLS.F
data.F.train.list = lapply(1:n_label, function(ix) data.frame(Y = Y.train.list[[ix]], 
                                                              F.train.list[[ix]][,-1]))
ml.lm.F.list = lapply(1:n_label, function(l) lm(Y~., data = data.F.train.list[[l]]))

# OLS.U
U.train = do.call(rbind, U.train.list)
HY.train.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%Y.train.list[[ix]])
HY.train = do.call(c, HY.train.list)
data.U.train = data.frame(Y = HY.train, U.train)
ml.lm.U = lm(Y~., data = data.U.train)

# compute weights from OLS.F and OLS.U 
Yhat.lm.U.train.list = lapply(1:n_label, function(ix) ml.lm.U$fitted.values[(ix.vec[ix]+1):ix.vec[ix+1]] 
                              + ml.lm.F.list[[ix]]$fitted.values)
sigma2 = sapply(1:n_label, function(l) mean((Y.train.list[[l]] - Yhat.lm.U.train.list[[l]])^2))
# new weights for U regression
w = do.call(c, lapply(1:n_label, function(l) rep(1/(sigma2[l]*(1-K.list[[l]]/n.train.vec[l])), n.train.vec[l]))) 

# WLS
OLS.WLS.U = lm(Y~., data = data.U.train, weight = w)
ridge.WLS.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0)
EN.WLS.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 0.5)
lasso.WLS.U = cv.glmnet(x = U.train, y = HY.train, weights = w, alpha = 1)

HYhat.test.OLS.WLS.list = lapply(1:n_label, function(ix) cbind(1, U.test.list[[ix]])%*%OLS.WLS.U$coefficients)
HYhat.test.ridge.WLS.list = lapply(1:n_label, function(ix) predict(ridge.WLS.U, s=ridge.WLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.EN.WLS.list = lapply(1:n_label, function(ix) predict(EN.WLS.U, s=EN.WLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.lasso.WLS.list = lapply(1:n_label, function(ix) predict(lasso.WLS.U, s=lasso.WLS.U$lambda.min, U.test.list[[ix]]))
PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)

Yhat.test.OLS.WLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.WLS.list[[ix]])
Yhat.test.ridge.WLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.WLS.list[[ix]])
Yhat.test.EN.WLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.WLS.list[[ix]])
Yhat.test.lasso.WLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.WLS.list[[ix]])

mse.OLS.WLS.list = compute.mse(Y.test.list, Yhat.test.OLS.WLS.list)
mse.ridge.WLS.list = compute.mse(Y.test.list, Yhat.test.ridge.WLS.list)
mse.EN.WLS.list = compute.mse(Y.test.list, Yhat.test.EN.WLS.list)
mse.lasso.WLS.list = compute.mse(Y.test.list, Yhat.test.lasso.WLS.list)

file.name = "result_overall.csv"
write.table(t(c(mse.lm.global, mse.ridge.global, mse.EN.global, mse.lasso.global, 
                mse.lm.WLS, mse.ridge.WLS, mse.EN.WLS, mse.lasso.WLS, 
                mse.ridge.X.class, mse.EN.X.class, mse.lasso.X.class, 
                mse.OLS.WLS.list[[2]], mse.ridge.WLS.list[[2]], mse.EN.WLS.list[[2]], mse.lasso.WLS.list[[2]], 
                myseed)),
            file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_class.csv"
write.table(t(c(mse.lm.global.vec, mse.ridge.global.vec, mse.EN.global.vec, mse.lasso.global.vec, 
                mse.lm.WLS.vec, mse.ridge.WLS.vec, mse.EN.WLS.vec, mse.lasso.WLS.vec, 
                mse.ridge.X.class.vec, mse.EN.X.class.vec, mse.lasso.X.class.vec, 
                mse.OLS.WLS.list[[1]], mse.ridge.WLS.list[[1]], mse.EN.WLS.list[[1]], mse.lasso.WLS.list[[1]], 
                myseed)), 
            file = file.name, sep = ',', append = T, col.names = F, row.names = F)





