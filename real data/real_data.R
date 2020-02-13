# ---------------------- reading shell command --------------------- 
args = (commandArgs(TRUE))
cat(args, "\n")
for (i in 1:length(args)) {
 eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------ 

library(caret)
library(glmnet)
require(methods)

path = "/nas/longleaf/home/peiyao/alpha/"

source(paste0(path, "/function/main_function.R"))
source(paste0(path, "/function/sim_function.R"))
load(paste0(path, "/data/ADNI2_clean3.RData"))

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

# # global lm
# data.train = data.frame(Y=Y.train, X.train)
# ml.lm.global = lm(Y~., data = data.train)
# Yhat.lm.global.test = predict(ml.lm.global, new = data.frame(X.test))
# mse.lm.global.vec = sapply(label.level, function(l) mean((Yhat.lm.global.test[label.test==l] - Y.test[label.test==l])^2))
# mse.lm.global = sum(mse.lm.global.vec*n.test.vec)/sum(n.test.vec)

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

# class ridge
ml.ridge.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], alpha = 0, standardize = F, lambda = exp(log(seq(p, .1, length.out = 50)))))
Yhat.ridge.X.class.test = lapply(1:n_label, function(ix) predict(ml.ridge.X.class[[ix]], s=ml.ridge.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.X.class = sum(mse.ridge.X.class.vec*n.test.vec)/sum(n.test.vec)

# class EN
ml.EN.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], standardize = F, alpha = 0.5))
Yhat.EN.X.class.test = lapply(1:n_label, function(ix) predict(ml.EN.X.class[[ix]], s=ml.EN.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.EN.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.EN.X.class = sum(mse.EN.X.class.vec*n.test.vec)/sum(n.test.vec)

# class lasso
ml.lasso.X.class = lapply(1:n_label, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]], standardize = F, alpha = 1))
Yhat.lasso.X.class.test = lapply(1:n_label, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.lasso.X.class.vec = sapply(1:n_label, function(ix) mean(((Yhat.lasso.X.class.test[[ix]])-(Y.test.list[[ix]]))^2))
mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)

# ------------------------------- ALPHA -----------------------------------------
X2U.list = lapply(1:n_label, function(ix) X2U1(X.train.list[[ix]], plot = F))
H.list = lapply(X2U.list, function(list) list$H)
K.list = lapply(X2U.list, function(list) list$K)
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

ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = F)
EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5, standardize = F)
lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = F)

HYhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.EN.OLS.list = lapply(1:n_label, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)

# Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
Yhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
Yhat.test.EN.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])

mse.ridge.OLS.list = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
mse.EN.OLS.list = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
mse.lasso.OLS.list = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)

#------------------------------------------ALPHA-0----------------------------------------------------
X2U.list = lapply(1:n_label, function(ix) X2U2(X.train.list[[ix]], K = 0, plot = F))
H.list = lapply(X2U.list, function(list) list$H)
K.list = lapply(X2U.list, function(list) list$K)
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

ridge.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0, standardize = F)
EN.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 0.5, standardize = F)
lasso.OLS.U = cv.glmnet(x = U.train, y = HY.train, alpha = 1, standardize = F)

# HYhat.test.OLS.list = lapply(1:n_label, function(ix) U.test.list[[ix]]%*%beta.OLS.U)
HYhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) predict(ridge.OLS.U, s=ridge.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.EN.OLS.list = lapply(1:n_label, function(ix) predict(EN.OLS.U, s=EN.OLS.U$lambda.min, U.test.list[[ix]]))
HYhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) predict(lasso.OLS.U, s=lasso.OLS.U$lambda.min, U.test.list[[ix]]))
PYhat.test.list = lapply(1:n_label, function(ix) F.test.list[[ix]]%*%ml.lm.F.list[[ix]]$coefficients)

# Yhat.test.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.OLS.list[[ix]])
Yhat.test.ridge.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.ridge.OLS.list[[ix]])
Yhat.test.EN.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.EN.OLS.list[[ix]])
Yhat.test.lasso.OLS.list = lapply(1:n_label, function(ix) PYhat.test.list[[ix]] + HYhat.test.lasso.OLS.list[[ix]])

# mse.OLS.list = compute.mse(Y.test.list, Yhat.test.OLS.list)
mse.ridge.OLS.list0 = compute.mse(Y.test.list, Yhat.test.ridge.OLS.list)
mse.EN.OLS.list0 = compute.mse(Y.test.list, Yhat.test.EN.OLS.list)
mse.lasso.OLS.list0 = compute.mse(Y.test.list, Yhat.test.lasso.OLS.list)

print(c(mse.ridge.global, mse.ridge.X.class, mse.ridge.OLS.list0$mse, mse.ridge.OLS.list$mse))
print(c(mse.EN.global, mse.EN.X.class, mse.EN.OLS.list0$mse, mse.EN.OLS.list$mse))
print(c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list0$mse, mse.lasso.OLS.list$mse))

file.name = "result_overall_ridge.csv"
write.table(t(c(mse.ridge.global, mse.ridge.X.class, mse.ridge.OLS.list0$mse, mse.ridge.OLS.list$mse)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_overall_EN.csv"
write.table(t(c(mse.EN.global, mse.EN.X.class, mse.EN.OLS.list0$mse, mse.EN.OLS.list$mse)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_overall_lasso.csv"
write.table(t(c(mse.lasso.global, mse.lasso.X.class, mse.lasso.OLS.list0$mse, mse.lasso.OLS.list$mse)), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_class_ridge.csv"
write.table(t(c(mse.ridge.global.vec, mse.ridge.X.class.vec, mse.ridge.OLS.list0[[1]], mse.ridge.OLS.list[[1]])), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_class_EN.csv"
write.table(t(c(mse.EN.global.vec, mse.EN.X.class.vec, mse.EN.OLS.list0[[1]], mse.EN.OLS.list[[1]])), file = file.name, sep = ',', append = T, col.names = F, row.names = F)

file.name = "result_class_lasso.csv"
write.table(t(c(mse.lasso.global.vec, mse.lasso.X.class.vec, mse.lasso.OLS.list0[[1]], mse.lasso.OLS.list[[1]])), file = file.name, sep = ',', append = T, col.names = F, row.names = F)




