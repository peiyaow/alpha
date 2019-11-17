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

#X.train = X.train[label.train < 4,]
#Y.train = Y.train[label.train < 4]
#X.test = X.test[label.test < 4,]
#Y.test = Y.test[label.test < 4]
#label.train = droplevels(label.train[label.train < 4])
#label.test = droplevels(label.test[label.test < 4])

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
X.train = do.call(rbind, X.train.list)
data.train.WLS = data.frame(Y = Y.train.WLS, X.train)
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
ml.ridge.WLS = cv.glmnet(x=X.train, y=Y.train.WLS, alpha = 0, weights = w)
Yhat.ridge.WLS.test = lapply(1:n_label, function(ix) predict(ml.ridge.WLS, s=ml.ridge.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.ridge.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.ridge.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
mse.ridge.WLS = sum(mse.ridge.WLS.vec*n.test.vec)/sum(n.test.vec)

# WLS EN
ml.EN.WLS = cv.glmnet(x=X.train, y=Y.train.WLS, alpha = 0.5, weights = w)
Yhat.EN.WLS.test = lapply(1:n_label, function(ix) predict(ml.EN.WLS, s=ml.EN.WLS$lambda.min, newx = X.test.list[[ix]]))
mse.EN.WLS.vec = sapply(1:n_label, function(ix) mean(((Yhat.EN.WLS.test[[ix]]+Y.train.mean[[ix]])-(Y.test.list[[ix]]))^2))
mse.EN.WLS = sum(mse.EN.WLS.vec*n.test.vec)/sum(n.test.vec)

# WLS lasso
ml.lasso.WLS = cv.glmnet(x=X.train, y=Y.train.WLS, alpha = 1, weights = w)
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




jive.list = myjive(X.train.list)
L = jive.list$L
K = jive.list$K
F.train.list = jive.list$F.list
U.train.list = jive.list$U.list
R.train.list = jive.list$R.list
P.list = jive.list$P.list






PY.list = lapply(1:n_label, function(l) P.list[[l]]%*%Y.train.list[[l]])
HY.train.list = lapply(1:n_label, function(l) Y.train.list[[l]] - PY.list[[l]])

FnU.test.list = lapply(1:n_label, function(ix) FnU.svd2(X.test.list[[ix]], L)) 
F.test.list = lapply(FnU.test.list, function(list) list$F) 
U.test.list = lapply(FnU.test.list, function(list) list$U) 

ncum.test.vec = c(0, cumsum(n.test.vec))
ncum.train.vec = c(0, cumsum(n.train.vec))
hehe = lm(Y.train.WLS~do.call(rbind, F.train.list))

Yhat = cbind(1, do.call(rbind, F.train.list))%*%as.matrix(hehe$coefficients)
Yhat.list = lapply(1:(n_label), function(l) Yhat[(ncum.train.vec[l]+1):ncum.train.vec[l+1],1])
mse.new.lm.WLS.vec1 = sapply(1:n_label, function(ix) mean((Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.train.list[[ix]])^2))

ww = do.call(c, sapply(1:n_label, function(ix) rep(1/mse.new.lm.WLS.vec1[ix], n.train.vec[ix])))
hehe.w = lm(Y.train.WLS~do.call(rbind, F.train.list), weights = ww)

Yhat = cbind(1, do.call(rbind, F.train.list))%*%as.matrix(hehe.w$coefficients)
Yhat.list = lapply(1:(n_label), function(l) Yhat[(ncum.train.vec[l]+1):ncum.train.vec[l+1],1])
mse.train.vec = sapply(1:n_label, function(ix) mean((Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.train.list[[ix]])^2))
HY.train.list = lapply(1:n_label, function(ix) Y.train.list[[ix]]-Yhat.list[[ix]]-Y.train.mean[[ix]])

###
FL.train = do.call(rbind, F.train.list)%*%L
ml.FL = cv.glmnet(x = FL.train, y = Y.train.WLS, alpha = 0, weights = ww)
###

FU.list = lapply(1:n_label, function(l) as.matrix(commonFL_given_K(U.train.list[[l]])$F))
KU.list = lapply(1:n_label, function(l) ncol(FU.list[[l]]))
L_U.list = lapply(1:n_label, function(l) as.matrix(commonFL_given_K(U.train.list[[l]])$L))

KU1.list = list()
for (l in 1:n_label){
  K = KU.list[[l]]
  for (k in 1:K){
    ml = lm(HY.train.list[[l]]~FU.list[[l]][,1:k])
    pvalue = summary(ml)$coefficients[-1,4][k]
    print(pvalue)
    if (pvalue >= 0.05){break}
  }
  KU1.list[[l]] = k-1
}


# testing
Yhat = cbind(1, do.call(rbind, F.test.list))%*%as.matrix(hehe.w$coefficients)
Yhat.list = lapply(1:(n_label), function(l) Yhat[(ncum.test.vec[l]+1):ncum.test.vec[l+1],1])
mse.new.lm.WLS.vec1 = sapply(1:n_label, function(ix) mean((Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.lm.WLS.vec
mse.new.lm.WLS.vec1

# version 3
lm.FU.list = list()
FU1.test.list = list()
for (k in 1:n_label){
  if (KU1.list[[k]] == 0){
    lm.FU.list[[k]] = lm(HY.train.list[[k]]~NULL)
    FU = NULL
  }else{
    lm.FU.list[[k]] = lm(HY.train.list[[k]]~FU.list[[k]][,1:KU1.list[[k]]])
    FU = FnU.svd2(U.test.list[[k]], L_U.list[[k]][1:KU1.list[[k]],])$F
  }
  FU1.test.list[[k]]=cbind(matrix(1, nrow = n.test.vec[k]), FU)
  
}

# lmlm
HYhat.list = lapply(1:n_label, function(l) FU1.test.list[[l]]%*%as.matrix(lm.FU.list[[l]]$coefficients))
mse.new.lm.WLS.vec03 = sapply(1:n_label, function(ix) mean((HYhat.list[[ix]] + Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))
mse.new.lm.WLS03 = sum(mse.new.lm.WLS.vec03*n.test.vec)/sum(n.test.vec)
mse.new.lm.WLS.vec1
mse.new.lm.WLS.vec03

Yhat.list1 = Yhat.list

# ridgelm
FL.test = do.call(rbind, F.test.list)%*%L
Yhat = predict(ml.FL, newx = FL.test, s = ml.FL$lambda.min)
Yhat.list = lapply(1:(n_label), function(l) Yhat[(ncum.test.vec[l]+1):ncum.test.vec[l+1],1])
mse.new.ridgelm.WLS03 = sapply(1:n_label, function(ix) mean((HYhat.list[[ix]] + Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))

# ridgelasso
ml.U.list = lapply(1:n_label, function(l) cv.glmnet(y = HY.train.list[[l]], x = U.train.list[[l]], a = 1))
HYhat.list = lapply(1:n_label, function(l) predict(ml.U.list[[l]], newx = U.test.list[[l]], s = ml.U.list[[l]]$lambda.min))
mse.new.ridgelasso.WLS03 = sapply(1:n_label, function(ix) mean((HYhat.list[[ix]] + Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))

# lmlasso
Yhat.list = Yhat.list1
mse.new.lmlasso.WLS03 = sapply(1:n_label, function(ix) mean((HYhat.list[[ix]] + Yhat.list[[ix]]+Y.train.mean[[ix]]-Y.test.list[[ix]])^2))

mse.new.lm.WLS.vec1
mse.ridge.WLS.vec
mse.new.lm.WLS.vec03
mse.new.ridgelm.WLS03
mse.new.ridgelasso.WLS03
mse.new.lmlasso.WLS03

