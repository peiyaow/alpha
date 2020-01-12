library(readr)
library(caret)
library(glmnet)

t = 1
X0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/X1a.txt", col_names = FALSE))
Y0 = as.matrix(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/Y5T.txt", col_names = FALSE))[, 4+5*(t-1)]
label0 = as.ordered(read_table("/Users/MonicaW/Documents/GitHub/LWPR/data/label1.txt", col_names = FALSE)$X1)
source('~/Documents/Research/coding/R/alpha/functions.R')

# remove NAs in Y
id.cut1 = !is.na(Y0)
X0 = X0[id.cut1, ]
Y0 = Y0[id.cut1]
label0 = label0[id.cut1]

# remove negatives in Y
id.cut2 = (Y0 > 0)
X0 = X0[id.cut2,]
Y0 = Y0[id.cut2]
label0 = label0[id.cut2]

# only select MRI+PET
X0 = X0[, 1:93]
X = X0
Y = Y0
label = label0

n = dim(X)[1]
p = dim(X)[2]

set.seed(1023)
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

X.train.mean = apply(X.train, 2, mean)
X.train.sd =  apply(X.train, 2, sd)
X.train = sweep(sweep(X.train, 2, X.train.mean), 2, X.train.sd, "/")
X.test = sweep(sweep(X.test, 2, X.train.mean), 2, X.train.sd, "/")

# overall model
# lm
data.train = data.frame(Y = Y.train, X.train)
data.test = data.frame(Y = Y.test, X.test)
ml.X = lm(Y~., data = data.train)
Yhat.X.test = predict(ml.X, newdata = data.test[, -1])
mse.X = mean((Yhat.X.test - Y.test)^2)

# lasso
ml.lasso.X = cv.glmnet(x=X.train, y=Y.train)
Yhat.X.lasso.test = predict(ml.lasso.X, s = ml.lasso.X$lambda.min, newx = X.test)
mse.lasso.X = mean((Yhat.X.lasso.test - Y.test)^2)

# seperate train
# lm
X.train.list = lapply(c(4,3,0), function(l) X.train[label.train == l,])
X.test.list = lapply(c(4,3,0), function(l) X.test[label.test == l,])
Y.train.list = lapply(c(4,3,0), function(l) Y.train[label.train == l])
Y.test.list = lapply(c(4,3,0), function(l) Y.test[label.test == l])
n.train.vec = sapply(X.train.list, nrow)
n.test.vec = sapply(X.test.list, nrow)
# standardize Y
# Y.train.mean = lapply(Y.train.list, mean)
# Y.train.sd = lapply(Y.train.list, sd)
# Y.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]])
# plot(do.call(c, Y.train.list))

# standardize X
# X.train.mean = lapply(X.train.list, colMeans)
# X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
# X.train.list = lapply(1:3, function(ix) sweep(sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
# X.test.list = lapply(1:3, function(ix) sweep(sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))

data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))

# lm
ml.X.class = lapply(1:3, function(ix) lm(Y~., data=data.train.list[[ix]]))
Yhat.X.class.test = lapply(1:3, function(ix) predict(ml.X.class[[ix]], newdata = data.test.list[[ix]][, -1]))
mse.X.class.vec = sapply(1:3, function(ix) mean((Yhat.X.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.X.class = sum(mse.X.class.vec*n.test.vec)/sum(n.test.vec)

# lasso
ml.lasso.X.class = lapply(1:3, function(ix) cv.glmnet(x=X.train.list[[ix]], y=Y.train.list[[ix]]))
Yhat.lasso.X.class.test = lapply(1:3, function(ix) predict(ml.lasso.X.class[[ix]], s=ml.lasso.X.class[[ix]]$lambda.min, newx = X.test.list[[ix]]))
mse.lasso.X.class.vec = sapply(1:3, function(ix) mean((Yhat.lasso.X.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.lasso.X.class = sum(mse.lasso.X.class.vec*n.test.vec)/sum(n.test.vec)

# ALPHA
# standardize Y
Y.train.mean = lapply(Y.train.list, mean)
Y.train.sd = lapply(Y.train.list, sd)
Y.train.list = lapply(1:3, function(ix) Y.train.list[[ix]] - Y.train.mean[[ix]])
Y.test.list = lapply(1:3, function(ix) Y.test.list[[ix]] - Y.train.mean[[ix]])

# standardize X
X.train.mean = lapply(X.train.list, colMeans)
# X.train.sd = lapply(X.train.list, function(X) apply(X, 2, sd))
X.train.list = lapply(1:3, function(ix) sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]))
X.test.list = lapply(1:3, function(ix) sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]))

data.train.list = lapply(1:3, function(ix) data.frame(Y=Y.train.list[[ix]], X.train.list[[ix]]))
data.test.list = lapply(1:3, function(ix) data.frame(Y=Y.test.list[[ix]], X.test.list[[ix]]))

# X.train.list = lapply(1:3, function(ix) sweep(sweep(X.train.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))
# X.test.list = lapply(1:3, function(ix) sweep(sweep(X.test.list[[ix]], 2, X.train.mean[[ix]]), 2, X.train.sd[[ix]], "/"))

H1 = X2U(X.train.list[[1]], plot = F)$H
H2 = X2U(X.train.list[[2]], plot = F)$H
H3 = X2U(X.train.list[[3]], plot = F)$H

U1 = H1%*%X.train.list[[1]]
U2 = H2%*%X.train.list[[2]]
U3 = H3%*%X.train.list[[3]]

Y_.1 = H1%*%Y.train.list[[1]]
Y_.2 = H2%*%Y.train.list[[2]]
Y_.3 = H3%*%Y.train.list[[3]]

PY.1 = Y.train.list[[1]] - Y_.1
PY.2 = Y.train.list[[2]] - Y_.2
PY.3 = Y.train.list[[3]] - Y_.3
  
U.train = rbind(U1, U2, U3)
Y_.train = c(Y_.1, Y_.2, Y_.3)
PY.train = c(PY.1, PY.2, PY.3)

plot(Y_.train)
plot(PY.train)  
  
data.U.train = data.frame(Y = Y_.train, U.train)
ml.U = lm(Y~., data = data.U.train)  
Yhat.U.class.test = lapply(1:3, function(ix) predict(ml.U, newdata = data.test.list[[ix]][, -1]))
mse.U.class.vec = sapply(1:3, function(ix) mean((Yhat.U.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.U.class = sum(mse.U.class.vec*n.test.vec)/sum(n.test.vec)

ml.lasso.U = cv.glmnet(x=U.train, y=Y_.train)
Yhat.lasso.U.class.test = lapply(1:3, function(ix) predict(ml.lasso.U, s=ml.lasso.U$lambda.min, newx = X.test.list[[ix]]))
mse.lasso.U.class.vec = sapply(1:3, function(ix) mean((Yhat.lasso.U.class.test[[ix]]-Y.test.list[[ix]])^2))
mse.lasso.U.class = sum(mse.lasso.U.class.vec*n.test.vec)/sum(n.test.vec)

#WLS
# n.train.vec = sapply(X.train.list, nrow)
# w = 1/c(rep(Y.train.sd[[1]], n.train.vec[[1]]), rep(Y.train.sd[[2]], n.train.vec[[2]]), rep(Y.train.sd[[3]], n.train.vec[[3]]))
# ml.U = lm(Y~., data = data.U.train, weights = w)  
# Yhat.U.class.test = lapply(1:3, function(ix) predict(ml.U, newdata = data.test.list[[ix]][, -1]))
# mse.U.class.vec = sapply(1:3, function(ix) mean((Yhat.U.class.test[[ix]]-Y.test.list[[ix]])^2))
# mse.U.class = sum(mse.U.class.vec*n.test.vec)/sum(n.test.vec)
# 
# 
# ml.lasso.U = cv.glmnet(x=U.train, y=Y_.train, weights = w)
# Yhat.lasso.U.class.test = lapply(1:3, function(ix) predict(ml.lasso.U, s=ml.lasso.U$lambda.min, newx = X.test.list[[ix]]))
# mse.lasso.U.class.vec = sapply(1:3, function(ix) mean((Yhat.lasso.U.class.test[[ix]]-Y.test.list[[ix]])^2))
# mse.lasso.U.class = sum(mse.lasso.U.class.vec*n.test.vec)/sum(n.test.vec)

  
  