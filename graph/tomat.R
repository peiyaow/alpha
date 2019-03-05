n = dim(X)[1]
p = dim(X)[2]
X.mean = apply(X, 2, mean)
X.sd = apply(X, 2, sd)
X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
# label = as.factor(label)
label.level = levels(label)
n_label = length(label.level)

X.list = lapply(label.level, function(l) X[label == l,])
Y.list = lapply(label.level, function(l) Y[label == l])
n.vec = sapply(X.list, nrow)
label.list = lapply(1:n_label, function(l) rep(label.level[l], n.vec[l]))

# X = do.call(rbind, X.list)
# Y = do.call(c, Y.list)
# label = do.call(c, label.list)
# 
# X = X[-c(10),]
# Y = Y[-c(10)]
# label = as.factor(label[-c(10)])
# 
# save(X, Y, label, file = "ADNI2_clean5.RData")

X.mean.list = lapply(X.list, colMeans)
X.list = lapply(1:n_label, function(ix) sweep(X.list[[ix]], 2, X.mean.list[[ix]]))

X2U.list = lapply(X.list, function(X)  X2U1(X, plot = T))

mycut = X2U4(X.list, plot = T)
X2U.list = lapply(X.list, function(X)  X2U.cut(X, mycut))

H.list = lapply(X2U.list, function(list) list$H)
P1.list = lapply(X2U.list, function(list) list$P1)
K.list = lapply(X2U.list, function(list) list$K)

U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.list[[ix]])

newX = do.call(rbind, X.list)
U = do.call(rbind, U.list)

writeMat("data4.mat", X = newX, U=U, n_vec = n.vec)


