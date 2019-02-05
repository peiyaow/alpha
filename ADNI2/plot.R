X1 = X[label == 3, ]
n = nrow(X1)
#  p = ncol(X)
PCA.res = eigen(X1%*%t(X1)/n)
plot(PCA.res$vectors[,1]*sqrt(n), PCA.res$vectors[,2]*sqrt(n))

plot(Y.train.list[[4]])

X2U.list = lapply(X.train.list, function(X)  X2U(X, plot = T))
F1 = c()
F2 = c()
for (l in 1:4){
  F1 = c(F1, X2U.list[[l]]$F1)
  F2 = c(F2, X2U.list[[l]]$F2)
  plot(X2U.list[[l]]$F1,  X2U.list[[l]]$F2)
}
par(mfrow = c(1,1))
plot(F1, F2)

L = NULL
L_res = NULL
for (l in 1:4){
  L = cbind(L, X2U.list[[l]]$L2)
  L_res = cbind(L_res, X2U.list[[l]]$L_res[3:4,])
  #plot(X2U.list[[l]]$L2[1,], X2U.list[[l]]$L2[2,])
}

#plot(L[1,], L[2,])
#label = c(rep("NC", 93), rep("SMC", 93), rep("eMCI", 93), rep("lMCI", 93))


data = data.frame(L1 = L[1,], L2 = L[2,], Lr1 = L_res[1,], Lr2 = L_res[2,], stage = c(rep("NC", 93), rep("SMC", 93), rep("eMCI", 93), rep("lMCI", 93)))
rownames(data) = NULL
library(ggplot2)
ggplot(data, aes(x=L1, y=L2, color = stage)) + geom_point()
ggplot(data, aes(x=Lr1, y=Lr2, color = stage)) + geom_point()

PY.list = list()
HY.list = list()
for (l in 1:4){
  PY.list[[l]] = X2U1(X.train.list[[l]])$P%*%Y.train.list[[l]]
  HY.list[[l]] = X2U1(X.train.list[[l]])$H%*%Y.train.list[[l]]
}

plot(do.call(c, PY.list))
plot(do.call(c, HY.list))


