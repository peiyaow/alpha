# n = dim(X)[1]
# p = dim(X)[2]
# X.mean = apply(X, 2, mean)
# X.sd = apply(X, 2, sd)
# X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
# # label = as.factor(label)
# label.level = levels(label)
# n_label = length(label.level)
# 
# X.list = lapply(label.level, function(l) X[label == l,])
# Y.list = lapply(label.level, function(l) Y[label == l])
# n.vec = sapply(X.list, nrow)
# label.list = lapply(1:n_label, function(l) rep(label.level[l], n.vec[l]))
# 
# X.mean.list = lapply(X.list, colMeans)
# X.list = lapply(1:n_label, function(ix) sweep(X.list[[ix]], 2, X.mean.list[[ix]]))
# 
# X2U.list = lapply(X.list, function(X)  X2U1(X, plot = T))
# 
# mycut = X2U4(X.list, plot = T)
# X2U.list = lapply(X.list, function(X)  X2U.cut(X, mycut))
# 
# H.list = lapply(X2U.list, function(list) list$H)
# P1.list = lapply(X2U.list, function(list) list$P1)
# K.list = lapply(X2U.list, function(list) list$K)
# 
# U.list = lapply(1:n_label, function(ix) H.list[[ix]]%*%X.list[[ix]])
library(R.matlab)

newX = do.call(rbind, X.list)
U = do.call(rbind, U.list)
FL = do.call(rbind, FL.list)

# apply(cbind(apply(newX, 1, min), apply(U,1,min), apply(FL,1, min)), 1, min)
# 
# writeMat("data5.mat", X = newX, U=U, FL = FL, n_vec = n.vec)
# order(apply(exp(FL),1, max)[1:173])
# # 19  85  44  18 162 122  78
# 570 + order(apply(exp(newX),1, max)[571:671])
# # 574 578 622 606 602 632 667 614
# 
# max(exp(newX[c(19,85,574,578),]))
# max(exp(FL[c(19,85,574,578),]))
# max(exp(U[c(19,85,574,578),]))

order(apply(exp(newX), 1, max)[1:173])
apply(exp(newX), 1, max)[1:173]
# 106  89  75   3  82 154 104
# 14  98  13  68 126  51 152 100  79 157 118  90  80 
apply(exp(newX), 1, max)[98]
apply(exp(FL), 1, max)[98]
apply(exp(U), 1, max)[98]

apply(exp(newX), 1, max)[100]
apply(exp(FL), 1, max)[100]
apply(exp(U), 1, max)[100]


570 + order(apply(exp(newX), 1, max)[571:697])
# 674 684 599 577 637 613 669 697 573
apply(exp(newX), 1, max)[669]
apply(exp(FL), 1, max)[669]
apply(exp(U), 1, max)[669]

apply(exp(newX), 1, max)[684]
apply(exp(FL), 1, max)[684]
apply(exp(U), 1, max)[684]

exp.select.X = exp(newX)[c(118, 157, 669, 684), ]
exp.select.FL = exp(FL)[c(118, 157, 669, 684), ]
exp.select.U = exp(U)[c(118, 157, 669, 684), ]

scaleX = (exp.select.X - min(exp.select.X))/(max(exp.select.X) - min(exp.select.X))*6
scaleFL = (exp.select.FL - min(exp.select.FL))/(max(exp.select.FL) - min(exp.select.FL))*6
scaleU = (exp.select.U- min(exp.select.U))/(max(exp.select.U) - min(exp.select.U))*6

scaleX[scaleX<0.1] = 0
scaleFL[scaleFL<0.1] = 0
scaleU[scaleU<0.1] = 0

writeMat("data7.mat", X = scaleX[,], U=scaleU[,], FL = scaleFL[,])

# min(exp(U[c(19,85,574,578),]))
# U[exp(U) < 1] = 0



