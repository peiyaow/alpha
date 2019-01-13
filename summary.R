# U separate select K
# U3 combine select K version1
# U4 combine select K version2

library(xtable)

getwd()
load("X2U3.RData")
mse.U3.result = mse.U.result
mse.U3.vec.result = mse.U.vec.result

load("X2U4.RData")
mse.U4.result = mse.U.result
mse.U4.vec.result = mse.U.vec.result

load("X2U4_scaleY.RData")
mse.U4_scaleY.result = mse.U.result
mse.U4_scaleY.vec.result = mse.U.vec.result

load("X2U.RData")
load("compare_X.RData")

# A
A.result = cbind(mse.U.result[,4:6], mse.U3.result[,4:6], mse.U4.result[,4:6])
colnames(A.result)[4:9] = c("lm.U3","lasso.U3","ridge.U3","lm.U4","lasso.U4","ridge.U4")
A.table = rbind(apply(A.result, 2, mean), apply(A.result, 2, sd))
row.names(A.table) = c("MSE", "SE")
xtable(t(A.table[, c("lm.U", "lm.U3", "lm.U4")]))
xtable(t(A.table[, c("lasso.U", "lasso.U3", "lasso.U4")]))
xtable(t(A.table[, c("ridge.U", "ridge.U3", "ridge.U4")]))

# B
B.result = cbind(mse.X.result[,], mse.U4.result[,4:6], mse.U4_scaleY.result[,4:6])
colnames(B.result)[15:17] = c("lm.U_scale","lasso.U_scale","ridge.U_scale")
B.result[,c("lm.U", "lm.U_scale")]
B.table = rbind(apply(B.result, 2, mean), apply(B.result, 2, sd))
row.names(B.table) = c("MSE", "SE")

xtable(t(B.table[,c(1:5,11)])) # groupwise
xtable(t(B.table[,6:10])) # global
xtable(t(B.table[,12:17])) # global


xtable(t(B.table[, c("lm.X.class", "lm.X", "lm.U", "lm.U_scale", "group_mean")]))
xtable(t(B.table[, c("lasso.X.class", "lasso.X", "lasso.U", "lasso.U_scale", "group_mean")]))
xtable(t(B.table[, c("ridge.X.class", "ridge.X", "ridge.U", "ridge.U_scale", "group_mean")]))

# conclusion: lm.X.class overfit

# C
C.result = cbind(mse.U.result[,1:3], mse.U4.result[,4:6])
#colnames(C.result)[7:9] = c("lm.U_scale","lasso.U_scale","ridge.U_scale")
C.table = rbind(apply(C.result, 2, mean), apply(C.result, 2, sd))
row.names(C.table) = c("MSE", "SE")
xtable(t(C.table))

C.AD.result = cbind(mse.U.vec.result[,,1][,1:3], mse.U4.vec.result[,,1][,4:6], mse.X.vec.result[,,1][,11])
colnames(C.AD.result) = c(colnames(C.result), "group_mean")
C.AD.table = rbind(apply(C.AD.result, 2, mean), apply(C.AD.result, 2, sd))
row.names(C.AD.table) = c("MSE", "SE")
xtable(t(C.AD.table))

C.MCI.result = cbind(mse.U.vec.result[,,2][,1:3], mse.U4.vec.result[,,2][,4:6], mse.X.vec.result[,,2][,11])
colnames(C.MCI.result) = c(colnames(C.result), "group_mean")
C.MCI.table = rbind(apply(C.MCI.result, 2, mean), apply(C.MCI.result, 2, sd))
row.names(C.MCI.table) = c("MSE", "SE")
xtable(t(C.MCI.table))

C.NC.result = cbind(mse.U.vec.result[,,3][,1:3], mse.U4.vec.result[,,3][,4:6], mse.X.vec.result[,,3][,11])
colnames(C.NC.result) = c(colnames(C.result), "group_mean")
C.NC.table = rbind(apply(C.NC.result, 2, mean), apply(C.NC.result, 2, sd))
row.names(C.NC.table) = c("MSE", "SE")
xtable(t(C.NC.table))


# D
# AD
D.AD.result = cbind(mse.X.vec.result[,,1], mse.U4.vec.result[,,1][,4:6], mse.U4_scaleY.vec.result[,,1][,4:6])
colnames(D.AD.result) = c(colnames(mse.X.result), colnames(mse.U4.result)[4:6], c("lm.U_scale","lasso.U_scale","ridge.U_scale"))
D.AD.table = rbind(apply(D.AD.result, 2, mean), apply(D.AD.result, 2, sd))
row.names(D.AD.table) = c("MSE", "SE")
xtable(t(D.AD.table[,c(1:5,11)]))
xtable(t(D.AD.table[,6:10]))
xtable(t(D.AD.table[,12:17]))
xtable(t(D.AD.table[, c("lm.X.class", "lm.X", "lm.U", "lm.U_scale", "group_mean")]))
xtable(t(D.AD.table[, c("lasso.X.class", "lasso.X", "lasso.U", "lasso.U_scale", "group_mean")]))
xtable(t(D.AD.table[, c("ridge.X.class", "ridge.X", "ridge.U", "ridge.U_scale", "group_mean")]))

# MCI
D.MCI.result = cbind(mse.X.vec.result[,,2], mse.U4.vec.result[,,2][,4:6], mse.U4_scaleY.vec.result[,,2][,4:6])
colnames(D.MCI.result) = c(colnames(mse.X.result), colnames(mse.U4.result)[4:6], c("lm.U_scale","lasso.U_scale","ridge.U_scale"))
D.MCI.table = rbind(apply(D.MCI.result, 2, mean), apply(D.MCI.result, 2, sd))
row.names(D.MCI.table) = c("MSE", "SE")
xtable(t(D.MCI.table[,c(1:5,11)]))
xtable(t(D.MCI.table[,6:10]))
xtable(t(D.MCI.table[,12:17]))
xtable(t(D.MCI.table[, c("lm.X.class", "lm.X", "lm.U", "lm.U_scale", "group_mean")]))
xtable(t(D.MCI.table[, c("lasso.X.class", "lasso.X", "lasso.U", "lasso.U_scale", "group_mean")]))
xtable(t(D.MCI.table[, c("ridge.X.class", "ridge.X", "ridge.U", "ridge.U_scale", "group_mean")]))

# NC 
D.NC.result = cbind(mse.X.vec.result[,,3], mse.U4.vec.result[,,3][,4:6], mse.U4_scaleY.vec.result[,,3][,4:6])
colnames(D.NC.result) = c(colnames(mse.X.result), colnames(mse.U4.result)[4:6], c("lm.U_scale","lasso.U_scale","ridge.U_scale"))
D.NC.table = rbind(apply(D.NC.result, 2, mean), apply(D.NC.result, 2, sd))
row.names(D.NC.table) = c("MSE", "SE")
xtable(t(D.NC.table[,c(1:5,11)]))
xtable(t(D.NC.table[,6:10]))
xtable(t(D.NC.table[,12:17]))

xtable(t(D.NC.table[, c("lm.X.class", "lm.X", "lm.U", "lm.U_scale", "group_mean")]))
xtable(t(D.NC.table[, c("lasso.X.class", "lasso.X", "lasso.U", "lasso.U_scale", "group_mean")]))
xtable(t(D.NC.table[, c("ridge.X.class", "ridge.X", "ridge.U", "ridge.U_scale", "group_mean")]))

# E








dim(B.table)

apply(haha[,c("lm.X", "lm.X.class", "lm.U")], 2, mean)
apply(haha[,c("lm.X", "lm.X.class", "lm.U")], 2, sd)


# groupwise result
apply(mse.U.vec.result[,,1], 2, mean)
apply(mse.U4.vec.result[,,3], 2, mean)
apply(mse.U4_scaleY.vec.result[,,3], 2, mean)

apply(mse.mean.vec.result, 2, mean)

apply(mse.U4_scaleY.result, 2, mean)
# compare U3 and U
apply(mse.U3.result, 2, mean)
apply(mse.U4.result, 2, mean)
apply(mse.U.result, 2, mean)

apply(mse.U3.result, 2, sd)
apply(mse.U4.result, 2, sd)
apply(mse.U.result, 2, sd)


U.result = cbind(mse.U3.result, mse.U.result)
U.result = U.result[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

U.result[U.result[,7] < U.result[,8], c(7,8)]

apply(mse.mean.vec.result, 2, mean)
apply(mse.U.vec.result[,,3], 2, mean)
apply(mse.U.result, 2, mean)

apply(mse.X.result, 2, mean)
apply(mse.X.vec.result[,,3], 2, mean)

rownames(mse.X.vec.list[[1]])

# 1/11
# class
load("/Users/MonicaW/Documents/Research/coding/R/alpha/compare_X_0.RData")
class1 = rbind(apply(mse.X.result[,1:4], 2, mean), apply(mse.X.result[,1:4], 2, sd))
load("/Users/MonicaW/Documents/Research/coding/R/alpha/svm_X_0.RData")
svm.class = rbind(apply(mse.X.result[,1:2], 2, mean), apply(mse.X.result[,1:2], 2, sd))
class.result = cbind(class1, svm.class)
xtable(t(class.result))

# global
load("/Users/MonicaW/Documents/Research/coding/R/alpha/compare_X_0.RData")
global1 = rbind(apply(mse.X.result[,6:9], 2, mean), apply(mse.X.result[,6:9], 2, sd))
load("/Users/MonicaW/Documents/Research/coding/R/alpha/svm_X_0.RData")
svm.global = rbind(apply(mse.X.result[,3:4], 2, mean), apply(mse.X.result[,3:4], 2, sd))
global.result = cbind(global1, svm.global)
xtable(t(global.result))
#U
load("/Users/MonicaW/Documents/Research/coding/R/alpha/scale_Y_X2U4_2.RData")
U.result = rbind(apply(mse.U.result, 2, mean), apply(mse.U.result, 2, sd))
xtable(t(U.result))
