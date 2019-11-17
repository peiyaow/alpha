source('~/Documents/Research/coding/R/alpha/functions.R')
load("~/Documents/Research/coding/R/realdata/ADNI1.RData")
load("~/Documents/Research/coding/R/realdata/myseeds.RData")
###

X = X[-c(112, 151, 167, 171, 770),]
label = label[-c(112, 151, 167, 171, 770)]
Y = Y[-c(112, 151, 167, 171, 770)]

X = X[-c(770),]
label = label[-c(770)]
Y = Y[-c(770)]
X = X[-c(112, 157, 167),]
label = label[-c(112, 157, 167)]
Y = Y[-c(112, 157, 167)]



###
# _0
X = X[-c(167,770),]
label = label[-c(167,770)]
Y = Y[-c(167,770)]

X = X[-151,]
label = label[-151]
Y = Y[-151]

X = X[-c(112,169),]
label = label[-c(112,169)]
Y = Y[-c(112,169)]

X = X[-c(70,126),]
label = label[-c(70,126)]
Y = Y[-c(70,126)]

# X = X[-c(112,151,170),]
# label = label[-c(112,151,170)]
# Y = Y[-c(112,151,170)]
ix = 1:length(Y)

# bc = BoxCoxTrans(Y)
# Y = predict(bc, Y)
# bcb = function(Y, bc) {return((Y*bc$lambda+1)^{1/bc$lambda})}

boxplot(Y)

par(mfrow = c(3,1))

ix.NC = ix[label==0]
X1 = X[label==0,]
Y1 = Y[label==0]
boxplot(Y1)
X1.mean = apply(X1, 2, mean)
X1 = sweep(X1, 2, X1.mean)

F1.2 = X2U2(X1, plot = T)$F2
plot(F1.2)
ix.NC[F1.2[,1]< (-4)]
summary(lm(Y1~F1))
plot(F1, Y1)

ix.MCI = ix[label==3]
X2 = X[label==3,]
Y2 = Y[label==3]
boxplot(Y2)
sort(Y2)

X2.mean = apply(X2, 2, mean)
X2 = sweep(X2, 2, X2.mean)
F2 = X2U2(X2, plot = T)$F_

F2.1 = X2U(X2, plot = T)$F1
F2.2 = X2U(X2, plot = T)$F2
boxplot(F2.1)
boxplot(F2.2)
plot(F2.1, F2.2)
plot(F2, Y2)
summary(lm(Y2~F2))
boxplot(Y2)


ix.AD = ix[label==4]
X3 = X[label==4,]
Y3 = Y[label==4]
boxplot(Y3)
plot(Y3)


X3.mean = apply(X3, 2, mean)
X3 = sweep(X3, 2, X3.mean)
F3 = X2U(X3, plot = T)$F_

F3.1 = X2U(X3, plot = T)$F1
F3.2 = X2U(X3, plot = T)$F2
plot(F3.1, F3.2)
ix.AD[abs(F3.2)>3]
ix.AD[abs(F3.2)>2&F3.1 > 2]


ix.AD[abs(F3.2)>2.5]
F3.2[abs(F3.2)>2&F3.1 > 1]
F3.1[abs(F3.2)>2&F3.1 > 1]

boxplot(F3.2)
boxplot(F3.1)
sort(Y3)
plot(F3, Y3)
summary(lm(Y3~F3))

boxplot(Y3)
