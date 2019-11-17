n = nrow(X)
ix = 1:n

ix1 = ix[label==0]
X1 = X[label==0,]
res1 = plotPC(X1)
par(mfrow = c(1,1))
plot(res1$F2[,1], res1$F2[,2])
ix1.remove = ix1[res1$F2[,2]< -3]
sort(res1$F2[,2])

ix1 = ix[label==1]
X1 = X[label==1,]
res1 = plotPC(X1)
par(mfrow = c(1,1))
plot(res1$F2[,1], res1$F2[,2])
ix2.remove = ix1[res1$F2[,2]< -2]

ix1 = ix[label==2]
X1 = X[label==2,]
res1 = plotPC(X1)
par(mfrow = c(1,1))
plot(res1$F2[,1], res1$F2[,2])
ix3.remove = ix1[res1$F2[,2]> 2]

ix1 = ix[label==3]
X1 = X[label==3,]
res1 = plotPC(X1)
par(mfrow = c(1,1))
plot(res1$F2[,1], res1$F2[,2])
ix4.remove = ix1[res1$F2[,2]< -3]

ix1 = ix[label==4]
X1 = X[label==4,]
res1 = plotPC(X1)
par(mfrow = c(1,1))
plot(res1$F2[,1], res1$F2[,2])
ix5.remove = ix1[res1$F2[,1]< -1]

ix.remove = c(ix1.remove, ix2.remove, ix3.remove, ix4.remove, ix5.remove)

X = X[-ix.remove, ]
Y = Y[-ix.remove]
label = label[-ix.remove]

summary(label)



n = length(Y)
ix = 1:n

ix1 = ix[label==0]
X1 = X[label==0,]
res1 = plotPC(X1)
ix1.remove = ix1[res1$F2[,1]< -1 & res1$F2[,2]>1]

X2 = X[label==1,]
ix2 = ix[label==1]
res2 = plotPC(X2)
ix2.remove = ix2[res2$F2[,2]< -2]

X3 = X[label==2,]
ix3 = ix[label==2]
res3 = plotPC(X3)
ix3.remove = ix3[abs(res3$F2[,1]) > 2]

X4 = X[label==3,]
ix4 = ix[label==3]
res4 = plotPC(X4)
ix4.remove = ix4[res4$F2[,1] < -1.5]

X5 = X[label==4,]
ix5 = ix[label == 4]
res5 = plotPC(X5)
ix5.remove = ix5[res5$F2[,2] < -3]

X = X[-ix5.remove, ]
Y = Y[-ix5.remove]
label = label[-ix5.remove]
summary(label)
save(X, Y, label, file = "ADNI2_clean3.RData")
