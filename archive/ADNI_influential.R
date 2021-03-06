load("/Users/MonicaW/Documents/GitHub/alpha/data/ADNI1.RData")
X1 = X
Y1 = Y
label1 = label
X1 = X1[label1!=4,]
Y1 = Y1[label1!=4]
label1 = label1[label1!=4]

load("/Users/MonicaW/Documents/GitHub/alpha/data/ADNI2.RData")
X2 = X
Y2 = Y
label2 = label
X2 = X2[label2!=4,]
Y2 = Y2[label2!=4]
label2 = label2[label2!=4]

n1 = nrow(X1)
n2 = nrow(X2)

X = rbind(X1, X2)
Y = c(Y1, Y2)
label = as.factor(c(rep(1, n1), rep(2, n2)))

X.mean = apply(X, 2, mean)
X.sd =  apply(X, 2, sd)
X = sweep(sweep(X, 2, X.mean), 2, X.sd, "/")
n = n1+n2
ix = 1:n

ix1 = ix[label==1]
X1 = X[label==1,]
res1 = plotPC(X1)
ix1.remove = ix1[res1$F2[,1]< -2.2]

X2 = X[label==2,]
ix2 = ix[label==2]
res2 = plotPC(X2)
ix2.remove = ix2[res2$F2[,1]> 2]

ix.remove = c(ix1.remove, ix2.remove)
X = X[-ix.remove, ]
Y = Y[-ix.remove]
label = label[-ix.remove]

summary(Y1)
summary(Y2)
summary(label)
plot(Y2)
plot(log(Y+1)~label)
plot(Y~label)



n = nrow(X)
ix = 1:n

ix1 = ix[label==1]
X1 = X[label==1,]
res1 = plotPC(X1)

X2 = X[label==2,]
ix2 = ix[label==2]
res2 = plotPC(X2)
ix2.remove = ix2[res2$F2[,2]> 3]

X = X[-ix2.remove, ]
Y = Y[-ix2.remove]
label = label[-ix2.remove]


n = nrow(X)
ix = 1:n

ix1 = ix[label==1]
X1 = X[label==1,]
res1 = plotPC(X1)

X2 = X[label==2,]
ix2 = ix[label==2]
res2 = plotPC(X2)
plot(res2$F_[,2],res2$F_[,3])
ix2.remove = ix2[abs(res2$F_[,3]) > 3]

X = X[-ix2.remove, ]
Y = Y[-ix2.remove]
label = label[-ix2.remove]


n = nrow(X)
ix = 1:n

ix1 = ix[label==1]
X1 = X[label==1,]
res1 = plotPC(X1)

X2 = X[label==2,]
ix2 = ix[label==2]
res2 = plotPC(X2)
plot(res2$F_[,2],res2$F_[,3])
ix2.remove = c(ix2[abs(res2$F_[,3]) > 3], ix2[res2$F_[,2] > 3])

save(X, Y, label, file = "ADNI_clean.RData")






