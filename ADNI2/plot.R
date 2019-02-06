library(ggplot2)

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
ggplot(data, aes(x=L1, y=L2, color = stage)) + geom_point() + ggtitle("Significant Factors")+ theme(plot.title = element_text(hjust = 0.5))
ggplot(data, aes(x=Lr1, y=Lr2, color = stage)) + geom_point() + ggtitle("Insignificant Factors")+ theme(plot.title = element_text(hjust = 0.5))

PY.list = list()
HY.list = list()
for (l in 1:4){
  PY.list[[l]] = X2U1(X.train.list[[l]])$P%*%Y.train.list[[l]]
  HY.list[[l]] = X2U1(X.train.list[[l]])$H%*%Y.train.list[[l]]
}


dat2 = data.frame(PY = do.call(c, PY.list), HY = do.call(c, HY.list), stage = c(rep("NC", n.train.vec[1]), rep("SMC", n.train.vec[2]), rep("eMCI", n.train.vec[3]), rep("lMCI", n.train.vec[4])))
index = 1:n.train
ggplot(dat2, aes(x= index, y=PY, color = stage)) + geom_point()+ggtitle("Plot of PY")+ theme(plot.title = element_text(hjust = 0.5))
ggplot(dat2, aes(x= index, y=HY, color = stage)) + geom_point()+ggtitle("Plot of HY")+ theme(plot.title = element_text(hjust = 0.5))


X2U1(as.matrix(X.test.list[[1]][14,] - X.train.mean[[1]]))$L2


