library(readr)
X = as.matrix(read_table("~/Documents/GitHub/alpha/data/X2.txt", col_names = F))
Y = read_table("~/Documents/GitHub/alpha/data/Y2.txt", col_names = F)$X1
label = as.ordered(read_table("~/Documents/GitHub/alpha/data/label2.txt", col_names = F)$X1)
summary(label)

X = X[!is.na(Y),]
label = label[!is.na(Y)]
Y = Y[!is.na(Y)]

save(X, Y, label, file = "ADNI2.RData")
