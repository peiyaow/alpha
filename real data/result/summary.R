result_overall_ridge <- read.csv("~/Documents/GitHub/alpha/real data/result/result_overall_ridge.csv", header=FALSE)
round(c(apply(result_overall_ridge, 2, mean), apply(result_overall_ridge, 2, sd)/10), digits = 3)[c(1, 5, 2, 6, 3, 7, 4, 8)]

result_overall_EN <- read.csv("~/Documents/GitHub/alpha/real data/result/result_overall_EN.csv", header=FALSE)
round(c(apply(result_overall_EN, 2, mean), apply(result_overall_EN, 2, sd)/10), digits = 3)[c(1, 5, 2, 6, 3, 7, 4, 8)]

result_overall_lasso <- read.csv("~/Documents/GitHub/alpha/real data/result/result_overall_lasso.csv", header=FALSE)
round(c(apply(result_overall_lasso, 2, mean), apply(result_overall_lasso, 2, sd)/10), digits = 3)[c(1, 5, 2, 6, 3, 7, 4, 8)]

result_class_ridge <- read.csv("~/Documents/GitHub/alpha/real data/result/result_class_ridge.csv", header=FALSE)
round(c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[1,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[2,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[2,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[3,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[3,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[4,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[4,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[5,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[4,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
# c(matrix(apply(result_class_ridge, 2, mean), ncol = 4)[5,], matrix(apply(result_class_ridge, 2, sd)/10, ncol = 4)[5,])[c(1, 5, 2, 6, 3, 7, 4, 8)]

result_class_EN <- read.csv("~/Documents/GitHub/alpha/real data/result/result_class_EN.csv", header=FALSE)
round(c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[1,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[2,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[3,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[4,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[5,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[2,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[2,])[c(1, 5, 2, 6, 3, 7, 4, 8)]
c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[3,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[3,])[c(1, 5, 2, 6, 3, 7, 4, 8)]
c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[4,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[4,])[c(1, 5, 2, 6, 3, 7, 4, 8)]
c(matrix(apply(result_class_EN, 2, mean), ncol = 4)[5,], matrix(apply(result_class_EN, 2, sd)/10, ncol = 4)[5,])[c(1, 5, 2, 6, 3, 7, 4, 8)]

result_class_lasso<- read.csv("~/Documents/GitHub/alpha/real data/result/result_class_lasso.csv", header=FALSE)
round(c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[1,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[2,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[3,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[4,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)
round(c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[5,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[1,])[c(1, 5, 2, 6, 3, 7, 4, 8)], digits = 3)

#c(matrix(apply(result_class_lasso, 2, mean), ncol = 4)[5,], matrix(apply(result_class_lasso, 2, sd)/10, ncol = 4)[5,])[c(1, 5, 2, 6, 3, 7, 4, 8)]
