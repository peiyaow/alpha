library(xtable)

setwd("~/Documents/GitHub/alpha/")
setwd("./bonferroni/")
result_overall <- read.csv("result_overall.csv", header=FALSE)
result = rbind(apply(result_overall[,-16], 2, mean), apply(result_overall[,-16], 2, sd)/sqrt(99))
summary1 = matrix(cbind(result[,1:8], NA, result[,9:15]), ncol = 4)
colnames(summary1) = c("Global", "WLS", "Groupwise", "proposed_bonf")
rownames(summary1) = c("OLS", " ", "ridge", " ", "EN", " ", "lasso", " ")

setwd("~/Documents/GitHub/alpha/")
setwd("./scheffe_forward/")
result_overall <- read.csv("result_overall.csv", header=FALSE)
result = rbind(apply(result_overall[,-16], 2, mean), apply(result_overall[,-16], 2, sd)/sqrt(99))
summary2 = matrix(cbind(result[,1:8], NA, result[,9:15]), ncol = 4)
colnames(summary2) = c("Global", "WLS", "Groupwise", "proposed_schef.f")
rownames(summary2) = c("OLS", " ", "ridge", " ", "EN", " ", "lasso", " ")

setwd("~/Documents/GitHub/alpha/")
setwd("./scheffe_backward/")
result_overall <- read.csv("result_overall.csv", header=FALSE)
result = rbind(apply(result_overall[,-16], 2, mean), apply(result_overall[,-16], 2, sd)/sqrt(99))
summary3 = matrix(cbind(result[,1:8], NA, result[,9:15]), ncol = 4)
colnames(summary3) = c("Global", "WLS", "Groupwise", "proposed_schef.b")
rownames(summary3) = c("OLS", " ", "ridge", " ", "EN", " ", "lasso", " ")

summary = cbind(summary1, summary2[,"proposed_schef.f"], summary3[,"proposed_schef.b"])
colnames(summary) = c("Global", "WLS", "Groupwise", "proposed_bonf", "proposed_schef.f", "proposed_schef.b")
xtable(summary, digits = 3)

# classwise
setwd("~/Documents/GitHub/alpha/")
setwd("./bonferroni/")
result_class <- read.csv("result_class.csv", header=FALSE)
result = rbind(apply(result_class[,-76], 2, mean), apply(result_class[,-76], 2, sd)/sqrt(99))
result.NC = result[,seq(1, 75, by = 5)]
result.SMC = result[,seq(2, 75, by = 5)]
result.eMCI = result[,seq(3, 75, by = 5)]
result.lMCI = result[,seq(4, 75, by = 5)]
result.AD = result[,seq(5, 75, by = 5)]

summary1.NC = matrix(cbind(result.NC[,1:8], NA, result.NC[,9:15]), ncol = 4)
summary1.SMC = matrix(cbind(result.SMC[,1:8], NA, result.SMC[,9:15]), ncol = 4)
summary1.eMCI = matrix(cbind(result.eMCI[,1:8], NA, result.eMCI[,9:15]), ncol = 4)
summary1.lMCI = matrix(cbind(result.lMCI[,1:8], NA, result.lMCI[,9:15]), ncol = 4)
summary1.AD = matrix(cbind(result.AD[,1:8], NA, result.AD[,9:15]), ncol = 4)

summary1.OLS = rbind(summary1.NC[1:2,], summary1.SMC[1:2,], summary1.eMCI[1:2,], summary1.lMCI[1:2,], summary1.AD[1:2,])
summary1.ridge = rbind(summary1.NC[3:4,], summary1.SMC[3:4,], summary1.eMCI[3:4,], summary1.lMCI[3:4,], summary1.AD[3:4,])
summary1.EN = rbind(summary1.NC[5:6,], summary1.SMC[5:6,], summary1.eMCI[5:6,], summary1.lMCI[5:6,], summary1.AD[5:6,])
summary1.lasso = rbind(summary1.NC[7:8,], summary1.SMC[7:8,], summary1.eMCI[7:8,], summary1.lMCI[7:8,], summary1.AD[7:8,])

summary1 = rbind(summary1.OLS, summary1.ridge, summary1.EN, summary1.lasso)


colnames(summary1) = c("Global", "WLS", "Groupwise", "proposed_bonf")
rownames(summary1) = rep(c("NC", " ", "SMC", " ", "eMCI", " ",  "lMCI", " ", "AD", " "), 4)
xtable(summary1, digits = 3)

