library(readr)
library(rlist)
library(xtable)
setwd('~/Documents/GitHub/alpha/sim/333/123same/')
path = "~/Documents/GitHub/alpha/"
load(paste0(path, "/sim/parameter/para_grid.RData"))

beta_U.grid = c(0, 1, 3)
sigma.vec.grid = list(c(1, 1, 1), c(1, 2, 3))
p.grid = 80*c(1, 2, 3, 4, 5)

res.list = list()
for (beta in beta_U.grid){
  for (sigma.vec in sigma.vec.grid){
    for (p in p.grid){
      para = list(beta = beta, sigma.vec = sigma.vec, p = p)
      dir_name = paste0(para$beta, "_", para$sigma.vec[3], "_", sprintf('%03d', para$p))
      result_overall <- read_csv(paste0(dir_name, "/result_overall.csv"), col_names = FALSE)
      res = round(t(matrix(t(cbind(matrix(apply(result_overall, 2, function(x) mean(x, na.rm = T)), ncol = 4, byrow = T), matrix(apply(result_overall, 2, function(x) sd(x, na.rm = T)), ncol = 4, byrow = T))), ncol = 12)), 3)
      colnames(res) = c("OLS", "ridge", "EN", "lasso")
      rownames(res) = c("global", "", "WLS", "", "class", "", "Factor", "", "Factor.LD", "",  "Factor.HD", "")
      res.list = list.append(res.list, res)
    }
  }
}

res.list = list()
for (beta in beta_U.grid){
  for (sigma.vec in sigma.vec.grid){
    for (p in p.grid){
      para = list(beta = beta, sigma.vec = sigma.vec, p = p)
      dir_name = paste0(para$beta, "_", para$sigma.vec[3], "_", sprintf('%03d', para$p))
      result_class <- read_csv(paste0(dir_name, "/result_class.csv"), col_names = FALSE)
      # result_overall = result_class[,seq(1, 72, by = 3)]
      # result_overall = result_class[,seq(2, 72, by = 3)]
      result_overall = result_class[,seq(3, 72, by = 3)]
      res = round(t(matrix(t(cbind(matrix(apply(result_overall, 2, function(x) mean(x, na.rm = T)), ncol = 4, byrow = T), matrix(apply(result_overall, 2, function(x) sd(x, na.rm = T)), ncol = 4, byrow = T))), ncol = 12)), 3)
      colnames(res) = c("OLS", "ridge", "EN", "lasso")
      rownames(res) = c("global", "", "WLS", "", "class", "", "Factor", "", "Factor.LD", "",  "Factor.HD", "")
      res.list = list.append(res.list, res)
    }
  }
}
res1.list = res.list
res2.list = res.list
res3.list = res.list

xtable(do.call(rbind, res.list[1:5]))
xtable(do.call(rbind, res.list[6:10]))
xtable(do.call(rbind, res.list[11:15]))
xtable(do.call(rbind, res.list[16:20]))
xtable(do.call(rbind, res.list[21:25]))
xtable(do.call(rbind, res.list[26:30]))

xtable(cbind(do.call(rbind, res1.list[16:20]), do.call(rbind, res2.list[16:20]), do.call(rbind, res3.list[16:20])))

xtable(cbind(do.call(rbind, res1.list[26:30]), do.call(rbind, res2.list[26:30]), do.call(rbind, res3.list[26:30])))

# col.names = c('lm.global', 'ridge.global', 'EN.global', 'lasso.global', 'lm.WLS', 'ridge.WLS', 'EN.WLS', 'lasso.WLS', 'lm.class', 'ridge.class', 'EN.class', 'lasso.class', 'OLS.U', 'ridge.U', 'EN.U', 'lasso.U', 'WLS.LD', 'ridge.LD', 'EN.LD', 'lasso.LD', 'OLS.HD', 'ridge.HD', 'EN.HD', 'lasso.HD')

setwd('~/Documents/GitHub/alpha/sim/333/111same/')
path = "~/Documents/GitHub/alpha/"
load(paste0(path, "/sim/parameter/para_grid.RData"))

beta_U.grid = c(0, 1, 3)
sigma.vec.grid = list(c(1, 1, 1), c(1, 2, 3))
p.grid = 80*c(1, 2, 3, 4, 5)

res.list = list()
for (beta in beta_U.grid){
  for (sigma.vec in sigma.vec.grid){
    for (p in p.grid){
      para = list(beta = beta, sigma.vec = sigma.vec, p = p)
      dir_name = paste0(para$beta, "_", para$sigma.vec[3], "_", sprintf('%03d', para$p))
      result_overall <- read_csv(paste0(dir_name, "/result_overall.csv"), col_names = FALSE)
      res = round(t(matrix(t(cbind(matrix(apply(result_overall, 2, function(x) mean(x, na.rm = T)), ncol = 4, byrow = T), matrix(apply(result_overall, 2, function(x) sd(x, na.rm = T)), ncol = 4, byrow = T))), ncol = 12)), 3)
      colnames(res) = c("OLS", "ridge", "EN", "lasso")
      rownames(res) = c("global", "", "WLS", "", "class", "", "Factor", "", "Factor.LD", "",  "Factor.HD", "")
      res.list = list.append(res.list, res[,])
    }
  }
}

xtable(do.call(rbind, res.list[11:15]))
xtable(do.call(rbind, res.list[16:20]))
# -----------------------------------------------------------------------------


ix.list = list()
ix.list = list.append(ix.list, 1:50, 51:100, 101:150, 151:200)

res.list = lapply(1:4, function(ix) rbind(apply(result_overall[ix.list[[ix]], ], 2, mean), apply(result_overall[ix.list[[ix]], ], 2, sd)))
res = round(do.call(rbind, res.list), 3)
colnames(res) = col.names 

# rbind(res[,1:4], res[,5:8], cbind(0, res[,9:11]), res[,12:15])
xtable(res[,c('lasso.global', 'lasso.WLS', 'lasso.class', 'lasso.proposed')])



# 80
rbind(apply(result_overall[1:50, ], 2, mean), apply(result_overall[1:50, ], 2, sd))

# 160
rbind(apply(result_overall[51:100, ], 2, mean), apply(result_overall[51:100, ], 2, sd))

# 240
rbind(apply(result_overall[101:150, ], 2, mean), apply(result_overall[101:150, ], 2, sd))

# 320
rbind(apply(result_overall[151:200, ], 2, mean), apply(result_overall[151:200, ], 2, sd))

# beta = 3
result_overall[13:16, ]
