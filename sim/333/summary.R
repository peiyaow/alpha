library(readr)
library(rlist)
library(xtable)
setwd('~/Documents/GitHub/alpha/sim/333')
result_overall <- read_csv("5/result_overall.csv", 
                           col_names = FALSE)

#lm.global, ridge.global, EN.global, lasso.global # 1 2 3 4
#lm.WLS, ridge.WLS, EN.WLS, lasso.WLS, # 5 6 7 8 
#ridge.X.class, EN.X.class, lasso.X.class, # 9 10 11
#OLS.WLS, ridge.WLS, EN.WLS, lasso.WLS # 12 13 14 15

col.names = c('lm.global', 'ridge.global', 'EN.global', 'lasso.global', 'lm.WLS', 'ridge.WLS', 'EN.WLS', 'lasso.WLS', 'ridge.class', 'EN.class', 'lasso.class', 'OLS.proposed', 'ridge.proposed', 'EN.proposed', 'lasso.proposed')

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
