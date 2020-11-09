library(fields)
file_path = "~/Documents/alpha_longleaf/sim_v2/"
penalty = 'ridge'
len_lambda = 50
len_C = 21
s = 4
file.name = paste0("result_", penalty, "_C_s=", format(s, nsmall = 1), ".csv")
MSE = as.matrix(read.csv(paste0(file_path, penalty, "_C/", file.name), header=FALSE))
# dim(MSE)

# list of nsim; per element is a matrix len_C by len_lambda
MSE.mtx.list = lapply(1:nrow(MSE), function(ix) 
  matrix(apply(matrix(MSE[ix,-1], ncol = 3), 1, mean), ncol = len_lambda))

# MSE.mtx.list = lapply(1:nrow(MSE), function(ix) 
#  matrix(apply(t(matrix(MSE[ix,-1], nrow = 3)), 1, mean), ncol = len_lambda))

# plot the matrix 
# row: C; column: lambda
image.plot(x = seq(0, 1, length.out = 21),
           y = seq(1, 50),
           z = Reduce("+", MSE.mtx.list)/length(MSE.mtx.list), 
           col = hcl.colors(50, palette = "Inferno", rev = T),
           horizontal = F,
           xlab = 'C', 
           ylab = 'lambda',
           main = paste0(penalty, "_s=", format(s, nsmall = 1))
          )


# fix lambda tune C
# fix C tune lambda
# tune both lambda and C


# CV
file.name = paste0("result_", penalty, "cv_C_s=", format(s, nsmall = 1), ".csv")
cvMSE = as.matrix(read.csv(paste0(file_path, penalty, "_C/", file.name), header=FALSE))
ix.mtx = t(apply(cvMSE, 1, function(row) row[1:2]))

C_path = seq(0, 1, length.out = 21) # 21
step_C = C_path[2] - C_path[1]
for (i in 1:nrow(ix.mtx)){
  C_ix = ix.mtx[i,1]
  C = C_path[C_ix]
  lambda = ix.mtx[i,2]
  C = runif(n = 1, min = C-step_C/2, max = C+step_C/2)
  lambda = runif(n = 1, min = lambda-1/2, max = lambda+1/2)
  points(x = C, y = lambda, pch = 20, cex = 0.7, col = "#040404")
}

cvMSE.mtx = t(apply(cvMSE, 1, function(row) apply(matrix(row[-(1:2)], ncol = 3), 1, mean))
)
apply(cvMSE.mtx, 2, mean)

