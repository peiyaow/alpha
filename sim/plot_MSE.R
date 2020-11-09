file_path = "~/Documents/alpha_longleaf/sim/"
penalty = 'ridge'
result_type = 'subgaussian'

s.vec = seq(0, 5, by = 0.1)
MSE.train.list = list()
MSE.test.list = list()
count = 1
for (s in s.vec){
  file.name = paste0("result_", penalty, "_", result_type, "_s=", format(s, nsmall = 1), ".RData")
  load(paste0(file_path, penalty, "_", result_type, "/", file.name))
  MSE.train.list[[count]] = t(sapply(MSE.train, function(MSE) apply(MSE, 1, mean)))
  MSE.test.list[[count]] = t(sapply(MSE.test, function(MSE) apply(MSE, 1, mean)))
  count = count + 1
}

MSE.train = do.call(rbind, MSE.train.list)
MSE.test = do.call(rbind, MSE.test.list)

# DIFF.beta
# DIFF.gamma
n_sim = 50
s.vec = s.vec[s.vec < s]
KK = length(s.vec)-1

MSE = MSE.test

MSE_mean = do.call(rbind, lapply(0:KK, function(k) apply(MSE[1:n_sim+k*n_sim, ], 2, mean)))
MSE_sd = do.call(rbind, lapply(0:KK, function(k) apply(MSE[1:n_sim+k*n_sim, ], 2, sd)))/sqrt(n_sim)

# global
i = 1
plot(s.vec, MSE_mean[,i], type = "l", 
     ylim = c(min(MSE_mean)-0.3, max(MSE_mean)+0.3), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(s.vec, rev(x)), c(MSE_mean[,i]+MSE_sd[,i], rev(MSE_mean[,i]-MSE_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(MSE_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(s.vec, rev(s.vec)), c(MSE_mean[,i]+MSE_sd[,i], rev(MSE_mean[,i]-MSE_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(s.vec, col = "red", MSE_mean[,i], lwd = 2, lty = 4)

# factor-0
i = 3
#plot(MSE_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(s.vec, rev(s.vec)), c(MSE_mean[,i]+MSE_sd[,i], rev(MSE_mean[,i]-MSE_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(s.vec, col = "blue", MSE_mean[,i], lwd = 2, lty = 3)

# factor
i = 4
#plot(MSE_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(s.vec, rev(s.vec)), c(MSE_mean[,i]+MSE_sd[,i], rev(MSE_mean[,i]-MSE_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(s.vec, col = "purple", MSE_mean[,i], lwd = 2, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lwd = 2, lty=c(5, 4, 3, 2), cex=1.2)

  
  
  
  
  