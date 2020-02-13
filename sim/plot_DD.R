#------------------------------------- groupwise model ---------------------------------------------
# ridge
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/ridge/DD.RData")

n = 50
x = seq(0, 5, by = 0.1)
# x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 2, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lwd = 2, lty=c(5, 4, 3, 2), cex=1.2)

# lasso
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/lasso/DD.RData")

n = 50
x = seq(0, 5, by = 0.1)
# x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 2, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lwd = 2, lty=c(5, 4, 3, 2), cex=1.2)

#------------------------------------- proposed model ---------------------------------------------
# ridge entire
load("~/Documents/GitHub/alpha/sim/results/proposed_model/ridge/DD.RData")

n = 50
x = seq(0, 5, by = 0.1)
x = x[x <= 3]
KK = length(x)-1
# DD = log(DD)

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean), max(AA_mean)), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 2, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lwd = 2, lty=c(5, 4, 3, 2), cex=1.2)

# lasso
load("~/Documents/GitHub/alpha/sim/results/proposed_model/lasso/DD.RData")

n = 50
x = seq(0, 5, by = 0.1)
x = x[x <= 2]
KK = length(x)-1
# DD = log(DD)

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean), max(AA_mean)), 
     lwd = 2, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 2, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 2, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 2, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lwd = 2, lty=c(5, 4, 3, 2), cex=1.2)

