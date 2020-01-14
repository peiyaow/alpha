#------------------------------------- groupwise model ---------------------------------------------
# ridge 111
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/ridge/DD.RData")

n = 50
x = seq(0, 2.5, by = 0.1)
x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

# ridge 112
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/112/ridge/DD.RData")

n = 50
x = seq(0, 2.5, by = 0.1)
x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)


# lasso 111
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/lasso/DD.RData")

n = 50
x = seq(0, 2.5, by = 0.1)
x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Groupwise", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

# lasso 112
load("~/Documents/GitHub/alpha/sim/results/groupwise_model/112/lasso/DD.RData")

n = 50
x = seq(0, 2.5, by = 0.1)
x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

#------------------------------------- proposed model ---------------------------------------------
# ridge entire
load("~/Documents/GitHub/alpha/sim/results/proposed_model/ridge/entire/DD.RData")

n = 50
x = seq(0, 3, by = 0.1)
#x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Groupwise", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

# ridge partial
load("~/Documents/GitHub/alpha/sim/results/proposed_model/ridge/partial/DD.RData")

n = 50
x = seq(0, 3, by = 0.1)
x = x[x <= 1.8]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

# lasso entire
load("~/Documents/GitHub/alpha/sim/results/proposed_model/lasso/entire/DD.RData")

n = 50
x = seq(0, 3, by = 0.1)
#x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Groupwise", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)

# lasso partial
load("~/Documents/GitHub/alpha/sim/results/proposed_model/lasso/partial/DD.RData")

n = 50
x = seq(0, 3, by = 0.1)
x = x[x <= 1.8]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd

# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# factor-0
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

# factor
i = 4
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "purple", AA_mean[,i], lwd = 1, lty = 2)

legend("topleft", inset=.1, legend=c("Global", "Group", "Factor-0", "Factor"),
       col=c("green", "red", "blue", "purple"), lty=c(5, 4, 3, 2), cex=1)



# proposed model
# lasso
load("/Users/MonicaW/Documents/GitHub/alpha/sim/results/proposed_model/lasso/both/DD.RData")
load("/Users/MonicaW/Documents/GitHub/alpha/sim/results/proposed_model/lasso/diff/DD.RData")
load("/Users/MonicaW/Documents/GitHub/alpha/sim/results/proposed_model/lasso/same/DD.RData")

# ridge
load("/Users/MonicaW/Documents/GitHub/alpha/sim/results/proposed_model/ridge/both/DD.RData")
load("/Users/MonicaW/Documents/GitHub/alpha/sim/results/proposed_model/ridge/diff/DD.RData")

n = 80
x = seq(0, 40, length.out = 25)/sqrt(p)
x = x[x <= 1]
KK = length(x)-1

DD_mean = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:KK, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)

AA_mean = DD_mean
AA_sd = DD_sd


# global
i = 1
plot(x, AA_mean[,i], type = "l", 
     ylim = c(min(AA_mean)-0.3, max(AA_mean)+0.3), 
     lwd = 1, lty = 5, col = "green", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)

# groupwise
i = 2
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 4)

# proposed
i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
lines(x, col = "blue", AA_mean[,i], lwd = 1, lty = 3)

legend("topleft", inset=.1, legend=c("Global", "Groupwise", "Proposed"),
       col=c("green", "red", "blue"), lty=c(5, 4, 3), cex=1)


