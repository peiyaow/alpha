gamma_g = function(gamma){return (1/gamma/4- 1/2)}
x = seq(50, 500)/200
y = gamma_g(x)
plot(x, y, type = "l")
x1 = seq(501, 1000)/200
y1 = -1/x1
plot(x1, y1, type = "l")

xx = c(x, x1)
yy = c(y, y1)
plot(xx, yy, type = "l")

AA[1:10+7*10, ]

AA_mean = do.call(rbind, lapply(0:4, function(k) apply(AA[1:5+k*5, ], 2, mean)))

AA_sd = do.call(rbind, lapply(0:4, function(k) apply(AA[1:5+k*5, ], 2, sd)))/sqrt(5)


plot(seq(0, 2, by = 0.5), AA_mean[,2], type = "l", ylim = c(9, 20), lwd = 1, lty = 4, col = "blue", xlab = "beta", ylab = "error")
polygon(c(seq(0, 2, by = 0.5), rev(seq(0, 2, by = 0.5))), c(AA_mean[,2]+AA_sd[,2], rev(AA_mean[,2]-AA_sd[,2])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
# lines(AA_mean[,2], lwd = 1, lty = 2)
#lines(AA_mean[,2]+AA_sd[,2], col = "red", lty = 2)
#lines(AA_mean[,2]-AA_sd[,2], col = "red", lty = 2)

i = 1
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(seq(0, 2, by = 0.5), rev(seq(0, 2, by = 0.5))), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(seq(0, 2, by = 0.5), col = "red", AA_mean[,i], lwd = 1, lty = 3)
#lines(AA_mean[,i]+AA_sd[,i], col = "red", lty = 2)
#lines(AA_mean[,i]-AA_sd[,i], col = "red", lty = 2)

legend("topleft", inset=.1, legend=c("Proposed", "Groupwise"),
       col=c("blue", "red"), lty=4:3, cex=1)






BB_mean = do.call(rbind, lapply(0:8, function(k) apply(BB[1:5+k*5, ], 2, mean)))
BB_sd = do.call(rbind, lapply(0:8, function(k) apply(BB[1:5+k*5, ], 2, sd)))/sqrt(5)
x = seq(0, 4, by = 0.5)[1:6]
AA_mean = BB_mean[1:6,]
AA_sd = BB_sd[1:6,]

plot(x, AA_mean[,2], type = "l", ylim = c(1, 4), lwd = 1, lty = 4, col = "blue", xlab = "beta blow factor", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,2]+AA_sd[,2], rev(AA_mean[,2]-AA_sd[,2])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
# lines(AA_mean[,2], lwd = 1, lty = 2)
#lines(AA_mean[,2]+AA_sd[,2], col = "red", lty = 2)
#lines(AA_mean[,2]-AA_sd[,2], col = "red", lty = 2)

i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 3)
#lines(AA_mean[,i]+AA_sd[,i], col = "red", lty = 2)
#lines(AA_mean[,i]-AA_sd[,i], col = "red", lty = 2)

i = 1
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "green", AA_mean[,i], lwd = 1, lty = 5)

legend("top", inset=.08, legend=c("Global", "Proposed", "Groupwise"),
       col=c("green", "blue", "red"), lty=c(5, 4, 3), cex=1)





# sigma 
CC_mean = do.call(rbind, lapply(0:9, function(k) apply(CC[1:5+k*5, ], 2, mean)))
CC_sd = do.call(rbind, lapply(0:9, function(k) apply(CC[1:5+k*5, ], 2, sd)))/sqrt(5)
x = seq(log(0.5), log(4), length.out = 10)

exp(seq(log(0.1), log(10), length.out = 10))[4:8]
AA_mean = CC_mean
AA_sd = CC_sd

plot(x, AA_mean[,2], type = "l", ylim = c(1, 4), lwd = 1, lty = 4, col = "blue", xlab = "log(Lambda blow factor)", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,2]+AA_sd[,2], rev(AA_mean[,2]-AA_sd[,2])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
# lines(AA_mean[,2], lwd = 1, lty = 2)
#lines(AA_mean[,2]+AA_sd[,2], col = "red", lty = 2)
#lines(AA_mean[,2]-AA_sd[,2], col = "red", lty = 2)

i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 3)
#lines(AA_mean[,i]+AA_sd[,i], col = "red", lty = 2)
#lines(AA_mean[,i]-AA_sd[,i], col = "red", lty = 2)

i = 1
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "green", AA_mean[,i], lwd = 1, lty = 5)

legend("top", inset=.08, legend=c("Global", "Proposed", "Groupwise"),
       col=c("green", "blue", "red"), lty=c(5, 4, 3), cex=1)





DD_mean = do.call(rbind, lapply(0:12, function(k) apply(DD[1:5+k*5, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:12, function(k) apply(DD[1:5+k*5, ], 2, sd)))/sqrt(5)
x = seq(0, 6, by = 0.5)

AA_mean = DD_mean
AA_sd = DD_sd

plot(x, AA_mean[,2], type = "l", 
     ylim = c(25, 70), 
     lwd = 1, lty = 4, col = "blue", xlab = "log(Lambda blow factor)", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,2]+AA_sd[,2], rev(AA_mean[,2]-AA_sd[,2])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)
# lines(AA_mean[,2], lwd = 1, lty = 2)
#lines(AA_mean[,2]+AA_sd[,2], col = "red", lty = 2)
#lines(AA_mean[,2]-AA_sd[,2], col = "red", lty = 2)

i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 3)
#lines(AA_mean[,i]+AA_sd[,i], col = "red", lty = 2)
#lines(AA_mean[,i]-AA_sd[,i], col = "red", lty = 2)

i = 1
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "green", AA_mean[,i], lwd = 1, lty = 5)

legend("top", inset=.08, legend=c("Global", "Groupwise", "Proposed"),
       col=c("green", "blue", "red"), lty=c(5, 4, 3), cex=1)


       