n = 50
DD_mean = do.call(rbind, lapply(0:12, function(k) apply(DD[1:n+k*n, ], 2, mean)))
DD_sd = do.call(rbind, lapply(0:12, function(k) apply(DD[1:n+k*n, ], 2, sd)))/sqrt(n)
x = seq(0, 6, by = 0.5)

AA_mean = DD_mean
AA_sd = DD_sd

plot(x, AA_mean[,2], type = "l", 
     ylim = c(28, 70), 
     lwd = 1, lty = 4, col = "blue", xlab = "h", ylab = "error")
polygon(c(x, rev(x)), c(AA_mean[,2]+AA_sd[,2], rev(AA_mean[,2]-AA_sd[,2])), col = rgb(0, 0, 255, 20, maxColorValue=255), border = FALSE)

i = 3
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(255, 0, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "red", AA_mean[,i], lwd = 1, lty = 3)

i = 1
#plot(AA_mean[,i], type = "l", ylim = c(9, 20))
polygon(c(x, rev(x)), c(AA_mean[,i]+AA_sd[,i], rev(AA_mean[,i]-AA_sd[,i])), col = rgb(0, 255, 0, 20, maxColorValue=255), border = FALSE)
lines(x, col = "green", AA_mean[,i], lwd = 1, lty = 5)

legend("topleft", inset=.1, legend=c("Global", "Groupwise", "Proposed"),
       col=c("green", "blue", "red"), lty=c(5, 4, 3), cex=1)


