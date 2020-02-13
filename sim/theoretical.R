n = 100
p = seq(10, 100000, by = 10)
mp = 1
wn = 1/sqrt(p)+sqrt(log(p)/n)
plot((1+wn*mp)*(sqrt(log(n)*log(p))/n+ 1/n^(1/4)/sqrt(p)) + mp/sqrt(n)/sqrt(p))
abline(h = 1/n)

# group vs proposed
n = 100
p = seq(10, 1000, by = 10)
delta = 1
mp = 1
wn = 1/sqrt(p)+sqrt(log(p)/n)
plot(p, sqrt(log(n)*log(p))/n + 1/n^(1/4)/sqrt(p) + mp*wn/sqrt(n), type = "l",
     ylab = "loss", col = "purple", lwd = 2, lty = 2)
abline(h=delta/sqrt(n)+1/n, col = "red", lwd = 2, lty = 4)
legend("topright", inset=.1, legend=c("Group", "Factor"),
       col=c("red", "purple"), lty=c(4, 2), lwd = 2, cex=1.2)


n = seq(100, 1000, by = 10)
p = 200
delta = 1
mp = 1
wn = 1/sqrt(p)+sqrt(log(p)/n)
plot(n, delta/sqrt(n)+1/n, ylim = c(0.02,.12),  type = "l", 
     col = "red", ylab = "loss", lwd = 2, lty = 4)
lines(n, sqrt(log(n)*log(p))/n + 1/n^(1/4)/sqrt(p) + mp*wn/sqrt(n), type = "l",
      ylab = "loss", col = "purple", lwd = 2, lty = 2)
legend("topright", inset=.1, legend=c("Group", "Factor"),
       col=c("red", "purple"), lwd = 2, lty=c(4, 2), cex=1.2)


n = 100
p = 200
delta = seq(1/sqrt(n), 2, by = 0.01)
mp = 1
wn = 1/sqrt(p)+sqrt(log(p)/n)
plot(delta, delta/sqrt(n)+1/n, type = "l", 
     col = "red", ylab = "loss", lwd = 2, lty = 4)
abline(h = sqrt(log(n)*log(p))/n + 1/n^(1/4)/sqrt(p) + mp*wn/sqrt(n), col = "purple", lwd = 2, lty = 2)
legend("topleft", inset=.1, legend=c("Group", "Factor"),
       col=c("red", "purple"), lwd = 2, lty=c(4, 2), cex=1.2)


