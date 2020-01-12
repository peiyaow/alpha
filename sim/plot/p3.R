# ridge
gamma = 0:15/10
s = 0

f1 = 1
f2 = 1/4+gamma/2
f3 = 1-gamma/2
f4 = 3/2-gamma
f5 = 2-gamma*(s+1)

nn = 1000
plot(gamma, apply(rbind(f1, f2, f3, f4, f5), 2, min), type = "l")

nn = 1000
plot(gamma, nn^(-apply(rbind(f1, f2, f3, f4, f5), 2, min)), type = "l", ylab = "Risk", main = "s = 0")

plot(gamma, apply(rbind(f1, f2, f3, f4, f5), 2, which.min))

# lasso
gamma = seq(0, 2, by = 0.1)
rs = 0
f1 = 1
f2 = 1/4+gamma/2
f3 = 3/2-gamma*(1/2+rs)
f4 = 1-gamma/2

nn = 1000
plot(gamma, apply(rbind(f1, f2, f3, f4), 2, min), type = "l")
nn = 1000
plot(gamma, nn^(-apply(rbind(f1, f2, f3, f4), 2, min)), type = "l", ylab = "Risk", main = "rs = 0")
plot(gamma, apply(rbind(f1, f2, f3, f4), 2, which.min))





# 


r1 = seq(3/4, 2, by = 0.01)
rs1 = 1/2/r1

plot(r1, rs1, xlab = "r", ylab = "rs", xlim = c(0, 2), ylim = c(0, 50),type = "l", lty = 2)

r2 = seq(0, 3/4, by = 0.01)
rs22 = 5/4/r2-1
lines(r2, rs22, lty = 6)

rr = c(r2, r1)
rs2 = 3/2/rr -1/2
lines(rr, rs2, lty = 5)
abline(v = 0.75, lty = 6)
abline(v = 2, lty = 3)
#text(x = 0.3, y = 1, labels = "1/4+r/2")
text(x = 0.35, y = 1, labels = "A")
#text(x = 1.35, y = 0.18, labels = "1-r/2")
text(x = 1.35, y = 0.16, labels = "B")
#text(x = 1, y = 0.9, labels = "3/2-r(1/2+rs)")
text(x = 1, y = 0.75, labels = "C")


rs = seq(0, 1, by = 0.05)
r1 = 1/2/rs
r2 = 3/2/(1/2+rs)

