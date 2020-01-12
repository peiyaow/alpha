x1 = rnorm(100)
x2 = rnorm(100)
x3 = rnorm(100)

F.all = X2F(cbind(x1,x2,x3))$F_
FF = X2U2(cbind(x1, x2), K = 1)$F_
U = X2U2(cbind(x1, x2), K = 1)$U

y = FF%*%c(1, 2) + U%*%c(1,1)
cor(y, F.all[,2:100])
cor(y, cbind(x1,x2))

lm(y~FF+U)
lm(y~x1+x2)


f1 = FF[,2]
f2 = FF[,3]
cov(y, x1)/sqrt(var(y))/sqrt(var(x1))
cov(y, x2)/sqrt(var(y))/sqrt(var(x2))
lm(y~x1+x2)
cov(y, f1)/sqrt(var(y))/sqrt(var(f1))
cov(y, f2)/sqrt(var(y))/sqrt(var(f2))
