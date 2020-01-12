library(MASS)
load("~/Documents/GitHub/alpha/para_333.RData")
n = 100
p = 240

sigma = 1
beta_U = c(rep(10, 3), rep(0, n_label), rep(0, p-6))
beta_F = c(rep(5, 3), rep(0, p-3))

L.list = list()
Sigma_U.list = list()
for (i in 1:(p/80)){
  L.list[[i]] = t(mvrnorm(80, mu_B, Sigma_B))
  Sigma_U.list[[i]] = Sigma_U
}
L = do.call(cbind, L.list)  
F_ = mvrnorm(n, mu = rep(0, K), diag(K))  
U = mvrnorm(n, mu = rep(0, p), bdiag(Sigma_U.list))

F_.test = mvrnorm(n, mu = rep(0, K), diag(K))
U.test = mvrnorm(n, mu = rep(0, p), bdiag(Sigma_U.list))




X = F_%*%L + U
X.test = F_.test%*%L + U.test
gamma = L%*%beta_F
eps = rnorm(n)*sigma^2
eps.test = rnorm(n)*sigma^2

Y = F_%*%gamma + U%*%beta_U + eps
Y.test = F_.test%*%gamma + U.test%*%beta_U + eps.test

X.mean = apply(X, 2, mean)
X = sweep(X, 2, X.mean)
X2U.list = X2U1(X, K = 10, plot = F)

# OLS
data = data.frame(Y=Y, X)
data.test = data.frame(Y=Y.test, X.test)
ml.lm.global = lm(Y~., data = data)
Yhat.lm.global.test = predict(ml.lm.global, new = data.test)
mse.lm.global = mean((Yhat.lm.global.test - Y.test)^2)

# global ridge
ml.ridge.global = cv.glmnet(x=X, y= Y, alpha = 0)
Yhat.ridge.global.test = predict(ml.ridge.global, s=ml.ridge.global$lambda.min, newx = X.test)
mse.ridge.global = mean((Yhat.ridge.global.test - Y.test)^2)

# global EN
ml.EN.global = cv.glmnet(x=X, y= Y, alpha = 0.5)
Yhat.EN.global.test = predict(ml.EN.global, s=ml.EN.global$lambda.min, newx = X.test)
mse.EN.global = mean((Yhat.EN.global.test - Y.test)^2)

# global lasso
ml.lasso.global = cv.glmnet(x=X, y= Y, alpha = 1)
Yhat.lasso.global.test = predict(ml.lasso.global, s=ml.lasso.global$lambda.min, newx = X.test)
mse.lasso.global = mean((Yhat.lasso.global.test - Y.test)^2)

# ALPHA lasso
X2U.list = X2U1(X, K = 10, plot = F)
PY = X2U.list$P%*%Y
data.F = data.frame(Y = PY, X2U.list$F_[,-1])
ml.lm.F = lm(Y~., data = data.F)
HY = X2U.list$H%*%Y
U = X2U.list$H%*%X

L_hat = matrix(X2U.list$L[-1,], ncol = p)
FnU.list = FnU.svd(X.test, L_hat)
F_hat.test = FnU.list$F_
U_hat.test = FnU.list$U

ridge.WLS.U = cv.glmnet(x = U, y = HY, alpha = 0)
EN.WLS.U = cv.glmnet(x = U, y = HY, alpha = 0.5)
lasso.WLS.U = cv.glmnet(x = U, y = HY, alpha = 1)

PYhat.test = F_hat.test%*%ml.lm.F$coefficients

HYhat.test.EN.WLS = predict(EN.WLS.U, s=EN.WLS.U$lambda.min, U_hat.test)
HYhat.test.lasso.WLS = predict(lasso.WLS.U, s=lasso.WLS.U$lambda.min, U_hat.test)
HYhat.test.ridge.WLS = predict(ridge.WLS.U, s=ridge.WLS.U$lambda.min, U_hat.test)

Yhat.test.lasso.WLS = HYhat.test.lasso.WLS + PYhat.test
Yhat.test.EN.WLS = HYhat.test.EN.WLS + PYhat.test
Yhat.test.ridge.WLS = HYhat.test.ridge.WLS + PYhat.test

mse.ridge.WLS = mean((Y.test - Yhat.test.ridge.WLS)^2)
mse.EN.WLS = mean((Y.test - Yhat.test.EN.WLS)^2)
mse.lasso.WLS = mean((Y.test - Yhat.test.lasso.WLS)^2)

