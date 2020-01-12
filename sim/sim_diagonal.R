Sigma_U_hat = POET(t(U.train), K=0, C= 0.5, matrix = "cor")$SigmaU

L.list[[2]]%*%solve(300*Sigma_U_hat)%*%t(U.train)%*%HY.train
L.list[[3]]%*%solve(300*diag(p))%*%t(U.train)%*%HY.train

F = P.list[[2]]%*%X.train.list[[2]]
Ut = t(X.train.list[[2]])%*%(diag(100) - P.list[[2]])

F = F.train.list[[1]]
X.train.list[[1]]%*%t(X.train.list[[1]])%*%F/100/F[,2]
t(X.train.list[[1]])[,]%*%F[,]/100/F[,2]
P_big = bdiag(P.list)
X_big = do.call(rbind, X.train.list)
F_big = P_big%*%X_big
U_big = (diag(300) - P_big)%*%X_big

haha = F_big%*%t(U_big)
t(U_big)%*%F_big

haha[1:100, 101:200]




# svd
SVD.res = svd(X.train.list[[1]]/sqrt(100), nu = 100, nv = 80)
PCA.res = eigen(X.train.list[[1]]%*%t(X.train.list[[1]])/100)
K = 3
new_diag = c(rep(1, K) , SVD.res$d[(K+1):p]^(-1)*SVD.res$d[K+1])
A = SVD.res$v%*%diag(new_diag, nrow = p, ncol = p)
diag(t(A)%*%A)
new_X = X.train.list[[1]]%*%A
SVD.newres = svd(new_X/sqrt(100))
SVD.newres$d

F_ = SVD.newres$u[,1:K]*sqrt(100)
t(F_)%*%F_
P = 1/n*F_%*%t(F_)
U = (diag(n) - P)%*%new_X
t(U)%*%U/100

F_ = SVD.res$u[,1:K]*sqrt(100)
t(F_)%*%F_
P = 1/n*F_%*%t(F_)
U = (diag(n) - P)%*%X.train.list[[1]]
t(U)%*%U/100 + lam*I

svd(t(U)%*%U/100)
c = 1
I = diag(p)
lam = 1
PCA.res = eigen(t(U)%*%U/100 + lam*I)
A = PCA.res$vectors%*%diag(1/sqrt(PCA.res$values))
t(A)%*%(t(U)%*%U/100 + lam*I)%*%(A)

t(SVD.res$v)%*%(t(U)%*%U/100 + c*I)%*%SVD.res$v


U_new = (SVD.res$u)%*%diag(c(rep(0,K), sqrt(100)*SVD.res$d[(K+1):p]), ncol = p, nrow = n)%*%t(SVD.res$v)

SVD.res$v%*%t(SVD.res$v)

t(U_new)%*%U_new

SVD.res$v%*%t(diag(c(rep(0,K), sqrt(100)*SVD.res$d[(K+1):p]), ncol = p, nrow = n))%*%diag(c(rep(0,K), sqrt(100)*SVD.res$d[(K+1):p]), ncol = p, nrow = n)

t(diag(c(rep(0,K), sqrt(100)*SVD.res$d[(K+1):p]), ncol = p, nrow = n))%*%diag(c(rep(0,K), sqrt(100)*SVD.res$d[(K+1):p]), ncol = p, nrow = n)


%*%t(SVD.res$v)


SVD.res$v%*%t(SVD.res$v)



a = svd(t(U)%*%U/100 + c*I)$v
b = svd(t(U)%*%U/100)$v




dim(SVD.res$v)
SVD.res$v%*%t(SVD.res$v)
dim(SVD.res$u)
t(SVD.res$u)%*%SVD.res$u[,1:3]
t(PCA.res$vectors)%*%PCA.res$vectors

SVD.res$d^2
PCA.res$values


haha = F%*%Ut
Ut%*%t(Ut)
Ut[,1]%*%t(F[1,])

Sigma_U_hat = POET(t(U.train.list[[2]]), K=0, C= 0.5, matrix = "vad")$SigmaU
L.list[[2]]%*%solve(100*Sigma_U_hat)%*%t(U.train.list[[2]])%*%Y.train.list[[2]]
%*%Y.train.list[[2]]
L.list[[2]]%*%solve(100*diag(p))%*%t(U.train.list[[2]])%*%Y.train.list[[2]]




