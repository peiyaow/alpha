p = 200
n = 100
p_share = 197
K = 3
n_label = 3
si = 5

# control parameter for beta
s = 5 # unique
h = 0.5 # share

ds = .3 # dilute parameter for spike
rho = 1/5
rrho = 0.01
d = 1/2

spike1 = c(20, 10, 10/3)*ds
spike2 = c(15, 7.5, 2.5)*ds
spike3 = c(10, 5, 5/3)*ds

c(sqrt(spike1), d)

1/sqrt(spike1)
1/sqrt(spike2)
1/sqrt(spike3)



sigma.vec = rep(si, n_label)

beta1_unique = c(0.5, 1, 1)*s
beta2_unique = c(1, 0.5, 1)*s
beta3_unique = c(1, 1, 0.5)*s

beta_share = rep(h, p_share)

beta1 = c(beta1_unique, beta_share, rep(0, p-p_share-3))
beta2 = c(beta2_unique, beta_share, rep(0, p-p_share-3))
beta3 = c(beta3_unique, beta_share, rep(0, p-p_share-3))

n.train.vec = c(n, n, n)
n.test.vec = c(n, n, n)
ix.vec = c(0, cumsum(n.train.vec))
label.test = as.factor(c(rep(1, n.test.vec[1]), rep(2, n.test.vec[2]), rep(3, n.test.vec[3])))
label.level = levels(label.test)

para1 = FactorModelPara(n, p, K, spike = spike1, d, du = spike1[K]*rho, rrho)
para2 = FactorModelPara(n, p, K, spike = spike2, d, du = spike2[K]*rho, rrho)
para3 = FactorModelPara(n, p, K, spike = spike3, d, du = spike3[K]*rho, rrho)

t(beta1)%*%t(para1$L)%*%para1$L%*%beta1
t(beta1)%*%para1$SigmaU%*%beta1

t(beta2)%*%t(para2$L)%*%para2$L%*%beta2
t(beta2)%*%para2$SigmaU%*%beta2

t(beta3)%*%t(para3$L)%*%para3$L%*%beta3
t(beta3)%*%para3$SigmaU%*%beta3

diag(para1$SigmaX)
