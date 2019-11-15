#our method
# ridge
n = 10^seq(3, 8, by = 0.5)
gamma = 0.7
p = n^gamma

m = 10
ng = n/m
F_risk = 1/ng + 1/ng^{1/4}/sqrt(p)
plot(F_risk, ylim = c(0, .1), type = "l")

mp = 0
mp = p^mp
wn = 1/sqrt(p) + sqrt(log(p)/n)
U_risk_ridge = sqrt(p*log(p))/n + (p*log(p)/n^(3/2))*(1 + mp*wn)
lines(U_risk_ridge)

# lasso
b1 = 1
rs = 0
s = p^rs
U_risk_lasso = (b1+1)*sqrt(p)*log(p)/n^{3/2}*sqrt(s) + sqrt(p*log(p))/n
lines(U_risk_lasso)

# groupwise
s_L = .1
s_d = .1
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
lines(F_risk_group_ridge)

# lasso
F_risk_group_lasso = sqrt(s)*sqrt(p/ng)*((Lg+sqrt(log(p)/ng))*dg + Lg/sqrt(ng) + sqrt(log(p)/ng))
lines(F_risk_group_lasso)

# global
# ridge
s_L = 0.5
s_d = 0.35
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d
F_risk_global_ridge = 1/sqrt(n)*Lg*dg + 1/n*Lg + log(p)/n
lines(F_risk_global_ridge)



