n = 10^seq(3, 8, by = 0.25)
gamma = 1.2
p = n^gamma
m = 10
ng = n/m


# groupwise
s_L = .1
s_d = .1
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
plot(F_risk_group_ridge, type="l", ylab = "risk", xlab = "n", lty = 2)

F_risk = 1/ng + 1/ng^{1/4}/sqrt(p)
lines(F_risk)


n = 10^seq(3, 8, by = 0.25)
gamma = 1.2
p = n^gamma
m = 10
ng = n/m

# groupwise case 1
s_L = 0.1
s_d = 0.1
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
plot(F_risk_group_ridge, type="l", ylim = c(0, 1.5), ylab = "risk", xlab = "n", lty = 2)

# groupwise case 2
s_L = .01
s_d = .8
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
lines(F_risk_group_ridge, type="l", ylab = "risk", xlab = "n", lty = 2)

# groupwise case 3
s_L = .1
s_d = .1
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
lines(F_risk_group_ridge, type="l", ylab = "risk", xlab = "n", lty = 2)

# my method 
mp = 0
mp = p^mp
wn = 1/sqrt(p) + sqrt(log(p)/n)
U_risk_ridge = sqrt(p*log(p))/n + (p*log(p)/n^(3/2))*(1 + mp*wn)
lines(U_risk_ridge)


n = 10^seq(3, 8, by = 0.25)
gamma = 1.1
p = n^gamma
m = 10
ng = n/m

# groupwise
s_L = .1
s_d = .5
Lg = sqrt(p)^s_L
dg = sqrt(p)^s_d

# ridge
F_risk_group_ridge = 1/sqrt(ng)*Lg*dg + 1/ng*Lg + log(p)/ng
plot(F_risk_group_ridge, type="l", ylim = c(0, 1.5),ylab = "risk", xlab = "n", lty = 2)

mp = 0
mp = p^mp
wn = 1/sqrt(p) + sqrt(log(p)/n)
U_risk_ridge = sqrt(p*log(p))/n + (p*log(p)/n^(3/2))*(1 + mp*wn)
lines(U_risk_ridge)
