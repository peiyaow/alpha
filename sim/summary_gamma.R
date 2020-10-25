result.gamma = rbind(apply(sqrt(do.call(rbind, DIFF.gamma)), 2, mean), 
              apply(sqrt(do.call(rbind, DIFF.gamma)), 2, sd)/sqrt(10))
row.names(result.gamma) = c("mean", "se")
apply(cbind(gamma1, gamma2, gamma3), 2, function(x) sqrt(mean(x^2)))
apply(cbind(gamma1, gamma2, gamma3), 2, function(x) sd(x))


rbind(apply(sqrt(do.call(rbind, DIFF.beta)), 2, mean), 
      apply(sqrt(do.call(rbind, DIFF.beta)), 2, sd)/sqrt(10))
sqrt(mean(beta_share^2))
sd(beta_share)

result.MSE = apply(array(unlist(MSE), c(5, 3, 10)), c(1,2), mean)
row.names(result.MSE) = c("global", "group", "alpha0", "alpha10", "alpha")
apply(array(unlist(MSE), c(5, 3, 10)), c(1,2), sd)/sqrt(10)
