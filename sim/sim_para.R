setwd("~/Documents/GitHub/alpha/sim/parameter/")
# save parameter
beta_U.grid = c(0, 1, 3)
sigma.vec.grid = list(c(1, 1, 1), c(1, 2, 3))
p.grid = 80*c(1, 2, 3, 4, 5)
para_grid.list = list()
i = 1
for (sigma.vec in sigma.vec.grid){
  for (p in p.grid){
    for (beta in beta_U.grid){
      para_grid.list[[i]] = list(beta = beta, sigma.vec = sigma.vec, p = p)
      i = i+1
    }
  }
}
save(para_grid.list, file = "para_grid.RData")

seeds = floor(runif(50)*10000)
save(seeds, file = "seeds.RData")
