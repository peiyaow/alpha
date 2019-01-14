nsim = 50
load("/nas/longleaf/home/peiyao/alpha/data/myseeds.RData")

for (i in 1:nsim){
  for (gamma in 10^(-3:0)){
    system(paste0('sbatch -o main_', -log10(gamma), '.out -t 5:00:00 --mem=4g --wrap="Rscript krr_U_longleaf.R gamma=', gamma, ' myseed=', myseeds[i], '"'))
  }
}