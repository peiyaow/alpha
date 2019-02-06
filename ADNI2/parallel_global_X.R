nsim = 50
load("/nas/longleaf/home/peiyao/alpha/data/myseeds.RData")

for (i in 1:nsim){
  system(paste0('sbatch -o main.out --mem=2g -t 1:00:00 --wrap="Rscript krr_X2U1_longleaf.R myseed=', myseeds[i], '"'))
}