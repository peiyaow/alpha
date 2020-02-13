nsim = 100
load("/nas/longleaf/home/peiyao/alpha/data/myseeds.RData")

for (i in 1:nsim){
  system(paste0('sbatch -o main.out --mem=2g -t 1:00:00 --wrap="Rscript real_data.R myseed=', myseeds[i], '"'))
}