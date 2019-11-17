load("/nas/longleaf/home/peiyao/alpha/data/myseeds.RData")
nsim = length(myseeds)

for (i in 1:nsim){
  system(paste0('sbatch -o main.out --mem=2g -t 1:00:00 --wrap="Rscript krr_X2U_longleaf_method2_single.R myseed=', myseeds[i], '"'))
}