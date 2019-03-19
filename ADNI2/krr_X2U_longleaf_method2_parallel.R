load("/nas/longleaf/home/peiyao/alpha/data/myseeds.RData")
nsim = length(myseeds)

# for (i in 1:nsim){
#   system(paste0('sbatch -o main.out --mem=2g -t 1:00:00 --wrap="Rscript krr_X2U_longleaf_method2.R myseed=', myseeds[i], '"'))
# }

for (i in 1:nsim){
  for (threshold in c(0.1, 0.15, 0.2, 0.25, 0.3)){
    system(paste0('sbatch -o main.out --mem=2g -t 1:00:00 --wrap="Rscript krr_X2U_longleaf_method2.R myseed=', myseeds[i], ' threshold=', threshold,'"'))
  }
}



