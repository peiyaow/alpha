set.seed(76)
seeds = floor(runif(80)*1e4)

for (ii in 1:80){
  system(paste0('sbatch -o main.out --mem=4g -t 1:00:00 --wrap="Rscript WLS_svd1.R myseed=', seeds[ii], '"'))
}
