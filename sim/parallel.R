load("~/alpha/data/myseeds0.RData")

for (i in 1:50){
  print(paste0('sbatch -o main.out --mem=4g -t 1:00:00 --wrap="Rscript sim.R myseed=', myseeds[i], '"'))
}
