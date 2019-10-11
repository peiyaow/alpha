for (ii in 1:30){
  system(paste0('sbatch -o main_', ii, '.out --mem=2g -t 2:00:00 --wrap="Rscript sim_main.R ii=', ii, '"'))
}