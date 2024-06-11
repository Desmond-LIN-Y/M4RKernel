
#!/bin/bash

#PBS -l select=1:ncpus=4:ompthreads=1:mem=50gb
#PBS -l walltime=08:00:00

module load anaconda3/personal
source activate M4R

Rscript scripts/hpc_test.R

