#!/bin/bash

#PBS -l select=1:ncpus=8:ompthreads=1:mem=80gb
#PBS -l walltime=08:00:00

#PBS -o /rds/general/user/dl1220/home/M4RKernel/output
#PBS -e /rds/general/user/dl1220/home/M4RKernel/output

module load anaconda3/personal
source activate M4R

Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/COVIMOD/COVIMOD_wave1_2dhsgp_m32_cd_vec.R "1234"


