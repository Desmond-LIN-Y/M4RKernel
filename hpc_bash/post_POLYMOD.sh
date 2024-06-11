#!/bin/bash

#PBS -l select=1:ncpus=8:ompthreads=1:mem=80gb
#PBS -l walltime=08:00:00
#PBS -o /rds/general/user/dl1220/home/M4RKernel/output
#PBS -e /rds/general/user/dl1220/home/M4RKernel/output

# Check if the input argument is provided:qw
module load anaconda3/personal
source activate M4R
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/POLYMOD.R "hsgp-m52-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/POLYMOD.R "hsgp-m32-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/POLYMOD.R "hsgp-m12-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/POLYMOD.R "hsgp-eq-rd"
