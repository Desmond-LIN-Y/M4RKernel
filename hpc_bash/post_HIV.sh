#!/bin/bash

#PBS -l select=1:ncpus=8:ompthreads=1:mem=80gb
#PBS -l walltime=08:00:00
#PBS -o /rds/general/user/dl1220/home/M4RKernel/output
#PBS -e /rds/general/user/dl1220/home/M4RKernel/output

# Check if the input argument is provided:qw
module load anaconda3/personal
source activate M4R


# Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/HIV_NEED_IN.R "HIV_2dHSGP_con_m52-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/HIV_NEED_IN.R "HIV_HSGP_con_rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/Post/HIV_NEED_IN.R "HIV_HSGP_con_cd"
