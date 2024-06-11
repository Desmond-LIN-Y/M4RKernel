#!/bin/bash

#PBS -l select=1:ncpus=8:ompthreads=1:mem=80gb
#PBS -l walltime=71:00:00
#PBS -o /rds/general/user/dl1220/home/M4RKernel/output
#PBS -e /rds/general/user/dl1220/home/M4RKernel/output

module load anaconda3/personal
source activate M4R

Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/HIV/HIV_2dhsgp_rd.R "HIV_2dHSGP_m32_m52"


