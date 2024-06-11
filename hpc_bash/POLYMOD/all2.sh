#!/bin/bash

#PBS -l select=1:ncpus=8:ompthreads=1:mem=80gb
#PBS -l walltime=72:00:00
#PBS -o /rds/general/user/dl1220/home/M4RKernel/output
#PBS -e /rds/general/user/dl1220/home/M4RKernel/output

module load anaconda3/personal
source activate M4R

Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-eq-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m12Neq"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m12Nm32"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m12-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m12Xeq"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m12Xm32"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m32-rd"
Rscript /rds/general/user/dl1220/home/M4RKernel/scripts/POLYMOD/POLYMOD_2dhsgp-rd.R "2dhsgp-m32Neq"