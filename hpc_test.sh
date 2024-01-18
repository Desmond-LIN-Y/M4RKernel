#!/bin/bash
########## EDIT THIS SECTION ##########
REPO_PATH=/rds/general/user/dl1220/home/M4RKernel
OUT_PATH=/rds/general/user/dl1220/home/M4RKernel/output
#######################################

# Create main script
cat > "$OUT_PATH/hpc_test.pbs" <<EOF
#!/bin/bash

#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=4:ompthreads=1:mem=50gb
#PBS -J 1-10

module load anaconda3/personal
source activate M4R

cd "$REPO_PATH"
Rscript hpc_test.R
EOF
cd $OUT_PATH
qsub hpc_test.pbs