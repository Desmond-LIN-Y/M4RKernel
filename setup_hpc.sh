#!/bin/bash
module load anaconda3/personal
conda activate M4R
conda create -n M4R r-base=4.2.3 -c conda-forge
source activate M4R
conda install r-ragg=1.2.6
Rscript -e "install-dependencies-hpc.R"
conda install conda-forge::r-cmdstanr
conda install conda-forge::r-posterior
conda install conda-forge::r-rlang
conda install conda-forge::r-vctrs
Rscript -e "cmdstanr::install_cmdstan()"
Rscript -e "cmdstanr::check_cmdstan_toolchain()"

