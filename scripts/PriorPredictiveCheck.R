library(optparse)
library(data.table)
library(yaml)
library(stringr)
library(cmdstanr)
library(tidyr)
if (Sys.info()["user"] != "DezLinYiming"){hpc = "/rds/general/user/dl1220/home/M4RKernel/"}else{
  hpc=""
  setwd("D:/Math/Y4/M4R/M4RKernel")
  set_cmdstan_path("C:/Users/DezLinYiming/Documents/.cmdstan/cmdstan-2.33.1")}

source("helper_func/2d_helper.R")
ind = make_S_matrix(30,30)
stan_data = list(A=85, age_idx_std = xxx,C1=1.5,C2=1.5,M1=30,M2=30,indices=ind)
model <- cmdstanr::cmdstan_model("stan_models/prior_predictive_check.stan", force_recompile = TRUE)

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains =4,
  iter_warmup = 500,
  iter_sampling = 1000,
  #max_treedepth = 13,
  adapt_delta = 0.9,
  refresh = 500,
  seed = 1,
  fixed_param = TRUE
)

samples = fit$draws(variables='f')
f11 =  as.array(samples)[, 1, 1]
