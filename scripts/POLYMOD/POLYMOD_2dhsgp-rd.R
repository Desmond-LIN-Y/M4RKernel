
# Run Stan models on COVIMOD data

cat("\n---------- Compiling and Running Stan Model ----------\n")

# Load libraries
library(optparse)
library(data.table)
library(yaml)
library(stringr)
library(cmdstanr)
library(tidyr)

##### ---------- I/O ---------- #####
if (Sys.info()["user"] != "DezLinYiming"){hpc = "/rds/general/user/dl1220/home/M4RKernel/"}else{
  hpc=""
  setwd("D:/Math/Y4/M4R/M4RKernel")
  set_cmdstan_path("E:/cmdstan/cmdstan-2.33.1")}

cat( "Configuring IO ...\n" )
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  args=c('2dgp-polynomial-cd')
}
# Use the first argument as a string parameter
model_name <- args[1]

cat(model_name)
config <- yaml::read_yaml(paste0(hpc,"config/covid_wave1_2dhsgp_cd.yml"))
model_path <- paste(hpc, "stan_models/", model_name, ".stan", sep="")
export_path <- paste(hpc, "output/stan_fits/POLYMOD/", model_name, ".rds", sep="")

## Configure Stan data
cat(" Configuring Stan data ...\n")
# change upto here
##################################################################################

# Load data
contact_data <- readRDS(paste0(hpc,"data/polymod-stratified-age.rds"))
contact_data <- list(contacts = contact_data$contacts, offsets = contact_data$offsets, population = contact_data$pop)

# Unpack data
dt_contacts <- contact_data$contacts
dt_offsets <- contact_data$offsets
dt_population <- contact_data$population

# Load helpers
source(paste0(hpc, "helper_func/make_stan_data.R"))
source(paste0(hpc, "helper_func/2d_helper.R"))

stan_data <- make_stan_data(A = config$data$num_participant_age_groups,
                            C = config$data$num_contact_age_groups,
                            dt_contacts = dt_contacts,
                            dt_offsets = dt_offsets,
                            dt_population = dt_population,
                            model_params = config$model,
                            single_wave = TRUE,
                            single_contact_age = config$data$single_contact_age)
stan_data$Xhsgp <- scale(make_nn_lowertri_idxset(85))
stan_data$S = make_S_matrix(30,30) # fix M
stan_data$sym_from_lowertri_idxset <- make_sym_from_lowertri_idx(85)
stan_data$A2 = 3655

stan_data$Xhsgp2 = scale(make_nn_idxset(85))

cat(" Compiling Stan model ...\n")
# Compile stan program
model <- cmdstanr::cmdstan_model(model_path, force_recompile = TRUE)

cat(" Running Stan model ...\n")

model_params <- config$model
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = model_params$warmup,
  iter_sampling = model_params$sampling,
  max_treedepth = model_params$max_treedepth,
  adapt_delta = model_params$adapt_delta,
  refresh = model_params$refresh,
  seed = 1234
)

cat(" Saving fitted model ...\n")
fit$save_object(file =export_path)

cat("\n Run Stan ALL DONE.\n")