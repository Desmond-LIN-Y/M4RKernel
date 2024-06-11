# Preamble: Generates diagnostic statistics and result plots
library(optparse)
library(data.table)
library(yaml)
library(cmdstanr)
library(bayesplot)
library(loo)
library(reshape2)
library(stringr)
library(ggplot2)
library(viridis)
library(pammtools)
cat("\n ---------- Begin Post-processing ---------- \n")
if (Sys.info()["user"] != "DezLinYiming"){hpc = "/rds/general/user/dl1220/home/M4RKernel/"}else{
  hpc=""
  setwd("D:/Math/Y4/M4R/M4RKernel")
  set_cmdstan_path("E:/cmdstan/cmdstan-2.33.1")}


args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  args=c('polynomial')
}


# Use the first argument as a string parameter
input_arg <- args[1]

export_path = paste0(paste0(paste0(hpc, "output/post/"), input_arg), "/")


# Print the input argument to verify it's being read correctly
cat(paste("The input argument is:", input_arg))


bayesplot::color_scheme_set(scheme = "mix-blue-pink")

model_path = paste0(paste0(hpc, paste0("output/stan_fits/", input_arg)), '.rds')
data_path = paste0(hpc, "data/polymod-stratified-age.rds")

##### ---------- Setup ---------- #####
cat(paste("\n Model path:", model_path))
cat(paste("\n Data path:", data_path), "\n\n")

fit <- readRDS(model_path)
data <- readRDS(data_path)

# Unpack data
dt_contacts <- data$contacts
dt_offsets <- data$offsets
dt_population <- data$population

# Load helpers
source(paste0(hpc, "helper_func/convergence_diagnostic_stats.R"))
source(paste0(hpc,"helper_func/posterior_predictive_check.R"))
source(paste0(hpc,"helper_func/posterior_log_contact_rates.R"))
source(paste0(hpc,"helper_func/posterior_contact_intensity.R"))


posterior_draws = readRDS(paste0(export_path, "posterior_draws.rds"))
dt_posterior <- posterior_log_contact_rates(posterior_draws)
saveRDS(dt_posterior, file = paste0(export_path, "dt_posterior.rds"))
dt_matrix <- posterior_contact_intensity(dt_posterior,
                                         dt_population,
                                         type = "matrix",
                                         outdir = export_path)

dt_margin <- posterior_contact_intensity(dt_posterior,
                                         dt_population,
                                         type = "marginal",
                                         outdir = export_path)
