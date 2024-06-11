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
 args = c("2dhsgp-eq-rd")
}


# Use the first argument as a string parameter
input_arg <- args[1]

export_path = paste0(paste0(paste0(hpc, "output/post/POLYMOD"), input_arg), "/")

if (dir.exists(export_path)) {
    cat("Directory exists:", export_path, "\n")
} else {
    cat("Create directory:", export_path, "\n")
    dir.create(export_path)
}
# Print the input argument to verify it's being read correctly
cat(paste("The input argument is:", input_arg))


bayesplot::color_scheme_set(scheme = "mix-blue-pink")

model_path = paste0(paste0(hpc, paste0("output/stan_fits/POLYMOD/", input_arg)), ".rds")
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

##### ---------- Assess convergence and mixing ---------- #####
cat(" Assess convergence and mixing ...\n")

# Make convergence diagnostic tables

fit_summary <- convergence_diagnostic_stats(fit, outdir = export_path)

# Make trace plots
cat("\n Making trace plots")

pars <- c('nu', 'gp_alpha', 'gp_rho_1', 'gp_rho_2')
pars_po <- fit$draws(pars)
p <- bayesplot::mcmc_trace(pars_po,
                           facet_args = list(nrow = 10, ncol = 1))
ggsave(file = file.path(export_path, 'mcmc_trace_parameters.png'),
       plot = p,
       h = 20,
       w = 20,
       limitsize = F)

# Make pairs plots
cat(" Making pairs plots\n")
p <- bayesplot::mcmc_pairs(pars_po,
                           off_diag_args = list(size = 0.3, alpha = 0.3))
ggsave(file = file.path(export_path, 'mcmc_pairs_parameters.png'),
       plot = p,
       h = 20,
       w = 20,
       limitsize = F)

###### ---------- Posterior predictive checks ---------- #####
cat(" Extracting posterior ...\n")
posterior_draws <- fit$draws(c("yhat_strata", "log_cnt_rate"),
                             inc_warmup = FALSE,
                             format = "draws_matrix")
saveRDS(posterior_draws, file = paste0(export_path, "posterior_draws.rds"))
cat(" Making posterior predictive checks ...\n")

dt_ppc <- posterior_predictive_check(posterior_draws,
                                     dt_contacts,
                                     single_contact_age = FALSE,
                                     outdir = export_path)

cat(" Extracting posterior contact intensities ...\n")

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


cat("\n DONE.\n")


