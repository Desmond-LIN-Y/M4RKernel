# To be used in conjunction with hpc_test.sh
# expected outcome: a file containing a string of numbers 
# appearing in appropriate directory

tryCatch({
  library(optparse)
  library(data.table)
  library(yaml)
  library(stringr)
  library(cmdstanr)
  library(tidyr)
  contact_data <- readRDS("/rds/general/user/dl1220/home/M4RKernel/data/COVIMOD_wave_1.rds")
  contact_data <- list(contacts = contact_data$contacts, offsets = contact_data$offsets, population = contact_data$pop)
  config <- yaml::read_yaml("/rds/general/user/dl1220/home/M4RKernel/config/covid_wave1_matern32cd.yml")
  # Unpack data

  dt_contacts <- contact_data$contacts
  dt_offsets <- contact_data$offsets
  dt_population <- contact_data$population
  outpath = "/rds/general/user/dl1220/home/M4RKernel"
  write.csv(dt_population, file.path(outpath, 'output', 'test1.csv'), row.names = F)
  1 = x
  # Path to model
  model_path <- "/rds/general/user/dl1220/home/M4RKernel/stan_models/hsgp-m32-cd.stan"
  
  # Export path
  export_path <- "/rds/general/user/dl1220/home/M4RKernel/output/stan_fits"
  
  
  # Load helpers
  
  source("/rds/general/user/dl1220/home/M4RKernel/helper_func/make_stan_data.R")
  write.csv(dt_population, file.path(outpath, 'output', 'test2.csv'), row.names = F)
  
  
stan_data <- make_stan_data(A = config$data$num_participant_age_groups,
                            C = config$data$num_contact_age_groups,
                            dt_contacts = dt_contacts,
                            dt_offsets = dt_offsets,
                            dt_population = dt_population,
                            model_params = config$model,
                            single_wave = TRUE,
                            single_contact_age = config$data$single_contact_age)
write.csv(stan_data$A, file.path(outpath, 'output', 'test3.csv'), row.names = F)}, 
error = function(e) {
  # Print the error message to the debug file
  debugfile <- file('/rds/general/user/dl1220/home/M4RKernel/output/bug.txt',open='wt')
  cat("\n ERROR: ", conditionMessage(e), "\n", file=debugfile)
}, finally = {
  # Close the debug file, even if there's an error
close(debugfile)
  cat("DONE")
})
