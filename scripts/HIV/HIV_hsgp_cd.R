# This is the script having equal 15-49 age range of participants and 15-69 aged partners
# we account for the age heter for female reported data, accounting for the under-reporting situation
# we did not apply the gender heter in this script

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)
library(cmdstanr)
#
# Define input arguments that can be changed by users

# # Define input arguments that can be changed by users----
# args <- list()
# args$seed <- 18L

if (Sys.info()["user"] != "DezLinYiming"){hpc = "/rds/general/user/dl1220/home/M4RKernel/"
}else{
  hpc=""
  setwd("D:/Math/Y4/M4R/M4RKernel")
  set_cmdstan_path("E:/cmdstan/cmdstan-2.33.1")}

account_missing = TRUE
community_name = "fishing" # 
source(paste0(hpc,'helper_func/GP-functions.R') )
source(paste0(hpc, "helper_func/2d_helper.R"))

args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  args=c('HIV_HSGP_con_cd')
}
model_name = args[1]
cat(model_name)
# change here!
model_path = paste(hpc, 'stan_models/', model_name, '.stan', sep="")
model <- cmdstanr::cmdstan_model(model_path, force_recompile = TRUE)
export_path <- paste(hpc, "output/stan_fits/", model_name, ".rds", sep = "")
#Change up to here
########################################################################################################
# load data----
cat("\n Loading the partnership data at individual level...")
# part.age: participant age
# part.sex: participant sex

# part.T: Census eligible individuals in the age/sex/round/comm of the participant
# Z number of sexual intercourse that the participant have had in the part year
# parter_age_X: Reported age of the partner in the Xth sexual intercourse that occurred in the past year.
# parter's age ranged from 15 to 69

# part.round: participant round
# part.comm: participant area = inland, fishing
# part.comm_id: participant community

# partner_living_same_household_1-4: Was the partner living in the same household (from rltnhh1) - can be NA, YES or NO.
# if the partner does not live in the same household as participants, then they will be asked about the community:
# partner_living_same_community_1-4: Was the partner living in the same community (from rltncm1) - can be NA, YES, NO or DK (don't know).
# sexp1out: Number of partners in the last 12 months outside the community - coded as Z, no NA
# part.hiv: Did the participant have HIV - can be NA (for 29 participants), N (negative) or P (positive)
#

dc <- as.data.table(readRDS(paste0(hpc,'data/data_age_mixing_within_comm.rds')))
dc <- dc[part.comm == community_name & part.round == "R019"]
str(dc)
tmp <- dc[, list(cp = sum(capped_cntcts, na.rm = TRUE)), by = c('part.sex', 'part.comm_id', 'part.round')]
tmp <- subset(tmp, cp == 0L, select = c('part.sex', 'part.comm_id'))
dc[, check := paste0(part.sex, '-', part.comm_id)]
tmp[, check := paste0(part.sex, '-', part.comm_id)]
dc <- dc[!(check %in% tmp$check)]
set(dc, NULL, 'check', NULL)
setkey(dc, part.comm_id, cont.sex, part.sex, cont.age, part.age)
tmp <- as.data.table(1:length(unique(dc$part.comm_id)))
tmp$part.comm_id <- unique(dc$part.comm_id)
setnames(tmp, 'V1', 'c')
dc <- merge(dc, tmp, by = 'part.comm_id')
if (account_missing)
{
  # apply rho or q.term (Q) to account for the missing reports
  tmp <- dc[,
            list(q.term = sum(capped_cntcts, na.rm = TRUE) / unique(sexp1in)),
            by = c("part.sex", "part.age", 'part.comm_id', 'part.round')
  ]
  dc <- merge(dc, tmp, by = c("part.sex","part.age", 'part.comm_id', 'part.round'), all.x = TRUE)
  tmp <- dc[,
            list(capped_cntcts.ma = sum(capped_cntcts, na.rm = TRUE) ,
                 sexp1in = unique(sexp1in)),
            by = c("part.sex", "part.age", 'part.comm_id', 'part.round')
  ]
  tmp <- tmp[, list(ma.cntcts.t = sum(capped_cntcts.ma, na.rm = TRUE),
                    sexp1in.t = sum(sexp1in, na.rm = TRUE)),
             by = c('part.age', 'part.sex')]
  tmp[, Q.term := ma.cntcts.t/sexp1in.t]
  dc <- merge(dc, tmp, by = c("part.sex","part.age"), all.x = TRUE)
  # set(dc, which(is.na(dc$q.term)), 'q.term', dc[is.na(q.term), Q.term])
  # if capped_cntcts == 0 while sexp1in !=0, we cannot scale the estimations, so we use Q.term
  
  # impute q.term
  set(dc, which(dc$ma.cntcts.t == 0 & dc$sexp1in.t != 0), 'q.term', dc[ma.cntcts.t == 0 & sexp1in.t != 0, Q.term])
  # set(dc, which(dc$q.term == 0), 'q.term', dc[q.term == 0, Q.term])
  # set(dc, which(is.na(dc$q.term) & dc$sexp1in == 0), 'q.term', 1)
  summary(dc$q.term)
  dc[q.term %in% c(Inf, 1, 0), q.term := Q.term]
  dc[is.na(q.term), q.term := Q.term]
  dc[is.na(q.term), q.term := 1]
  
  set(dc, NULL, c('ma.cntcts.t', 'sexp1in.t'), NULL)
}
tmp <- as.data.table(expand.grid(part.sex = c("M"), part.age = unique(dc$part.age), cont.age = unique(dc$cont.age), part.comm_id = unique(dc$part.comm_id)))
tmp[, cont.sex := 'F']
setkey(tmp,part.comm_id, cont.sex, part.sex, cont.age, part.age)

tmp[, fm_mat_idx := seq_len(length(part.sex)), by = c('part.comm_id')]
tmp2 <- copy(tmp)
setnames(tmp2, c('cont.age','part.age','cont.sex','part.sex'), c('part.age','cont.age','part.sex','cont.sex'))
tmp <- rbind(tmp, tmp2)

dc <- merge(dc, tmp, by = c('cont.sex','part.sex','cont.age','part.age', 'part.comm_id'), all.x = TRUE)
setkey(dc,part.comm_id, cont.sex, part.sex, cont.age, part.age)

# age-age idx related to community: suitable for long vector input data, i.e. age heaping groups
dc[, age_age_idx := seq_len(length(part.age)), by = c('part.sex')]
setnames(dc, 'capped_cntcts', 'capped_cntcts_in')


stan_data <- list()
stan_data$A1 <- length(unique(dc$part.age))
stan_data$A2 <- length(unique(dc$cont.age))
# stan_data$age1 <- unique(dc$part.age)
# stan_data$age2 <- unique(dc$cont.age)

# community label
stan_data$Nc <- length(unique(dc$c))
stan_data$community_index_mf_all <- dc[part.sex == 'M' & cont.sex == 'F', c]
stan_data$community_index_fm_all <- dc[part.sex == 'F' & cont.sex == 'M', c]

# contact data
stan_data$Nmf <- nrow(dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L,])
stan_data$Nmf_all <- nrow(dc[part.sex == 'M' & cont.sex == 'F',])
stan_data$ymf <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, capped_cntcts_in]
stan_data$ymf_rowmajor_matrix_index_all <- dc[part.sex == 'M' & cont.sex == 'F', fm_mat_idx]
stan_data$ymf_idx <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, age_age_idx]

stan_data$Nfm <- nrow(dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L,])
stan_data$yfm <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, capped_cntcts_in]
stan_data$yfm_rowmajor_matrix_index_all <- dc[part.sex == 'F' & cont.sex == 'M', fm_mat_idx]
stan_data$yfm_idx <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, age_age_idx]

# offsets
stan_data$log_offset_mf <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, log(part) + log(q.term)]
stan_data$log_offset_fm <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, log(part) + log(q.term)]

# for age heaping
dc[pop == 0, pop := 1]
stan_data$log_pop_mf_all <- dc[part.sex == 'M' & cont.sex == 'F', log(pop)]
stan_data$log_pop_fm_all <- dc[part.sex == 'F' & cont.sex == 'M', log(pop)]
stan_data$age_std = scale(1:55)[,1]

# HSGP input variables
stan_data$c_age1 <- 1.5
stan_data$c_age2 <- 1.5
stan_data$M_age1 <- 30
stan_data$M_age2 <- 30
stan_data$index_mf <- unique(dc[part.sex == 'M', list(part.age,cont.age,fm_mat_idx)])$fm_mat_idx
stan_data$index_fm <- unique(dc[part.sex == 'F', list(part.age,cont.age,fm_mat_idx)])$fm_mat_idx

# standardise the age with the full combination


indices <- matrix(NA, 30 * 30, 2)
mm <- 0;
for (r in 1:30){
  for (s in 1:30){
    mm <- mm + 1
    indices[mm,] <- c(r, s)
  }
}
stan_data$S <- indices
model_inits <- list(
  list(
    log_random_effect_community = rnorm(stan_data$Nc, -9.5, 0.2),
    log_effect_under_report_baseline_f = rnorm(1, -0.5, 0.2)
  ),
  list(log_random_effect_community = rnorm(stan_data$Nc, -9.5, 0.2),
       log_effect_under_report_baseline_f = rnorm(1, -0.5, 0.2)
  ),
  list(log_random_effect_community = rnorm(stan_data$Nc, -9.5, 0.2),
       log_effect_under_report_baseline_f = rnorm(1, -0.5, 0.2)
  ),
  list(log_random_effect_community = rnorm(stan_data$Nc, -9.5, 0.2),
       log_effect_under_report_baseline_f = rnorm(1, -0.5, 0.2)
  )
)





# Run Stan model---

cat("\n Run Stan model ...")

model_fit <- model$sample(
  data = stan_data,
  seed = 18,
  chains = 4,
  parallel_chains = 4,
  refresh = 300,
  iter_warmup = 300,
  iter_sampling = 2000,
  
  #iter_warmup = 1e3,
  #iter_sampling = 5e3,
  max_treedepth = 12,
  # adapt_delta = 0.95,
  # 4 chains in parallel with 4 threads each
  # threads_per_chain = 4,
  save_warmup = TRUE,
  init = model_inits
)


cat("\n Save fitted data to file \n")
model_fit$save_object(file = export_path)
cat("\nDone\n.")
