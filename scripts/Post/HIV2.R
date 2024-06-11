require(data.table)
require(ggplot2) 
require(cmdstanr)
require(bayesplot)
require(posterior)
require(tidyverse)

if (Sys.info()["user"] != "DezLinYiming"){hpc = "/rds/general/user/dl1220/home/M4RKernel/"}else{
  hpc=""
  setwd("D:/Math/Y4/M4R/M4RKernel")
  set_cmdstan_path("E:/cmdstan/cmdstan-2.33.1")}
cat( "Configuring IO ...\n" )
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) == 0) {
  args = c("HIV_2dHSGP_m32_m52-cd", "1")
}
model_name = args[1]

source( paste0(hpc,'helper_func/postprocessing_function.R') )
dc <- readRDS(paste0(hpc, "data/data_age_mixing_within_comm.rds"))
dc = pre_process_dc(dc)

model_fit <- readRDS(paste(hpc, "output/stan_fits/", model_name, ".rds", sep=""))
set.seed(18)


Nmf_all <- nrow(dc[part.sex == 'M' & cont.sex == 'F',])
Nmf <- nrow(dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L,])
Nfm <- nrow(dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L,])
yfm <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, capped_cntcts_in]
ymf <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, capped_cntcts_in]
log_offset_fm <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, log(part) + log(q.term)]
log_offset_mf <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, log(part) + log(q.term)]
ymf_idx <- dc[part.sex == 'M' & cont.sex == 'F' & part > 0L & pop > 0L, age_age_idx]
log_m_mf = model_fit$draws("log_m_mf")
log_mu_mf_pred = log_m_mf[,,ymf_idx] + log_offset_mf
log_m_mf <- NULL
log_under_reported_m_fm = model_fit$draws("log_under_reported_m_fm")
yfm_idx <- dc[part.sex == 'F' & cont.sex == 'M' & part > 0L & pop > 0L, age_age_idx]
log_mu_fm_pred = log_under_reported_m_fm[,,yfm_idx] + log_offset_fm;
log_under_reported_m_fm <- NULL
disp = model_fit$draws(variables='overdispersion')
log_likelihoods1 <- array(NA, dim = c(2000, 4))
log_likelihoods2 <- array(NA, dim = c(2000, 4))
for (i in 1:2000) {
  for (j in 1:4) {
    log_likelihoods1[i, j] <- sum(dnbinom(ymf, mu = exp(log_mu_mf_pred[i, j, ] + 1e-13), size = as.numeric(disp[i, j, 1]), log = TRUE))
    log_likelihoods2[i, j] <- sum(dnbinom(yfm, mu = exp(log_mu_fm_pred[i, j, ] + 1e-13), size = as.numeric(disp[i, j, 1]), log = TRUE))
  }
}



# Compute total ELPD by summing ELPD for mf and fm data
elpd <- mean(log_likelihoods1 + log_likelihoods2)
cat("ELPD:", elpd)
variable_names = load_target_vars_diagnostic()
su = as.data.table(posterior::summarise_draws(
  model_fit$draws(
    variables = variable_names,
    inc_warmup = FALSE)
))
cat("\nThe minimal nESS is: ",su[,min(ess_bulk)], "\nThe maximal rhat is: ",su[,max(rhat)])

# Diagnostics----
diagnostics <- posterior::summarise_draws(posterior::as_draws_df(
  model_fit$sampler_diagnostics()), ~quantile(.x, probs = c(0.01, 0.5, 0.99)
  ))
diagnostics <- rbind(diagnostics, c('min_ess_bulk', su[,min(ess_bulk)],su[,min(ess_bulk)],su[,min(ess_bulk)]))
diagnostics <- rbind(diagnostics, c('max_rhat', su[,max(rhat)],su[,max(rhat)],su[,max(rhat)]))

su.target.vars <- c('y_pred_fm', 'y_pred_mf')
pd <- model_fit$draws(variables = su.target.vars, inc_warmup = FALSE)

select.chains <- seq_along(dimnames(pd)[['chain']])
iters <- seq(from = 1, to = model_fit$metadata()[['iter_sampling']], length.out = ceiling(1e4 / length(select.chains)))
# iters <- seq_len(model_fit$metadata()[['iter_sampling']])
po <- list()


tmp <- pd[,,which(grepl('y_pred_fm',dimnames(pd)[[3]]))]
po$y_pred_fm <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('y_pred_mf',dimnames(pd)[[3]]))]
po$y_pred_mf <- (unname(apply(tmp[iters,select.chains,], 3, rbind)))

gc()

# PPC check
y.fm <- as.data.table(reshape2::melt(po$y_pred_fm))
setnames(y.fm, 1:3, c('iterations', 'id', 'value'))
pds.quantiles <- c(.025,.25,.5,.75,.975)
pds.quantilelabels <- c('CL','IL','M','IU','CU')
# get the ur effects of women
y.fm <- y.fm[,
             list(sd = sd((value)),
                  value_m = quantile((value), p = pds.quantiles, na.rm = TRUE),
                  stat = pds.quantilelabels),
             by = c('id')]

y.fm <- dcast.data.table(y.fm, id + sd~stat, value.var = 'value_m')
tmp2 <- subset(dc, part.sex == 'F' & part > 0 & q.term > 0, select = c(part.sex,part.age,cont.age,part.comm_id,capped_cntcts_in)) #, gamma))
y.fm <- cbind(y.fm, tmp2)
tmp2 <- subset(dc, part.sex == 'F', select = c(part.sex,part.age,cont.age,part.comm_id,part, pop)) #, gamma))
y.fm <- merge(y.fm, tmp2, c('part.sex', 'part.age', 'cont.age', 'part.comm_id'), all = T)
# repeat for mf
y.mf <- as.data.table(reshape2::melt(po$y_pred_mf))
setnames(y.mf, 1:3, c('iterations', 'id', 'value'))

# get the ur effects of women
y.mf <- y.mf[,
             list(sd = sd((value)),
                  value_m = quantile((value), p = pds.quantiles, na.rm = TRUE),
                  stat = pds.quantilelabels),
             by = c('id')]

y.mf <- dcast.data.table(y.mf, id + sd~stat, value.var = 'value_m')
tmp2 <- subset(dc, part.sex == 'M' & part > 0 & q.term > 0, select = c(part.sex,part.age,cont.age,part.comm_id,capped_cntcts_in)) #, gamma))
y.mf <- cbind(y.mf, tmp2)
tmp2 <- subset(dc, part.sex == 'M', select = c(part.sex,part.age,cont.age,part.comm_id,part, pop)) #, gamma))
y.mf <- merge(y.mf, tmp2, c('part.sex', 'part.age', 'cont.age', 'part.comm_id'), all = T)
pds <- rbind(y.fm, y.mf, use.names = T, fill = T)


# PPC check diagnostic ----
# if the pop size is 0, then estimated contacts is 0

set(pds, which(pds$pop == 0), c('CL', 'CU', 'IL', 'IU', 'M'), 0)

pds[, if_out := (capped_cntcts_in < CL | CU < capped_cntcts_in)]

tmp <- mean(pds[, as.integer(if_out)], na.rm = TRUE)
tmp <- 1 - round(tmp, 6)
tmp <- paste0(tmp * 100, '%')
cat("PPC is ", tmp)
