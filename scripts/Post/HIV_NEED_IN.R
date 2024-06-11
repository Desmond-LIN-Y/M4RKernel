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
args = c("HIV_HSGP_m52_rd", "1")
}
model_name = args[1]
source( paste0(hpc,'helper_func/plot_functions.R') )
source( paste0(hpc,'helper_func/postprocessing_function.R') )
dc <- readRDS(paste0(hpc, "data/data_age_mixing_within_comm.rds"))
dc = pre_process_dc(dc)

model_fit <- readRDS(paste(hpc, "output/stan_fits/",model_name, ".rds", sep=""))
set.seed(18)



# Extract Monte Carlo samples for postprocessing----
cat("\nExtract Monte Carlo samples for postprocessing ...\n")
# get the posterior median of the under-reporting effects of women
# for the external community sexual network computation

  # log_m_mf are the estimated log of contact intensities
  su.target.vars <- load_target_vars_posterior()
  pd <- model_fit$draws(variables = su.target.vars, inc_warmup = FALSE)

  select.chains <- seq_along(dimnames(pd)[['chain']])
  iters <- seq(from = 1, to = model_fit$metadata()[['iter_sampling']], length.out = ceiling(1e4 / length(select.chains)))
  # iters <- seq_len(model_fit$metadata()[['iter_sampling']])
  po <- list()


tmp <- pd[,,which(grepl('log_m_mf',dimnames(pd)[[3]]))]
po$log_m_mf <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_m_fm',dimnames(pd)[[3]]))]
po$log_m_fm <- (unname(apply(tmp[iters,select.chains,], 3, rbind)))
# extract the 2dgp
tmp <- pd[,,which(grepl('f_fm',dimnames(pd)[[3]]))]
po$log_effect_rf <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_effect_age_under_report_f',dimnames(pd)[[3]]))]
po$log_effect_age_under_report_f <- unname(apply(tmp[iters,select.chains,], 3, rbind))

pd <- NULL
gc()

pds.quantiles <- c(.025,.25,.5,.75,.975)
pds.quantilelabels <- c('CL','IL','M','IU','CU')

# Viz the random functions ----
pd.rf <- as.data.table(reshape2::melt(po$log_effect_rf))
setnames(pd.rf, 1:3, c('iterations', 'fm_mat_idx', 'value'))
tmp <- pd.rf[,
                     list(sd = sd((value)),
                          value_m = quantile((value), p = pds.quantiles, na.rm = TRUE),
                          stat = pds.quantilelabels),
                     by = c('fm_mat_idx')]
tmp <- dcast.data.table(tmp, fm_mat_idx + sd~stat, value.var = 'value_m')
tp <- merge(tmp, unique(dc[part.sex == 'M', list(part.age,cont.age,part.sex,cont.sex,fm_mat_idx)]),
            by = c('fm_mat_idx'), all = T)
write.csv(tp, file = paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_posterior_within_commu_rf.csv", sep=""), row.names = F)


# Extract age and gender specific contact intensities m_ab^{cgh} ----
cat("\nExtract the under-reporting effects of women ...\n")
pd.urreport.f <- as.data.table(reshape2::melt(po$log_effect_age_under_report_f))
setnames(pd.urreport.f, 1:3, c('iterations', 'fm_mat_idx', 'value'))

# get the ur effects of women
tmp <- pd.urreport.f[,
                     list(sd = sd((value)),
                          value_m = quantile((value), p = pds.quantiles, na.rm = TRUE),
                          stat = pds.quantilelabels),
                     by = c('fm_mat_idx')]

ur.f.ui <- dcast.data.table(tmp, fm_mat_idx + sd~stat, value.var = 'value_m')

# save the quantiles of the under-reporting effects of women to file
cat("\nWrite posterior of under-reporting effects to file")


  write.csv(ur.f.ui, file =paste(hpc, "output/post/HIV/", model_name, "/", model_name,"_posterior_within_commu_ur_women_fix.csv", sep=""), row.names = F)


tmp2 <- subset(dc, part.sex == 'F', select = c(fm_mat_idx, age_age_idx)) #, gamma))
pd.urreport.f <- merge(pd.urreport.f, tmp2, by = c('fm_mat_idx'), allow.cartesian = TRUE)

pd.urreport.f <- pd.urreport.f[, list(iterations,age_age_idx,value)]

cat("\nExtract age and gender specific contact intensities ...\n")
pd.f <- as.data.table(reshape2::melt(po$log_m_fm))
setnames(pd.f, 1:3, c('iterations', 'age_age_idx', 'value'))

tmp2 <- subset(dc, part.sex == 'F', select = c(age_age_idx, pop, part.comm_id)) #, gamma))
pd.f <- merge(pd.f, tmp2, by = c('age_age_idx'), allow.cartesian = TRUE)
# merge the under-reporting effects and fm contact intensities with the ur effects
pd.f <- merge(pd.f, pd.urreport.f, by = c('age_age_idx', 'iterations'))
pd.f[, value_m := exp(value.x + value.y)]

pd.f[, value.wo.ur_m := exp(value.x)]
# pd.f[, total.wo.ur_m := exp(value.x + log(pop))]

# get the fm contact intensities with the ur effects
tmp <- pd.f[,
            list(sd = sd(value_m),
                 value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                 stat = pds.quantilelabels),
            by = c('age_age_idx', 'part.comm_id')
]
tmp <- dcast.data.table(tmp, age_age_idx + part.comm_id + sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

# get the fm contact intensities without the ur effects
tmp.wo.ur <- pd.f[,
                  list(sd = sd(value.wo.ur_m),
                       value.wo.ur_m = quantile(value.wo.ur_m, p = pds.quantiles, na.rm = TRUE),
                       stat = pds.quantilelabels),
                  by = c('age_age_idx', 'part.comm_id')
]
tmp.wo.ur <- dcast.data.table(tmp.wo.ur, age_age_idx + part.comm_id + sd~stat, value.var = 'value.wo.ur_m')
tmp.wo.ur$label <- 'without_under_reporting'

tmp.fm <- rbind(tmp.wo.ur, tmp)
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pd.m <- as.data.table(reshape2::melt(po$log_m_mf))
po = NULL
gc()
setnames(pd.m, 1:3, c('iterations', 'age_age_idx', 'value'))

tmp2 <- subset(dc, part.sex == 'M', select = c(age_age_idx, pop, part.comm_id)) #, gamma))
pd.m <- merge(pd.m, tmp2, by = c('age_age_idx'), allow.cartesian = TRUE)
pd.m[, value_m := exp(value)]


# get the mf contact intensities
tmp <- pd.m[,
            list(sd = sd(value_m),
                 value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                 stat = pds.quantilelabels),
            by = c('age_age_idx', 'part.comm_id')
]
tmp <- dcast.data.table(tmp, age_age_idx + part.comm_id + sd~stat, value.var = 'value_m')
tmp$label <- 'without_under_reporting'
tmp[, part.sex := 'M']
cat("\nDone for men contacting women\n")
tmp.fm[, part.sex := 'F']

pds <- rbind(tmp.fm, tmp, use.names = TRUE, fill = TRUE)
str(pds)

# get the empirical contact intensities
str(dc)
emp <- dc[, list(part.sex,age_age_idx,part.age,cont.age,capped_cntcts_in,part,q.term,part.comm_id,sexp1in)]
emp[, cntct_intensity_rho_scaled_empirical := capped_cntcts_in/part/q.term]
emp[, cntct_intensity_empirical := capped_cntcts_in/part]

# combine the empiricals and the estimates
pds <- merge(pds, emp, by = c('part.comm_id', 'age_age_idx', 'part.sex'), all = T)

tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_contact-intensities-within-comm_age_sex.rds", sep="")
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)

# Extract age and gender specific marginal contact intensities----
cat("\nExtract age and gender specific marginal contact intensities ...\n")
# women-specific contact intensities
pd.f <- merge(pd.f, unique(dc[part.sex == 'F', list(age_age_idx, part.age, cont.age)]), by = c('age_age_idx'), all.x = T)
pd.f.gender <- pd.f[, list(value_m = sum(value.wo.ur_m, na.rm = T)),
                    by = c('part.age', 'iterations', 'part.comm_id')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age', 'part.comm_id')
]
tmp.wo.ur <- dcast.data.table(tmp, part.age + part.comm_id + sd~stat, value.var = 'value_m')
tmp.wo.ur$label <- 'without_under_reporting'

# with under-reporting effect
pd.f.gender <- pd.f[, list(value_m = sum(value_m, na.rm = T)),
                    by = c('part.age', 'iterations', 'part.comm_id')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age', 'part.comm_id')
]
tmp <- dcast.data.table(tmp, part.age + part.comm_id + sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

tmp <- rbind(tmp.wo.ur, tmp, use.names = T, fill = T)
tmp[, part.sex := 'F']
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pd.m <- merge(pd.m, unique(dc[part.sex == 'M', list(age_age_idx, part.age, cont.age)]), by = c('age_age_idx'), all.x = T)
pd.m.gender <- pd.m[, list(value_m = sum(value_m, na.rm = T)),
                    by = c('part.age', 'iterations', 'part.comm_id')]
tmp.m <- pd.m.gender[,
                     list(sd = sd(value_m),
                          value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                          stat = pds.quantilelabels),
                     by = c('part.age', 'part.comm_id')
]
tmp.m <- dcast.data.table(tmp.m, part.age + part.comm_id + sd~stat, value.var = 'value_m')
tmp.m$label <- 'without_under_reporting'
tmp.m[, part.sex := 'M']

pds <- rbind(tmp, tmp.m, use.names = T, fill = T)

tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_marginal-contact-intensities-with-comm_age_sex.rds", sep="")
cat("\nWrite marginal contact intensities to file ", tmp)
saveRDS(pds, file = tmp)


# Location level ----
# Extract gender specific contact intensities----
cat("\nExtract gender specific contact intensities ...\n")
# women-specific contact intensities
pop.f <- unique(dc[cont.sex == 'F', list(cont.age,cont.sex,part.comm_id,pop)])
setnames(pop.f, c('cont.age', 'cont.sex', 'pop'), c('part.age', 'part.sex', 'pop.part'))
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value.wo.ur_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.comm_id')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.comm_id')
]
pop.f.all <- pop.f[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.comm_id')]
tmp <- merge(tmp, pop.f.all, by = c('part.comm_id'), all = T)
tmp[, value_m := value_m/pop.all]
tmp.wo.ur <- dcast.data.table(tmp, part.comm_id + sd~stat, value.var = 'value_m')
tmp.wo.ur$label <- 'without_under_reporting'

# with under-reporting effect
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.comm_id')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.comm_id')
]
tmp <- merge(tmp, pop.f.all, by = c('part.comm_id'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.comm_id + sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

tmp.f <- rbind(tmp.wo.ur, tmp, use.names = T, fill = T)
tmp.f[, part.sex := 'F']
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pop.m <- unique(dc[cont.sex == 'M', list(cont.age,cont.sex,part.comm_id,pop)])
setnames(pop.m, c('cont.age', 'cont.sex', 'pop'), c('part.age', 'part.sex', 'pop.part'))
pd.m.gender <- merge(pd.m, pop.m, by = c('part.age', 'part.comm_id'), all = T)
pd.m.gender <- pd.m.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.comm_id')]
tmp <- pd.m.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.comm_id')
]
pop.m.all <- pop.m[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.comm_id')]
tmp <- merge(tmp, pop.m.all, by = c('part.comm_id'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.comm_id + sd~stat, value.var = 'value_m')

tmp$label <- 'without_under_reporting'
tmp[, part.sex := 'M']

pds <- rbind(tmp.f, tmp, use.names = T, fill = T)

cat("\nDone for men contacting women\n")
tmp <-  paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_marginal-contact-intensities_within-comm_sex.rds", sep="")
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)

# contact intensities by age-age and gender
# women-specific contact intensities
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value.wo.ur_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age','cont.age')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age','cont.age')
]
pop.f.all <- pop.f[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.age')]
tmp <- merge(tmp, pop.f.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp.wo.ur <- dcast.data.table(tmp, part.age + cont.age + sd~stat, value.var = 'value_m')
tmp.wo.ur$label <- 'without_under_reporting'

# with under-reporting effect
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age','cont.age')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age','cont.age')
]
tmp <- merge(tmp, pop.f.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.age + cont.age + sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

tmp.f <- rbind(tmp.wo.ur, tmp, use.names = T, fill = T)
tmp.f[, part.sex := 'F']
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pd.m.gender <- merge(pd.m, pop.m, by = c('part.age', 'part.comm_id'), all = T)
pd.m.gender <- pd.m.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age','cont.age')]
tmp <- pd.m.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age','cont.age')
]
pop.m.all <- pop.m[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.age')]
tmp <- merge(tmp, pop.m.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.age + cont.age + sd~stat, value.var = 'value_m')

tmp$label <- 'without_under_reporting'
tmp[, part.sex := 'M']

pds <- rbind(tmp.f, tmp, use.names = T, fill = T)
cat("\nDone for men contacting women\n")

# get the empirical contact intensities
str(dc)
emp <- dc[, list(part.sex,part.age,cont.age,capped_cntcts_in,part,q.term,part.comm_id,sexp1in)]
emp.part <- unique(emp[, list(part.comm_id, part, part.sex, part.age)])
emp.part <- emp.part[, list(part = sum(part, na.rm = T)), by = c('part.age', 'part.sex')]
emp <- emp[, list(capped_cntcts_in = sum(capped_cntcts_in, na.rm = T)), by = c('part.age', 'cont.age', 'part.sex')]
emp <- merge(emp, emp.part, by = c('part.age', 'part.sex'), all = T)
# empiricals are with under-reporting issue (but for men, there is no under-reporting issues)
emp[, cntct_intensity_empirical := capped_cntcts_in/part]

# combine the empiricals and the estimates
pds <- merge(pds, emp, by = c('part.age', 'cont.age', 'part.sex'), all = T)

tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_contact-intensities_within_age_age_sex.rds", sep="")
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)



# Extract age and gender specific marginal contact intensities----
cat("\nExtract age and gender specific marginal contact intensities ...\n")
# at the area level
# marginal contact intensities by age and gender
# women-specific contact intensities
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value.wo.ur_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age')
]
pop.f.all <- pop.f[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.age')]
tmp <- merge(tmp, pop.f.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp.wo.ur <- dcast.data.table(tmp, part.age + sd~stat, value.var = 'value_m')
tmp.wo.ur$label <- 'without_under_reporting'

# with under-reporting effect
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age')
]
tmp <- merge(tmp, pop.f.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.age + sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

tmp.f <- rbind(tmp.wo.ur, tmp, use.names = T, fill = T)
tmp.f[, part.sex := 'F']
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pd.m.gender <- merge(pd.m, pop.m, by = c('part.age', 'part.comm_id'), all = T)
pd.m.gender <- pd.m.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations', 'part.age')]
tmp <- pd.m.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels),
                   by = c('part.age')
]
pop.m.all <- pop.m[, list(pop.all = sum(pop.part, na.rm = T)), by = c('part.age')]
tmp <- merge(tmp, pop.m.all, by = c('part.age'), all = T)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, part.age + sd~stat, value.var = 'value_m')

tmp$label <- 'without_under_reporting'
tmp[, part.sex := 'M']

pds <- rbind(tmp.f, tmp, use.names = T, fill = T)
cat("\nDone for men contacting women\n")
tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_marginal-contact-intensities_within_age_sex.rds", sep="")
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)


# marginal contact intensities by gender
# women-specific contact intensities
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value.wo.ur_m * pop.part, na.rm = T)),
                           by = c('iterations')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels)]

pop.f.all <- pop.f[, list(pop.all = sum(pop.part, na.rm = T))]
tmp <- cbind(tmp, pop.f.all)
tmp[, value_m := value_m/pop.all]
tmp.wo.ur <- dcast.data.table(tmp, sd~stat, value.var = 'value_m')
tmp.wo.ur$label <- 'without_under_reporting'

# with under-reporting effect
pd.f.gender <- merge(pd.f, pop.f, by = c('part.age', 'part.comm_id'), all = T)
pd.f.gender <- pd.f.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations')]
tmp <- pd.f.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels)]
tmp <- cbind(tmp, pop.f.all)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, sd~stat, value.var = 'value_m')
tmp$label <- 'with_under_reporting'

tmp.f <- rbind(tmp.wo.ur, tmp, use.names = T, fill = T)
tmp.f[, part.sex := 'F']
cat("\nDone for women contacting men\n")

# and only now we repeat for mf
pd.m.gender <- merge(pd.m, pop.m, by = c('part.age', 'part.comm_id'), all = T)
pd.m.gender <- pd.m.gender[, list(value_m = sum(value_m * pop.part, na.rm = T)),
                           by = c('iterations')]
tmp <- pd.m.gender[,
                   list(sd = sd(value_m),
                        value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
                        stat = pds.quantilelabels)]
pop.m.all <- pop.m[, list(pop.all = sum(pop.part, na.rm = T))]
tmp <- cbind(tmp, pop.m.all)
tmp[, value_m := value_m/pop.all]
tmp <- dcast.data.table(tmp, sd~stat, value.var = 'value_m')

tmp$label <- 'without_under_reporting'
tmp[, part.sex := 'M']

pds <- rbind(tmp.f, tmp, use.name = T, fill = T)
cat("\nDone for men contacting women\n")
tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_marginal-contact-intensities_within_sex.rds", sep="")
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)

# Posterior predictive check on age and gender specific contact intensities----
cat("\nPosterior predictive check on age and gender specific contact intensities ...\n")

  # log_m_mf are the estimated log of contact intensities
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
pd <- NULL
gc()

# PPC check
y.fm <- as.data.table(reshape2::melt(po$y_pred_fm))
setnames(y.fm, 1:3, c('iterations', 'id', 'value'))

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
cat('\nDone for PPC for women reported data ...\n')

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
cat('\nDone for PPC for men reported data ...\n')

pds <- rbind(y.fm, y.mf, use.names = T, fill = T)

tmp <-paste(hpc, "output/post/HIV/", model_name, "/", model_name,"_posterior-predictions.rds", sep="")
cat("\nWrite posterior predictions to file ", tmp)
saveRDS(pds, file = tmp)


# PPC check diagnostic ----
# if the pop size is 0, then estimated contacts is 0
diagnostics = give_diagnostics(dc, model_fit, hpc, model_name)
set(pds, which(pds$pop == 0), c('CL', 'CU', 'IL', 'IU', 'M'), 0)

pds[, if_out := (capped_cntcts_in < CL | CU < capped_cntcts_in)]

tmp <- mean(pds[, as.integer(if_out)], na.rm = TRUE)
tmp <- 1 - round(tmp, 2)
diagnostics <- rbind(diagnostics, c('PPC_all', tmp))

tmp <- paste0(tmp * 100, '%')
cat("\nProportion of data points inside the 95% posterior prediction interval = ", tmp)

for (i in unique(dc$part.comm_id))
{
  tmp <- mean(pds[part.comm_id == i, as.integer(if_out)], na.rm = TRUE)
  tmp <- 1 - round(tmp, 2)
  diagnostics <- rbind(diagnostics, c(paste0('PPC_comm_',i) , tmp))
  tmp <- paste0(tmp * 100, '%')
  cat('\nCommunity', i , ' : ',tmp, 'data are inside the 95% P.I.\n')
}


# plot if the contacts are outside the UI
if (1)
{
  pds.save <- copy(pds)
  pds <- update_sex_facet_name(pds)
  # pds.tp[, ]
  pm <- ggplot(pds[part.age <= 65 & cont.age <= 65 ]) +
    geom_tile( aes(x = part.age, y = cont.age, fill = (if_out))) +
    coord_equal() +
       labs(fill = 'If the reports outside the estimated UI. in the previous year (Jun 2018-Nov 2020)') +
    facet_grid(part.comm_id ~ paste0(plt.part.sex, ' contacting ', plt.cont.sex)) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_colour_manual(values = (c(
      # "white",
      '#ffffcc',
      # "#ffffd9",
      # "#edf8b1",
      # "#c7e9b4",
      "#7fcdbb"
      # "#41b6c4",
      # "#1d91c0",
      # "#225ea8",
      # "#253494"
      # "#081d58"
    )),
    # values = scales::rescale( c(0,0.02, 0.04, 0.08, 0.12, 0.15, 0.3, 0.5, 1)),
    na.value = "grey90") +
    xlab("Age of contacting individuals") +
    ylab("Age of contacted individuals") +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 5, family = 'sans'),
          # axis.text.x = element_blank(),
          legend.text = element_text(size = 5, family = 'sans'),
          strip.text = element_text(size = 7, family = 'sans'),
          legend.title = element_text(size = 6, family = 'sans'),
          axis.title = element_text(size = 7, family = 'sans')
    )

  pm
}

ggsave(file = paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_PPC.png" , sep=""), pm,width = 10, height = 20)
cat('\nDone for the postprocessing!\n')







