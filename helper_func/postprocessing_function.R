# Postprocessing functions for sexual network ----
# community level
# area level

# within the under-reporting effects of women
# without the under-reporting effects of women (from the model directly and the total number would be the same)
#

# basic ----
# add date time
add_time <- function()
{
  dt <- data.table(round = c('R015', 'R019'),
                   round.plt = c('R15', 'R19'),
                   date = c('Aug 2011-Aug 2013', 'Jun 2018-Nov 2020'))
  return(dt)
}

# update the sex name
update_sex_facet_name <- function(dt)
{
    dt[, plt.part.sex := ifelse(part.sex == 'F', 'Women', 'Men')]
    if( 'cont.sex' %in%  colnames(dt))
    {
        dt[, plt.cont.sex := ifelse(cont.sex == 'F', 'Women', 'Men')]
    }else{
        dt[, plt.cont.sex := ifelse(plt.part.sex == 'Men', 'Women', 'Men')]
    }
    return(dt)
}

# update the wording
update_model_facet_name <- function(dt)
{
  dt[grepl('all', model), model := 'Posterior median sexual partnerships intensitiy']
  dt[grepl('within', model), model := 'Posterior median sexual partnerships intensitiy\nin the same community']
  dt[grepl('external', model), model := 'Posterior median sexual partnerships intensitiy\noutside the community']
  return(dt)
}

update_label_facet_name <- function(dt)
{
  dt[grepl('without_', label), label := 'With under reporting adjustment']
  dt[grepl('with_', label), label := 'Without under reporting adjustment']
  return(dt)
}

update_location_facet_name <- function(dt)
{
  dt[grepl('fishing', part.comm), plt.comm := 'Fishing']
  dt[grepl('inland', part.comm), plt.comm := 'Inland']

  return(dt)
}

# posteriors ----
# load the samplers of the paramteres of interest
load_target_vars_diagnostic <- function()
{

    su.target.vars <- c('log_effect_community','log_effect_under_report_baseline_f', 'overdispersion','gp_rho_age1','gp_rho_age2','gp_sigma1','gp_sigma2','gp_rho_agef1', 'gp_rho_agef2', 'gp_sigma_f1','gp_sigma_f2','z', 'z_f', 'f_fm', 'log_effect_age_under_report_f')

    su.target.vars <- c('log_effect_community','overdispersion','lscale_1','lscale_2','gpscale','f_fm', 'z')
   
   
        su.target.vars <- c(su.target.vars,  'lscale_ur_1', 'lscale_ur_2', 'gpscale_ur', 'z_ur')
    

      su.target.vars <- c(su.target.vars, 'log_effect_age_under_report_f')


return(su.target.vars)
}

load_target_vars_posterior <- function()
{
  
    su.target.vars <- c('f_fm', 'log_effect_age_under_report_f', 'log_m_mf', 'log_m_fm')
    return(su.target.vars)
}


pre_process_dc <- function(dc){
  account_missing = TRUE
  community_name = "fishing" # 
  source(paste0(hpc,'helper_func/GP-functions.R') )
  source(paste0(hpc, "helper_func/2d_helper.R"))
  
 
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
  
  dc <- dc[part.comm == community_name & part.round == "R019"]
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
}

give_diagnostics <- function(dc, model_fit, hpc, model_name){
  
  variable_names = load_target_vars_diagnostic()
  su = as.data.table(posterior::summarise_draws(
    model_fit$draws(
      variables = variable_names,
      inc_warmup = FALSE)
  ))
  
  cat("\nThe minimal nESS is: ",su[,min(ess_bulk)], "\nThe maximal rhat is: ",su[,max(rhat)])
  write.csv(su, file =paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_convergence.csv", sep="") , row.names = TRUE)
  
  # Diagnostics----
  diagnostics <- posterior::summarise_draws(posterior::as_draws_df(
    model_fit$sampler_diagnostics()), ~quantile(.x, probs = c(0.01, 0.5, 0.99)
    ))
  diagnostics <- rbind(diagnostics, c('min_ess_bulk', su[,min(ess_bulk)],su[,min(ess_bulk)],su[,min(ess_bulk)]))
  diagnostics <- rbind(diagnostics, c('max_rhat', su[,max(rhat)],su[,max(rhat)],su[,max(rhat)]))
  write.csv(diagnostics, file = paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_chains_diagnostics.csv", sep="") , row.names = TRUE)
  
  tmp <- su[order(ess_bulk),][1:3,][,variable]
  bayesplot::color_scheme_set("mix-blue-red")
  p <- bayesplot::mcmc_trace(
    model_fit$draws(variables = tmp, inc_warmup = TRUE),
    pars = tmp,
    n_warmup = 300, # TODO please note
    facet_args = list(ncol = 1, strip.position = "top")
  ) +
    theme(legend.position = 'bottom')
  
  tmp <-paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_worst3traces.png", sep="")
  cat("\nWrite worst 3 traces to file ", tmp)
  ggsave(file = tmp, p, width = 10, height = 20)
  
  p <- model_fit$draws() %>%
    as_draws_df() %>%
    as_tibble() %>%
    select(starts_with(c( 'log_effect_community')))  %>%
    mcmc_intervals()
  
  tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_posterior-comm-inter.png", sep="")
  cat("\nSave posterior distribution to file ", tmp)
  ggsave(file = tmp, p, width = 4, height = 1 * (length(unique(dc$part.comm_id))  ), limitsize = FALSE)
  
  compute.time <- model_fit$metadata()$time
  tmp <- paste(hpc, "output/post/HIV/", model_name, "/", model_name, "_computing-time.csv", sep="")
  cat("\nWrite computing time to file ", tmp)
  write.csv(compute.time, file = tmp, row.names = F)
  return(diagnostics)
}

