# plot for paper
naturemed_reqs <- function()
{
  # call this before doing your plots
  regs <<- theme(axis.text = element_text(size = 5, family = 'sans'),
                 text = element_text(size = 7,family = 'sans'),
                 legend.text = element_text(size = 7, family = 'sans'))
}

ggarrange_nature <- function(
  ...,
  plotlist = NULL,
  ncol = NULL,
  nrow = NULL,
  labels = NULL,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  align = c("none", "h", "v", "hv"),
  widths = 1,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL,
  add_reqs=TRUE
){
  if (add_reqs)
  {
       reqs <- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))

  }

  plots <- c(list(...), plotlist)

  if (add_reqs)
  {    plots <- lapply(plots, function(p){p + reqs})


  }

  out <- ggarrange(plotlist = plots,
                   ncol = ncol,
                   nrow = nrow,
                   labels = labels,
                   label.x = label.x,
                   label.y = label.y,
                   hjust = hjust,
                   vjust = vjust,
                   font.label = list(size = 8, color = "black", face = "bold", family = 'sans'),
                   align = align,
                   widths = widths,
                   heights = heights,
                   legend = legend,
                   common.legend = common.legend,
                   legend.grob = legend.grob)
  return(out)
}

ggsave_nature <- function(filename, p, w = 18,h = 24, add_reqs=TRUE)
{
  # check size
  tmp <- sort(c(w,h))
  if (tmp[1] > 18 | tmp[2] > 24)
    warning('Plot is bigger than allowed for EDFs. Maximum size is 18cm x 24cm\n')
  if ( tmp[1] < 10)
    warning('w and h represent cm units, not inches. Are you sure you want to save such a small plot?\n')

  # apply changes: for the moment only works for simple plots, but not for ggarrange.
  # Let's see if anyone answers this:
  # https://stackoverflow.com/questions/74379207/simultaneously-applying-same-modification-to-all-ggarrange-subplots
  if (add_reqs)
  {
    reqs <- theme(axis.text = element_text(size = 5, family = 'sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))
    p <- p + reqs
  }

  # save
  ggsave(filename=filename, plot=p, width=w, height=h, units='cm', dpi=310)
}



# Function for plotting contact intensities patterns and marginal age contact intensities posterior distribution----

# plot four panels contact intensities patterns----


plot_contact_intensity_matrix = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    facet_grid( ~ paste0("Participants' sex: ",  part.sex )) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +

    xlab("Age of contacting individuals") +
    ylab("Age of contacted individuals") +
    theme(legend.position = 'bottom')
}
plot_contact_intensity_sd_matrix = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = sd)) +
    coord_equal() +
    # geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    facet_grid( ~ paste0("Participants' sex: ",  part.sex )) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    xlab("Age of contacting individuals") +
    ylab("Age of contacted individuals") +
    theme(legend.position = 'bottom')
}

plot_contact_intensity_matrix_empirical = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    facet_grid( ~ factor(label, levels = unique(d22$label)) + paste0("Participants' sex: ", part.sex )) +
    viridis::scale_fill_viridis(na.value = "white",  option = 'H') +

    # scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    xlab("Age of contacting individuals") +
    ylab("Age of contacted individuals") +
    theme(legend.position = 'bottom')
}

# plot four panels contact rate patterns----
plot_contact_rate_matrix = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    paste0("Participants' sex: ",  part.sex ) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

plot_contact_rate_sd_matrix = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = sd)) +
    coord_equal() +
    # geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    paste0("Participants' sex: ",  part.sex ) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    theme(legend.position = 'bottom')
}
plot_contact_rate_matrix_empirical = function(d22)
{

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    geom_abline(slope = 1, intercept = 0, colour = 'white', size = 1.5, linetype = 2) +

    facet_grid( ~ factor(label, levels = unique(d22$label)) + paste0("Participants' sex: ",  part.sex )) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    theme(legend.position = 'bottom')
}


# plot posterior age contact intensities distribution----
plot_age_contact_intensity_MF = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                           colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_intensity_FM = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_intensity_MF_empirical = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    geom_point(aes(y = cntct_intensity_empirical), col = "black", size = 2, shape = "*") +
    geom_point(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", size = 2, shape = "+") +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_intensity_FM_empirical = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    geom_point(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", size = 2, shape = "+") +
    geom_point(aes(y = cntct_intensity_empirical), col = "black", size = 2, shape = "*") +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

# plot posterior age contact rate distrbution ----
plot_age_contact_rate_MF = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_rate_FM = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_rate_MF_empirical = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    geom_point(aes(y = cntct_rate_empirical), col = "black", shape = "*") +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_rate_FM_empirical = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("Participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("Participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    geom_point(aes(y = cntct_rate_empirical), col = "black", shape = "*") +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}
