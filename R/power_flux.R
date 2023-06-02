library(data.table)
library(ggplot2)
library(units)
library(mgcv)
library(caTools)
library(fitdistrplus)
# library(boot)
# library(MBESS)

# constants
# define the conversion unit between g N and moles of N2O
install_unit("mol_n2o", "28 g", "mol wt of N in N2O")
secsPerDay <- 24*60*60
rho <- set_units(41.742, mol / m^3) # air density at STP, mol m-3

# Function that returns the mode of a vector
Mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# Function that returns Root Mean Squared Error
rmse <- function(error){
    sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error){
    mean(abs(error))
}

# Function fits gam and lognormal model to time series of fluxes
# and calculates cumulative flux by these and trapezoidal method.
# Returns model parameters and predictions of time series and cumulative fluxes
fit_gam_logn <- function(df, t, Fn2o){
#df <- subset(mdf, site == "Edinburgh Grass"); t = "secs_since_fert"; Fn2o = "Fn2o_nmol_m2_s"
#df = subset(ch_df, fertID == "EBS_2006-08-08"); t = "secs_since_fert"; Fn2o = "Fn2o"
   
  df <- drop_units(df)
  # get arithmetic mean of flux at each time point
  Fn2o_mean <- aggregate(df[[Fn2o]], by = list(df[[t]]), FUN = mean)
  # rename agg outputs for time series of ch means
  t_mean <- Fn2o_mean[["Group.1"]]
  Fn2o_mean <- Fn2o_mean[["x"]]

  # time point to accumulate over 
  ####### check units ######
  t_cum <- 25 * secsPerDay # 
#  t_cum <- set_units(25 * secsPerDay, secs) # 

  # 1. Calculate cumulative flux by trapezoidal method
  Fcum_trap <- trapz(t_mean, Fn2o_mean)
  Fcum_trap <- trapz(c(t_mean, t_cum), c(Fn2o_mean, 0))
    
  # 2. Calculate cumulative flux for dndc output by gam method
  ########### need to add a try catch here - sometimes fails
  x <- t_mean 
  result <- try(mod.gam <- gam(Fn2o_mean ~ s(x, bs = "cr")))
  # fitted model predictions
  adf <- data.frame(x = t_cum)
  #flux_gam <- predict.gam(mod.gam, df, se=FALSE, newdata.guaranteed=TRUE)
  if (class(result) == "try-error"){ # fit gam didn't work
    # use lm instead
    mod.lm <- lm(Fn2o_mean ~ x)
    Fpred_gam <- predict(mod.lm, se=FALSE)
    Fcum_gam  <- predict(mod.lm, adf, se=FALSE)
  } else {   # use fit dist parameters
    Fpred_gam <- predict.gam(mod.gam, se=FALSE)
    Fcum_gam  <- predict.gam(mod.gam, adf, newdata.guaranteed=TRUE, se=FALSE)
  }
  
  # 3. Calculate cumulative flux for dndc output by lognormal method
  result <- try(fit <- fitdist(t_mean, "lnorm", weights = as.integer(Fn2o_mean*1e6+1)))
  if (class(result) == "try-error"){ # fitdist didn't work
    # use stats & Hmisc functions
    delta <- weighted.mean(log(t), (Fn2o_mean+1))
    k     <-  sqrt(wtd.var(log(t), (Fn2o_mean+1)))
  } else {   # use fit dist parameters
    delta <- as.numeric(fit$estimate["meanlog"])
    k     <- as.numeric(fit$estimate["sdlog"])
  }
  Fpred_logn <- dlnorm(t_mean, delta, k)
  fit <- lm(Fn2o_mean ~ 0 + Fpred_logn)
  alpha <- as.numeric(coef(fit))
  # fitted model predictions
  Fpred_logn <- Fpred_logn * alpha
  # predicted cumulative flux
  Fcum_logn <- plnorm(t_cum, delta, k) * alpha

  return(list(t_mean=t_mean, Fn2o_mean=Fn2o_mean, Fcum_trap=Fcum_trap, Fpred_gam=Fpred_gam, Fcum_gam=Fcum_gam, 
    delta=delta, k=k, alpha=alpha, Fpred_logn=Fpred_logn, Fcum_logn=Fcum_logn))
}

# df_si <- get_SI_prefix_table()
get_SI_prefix_table <- function() {
  df <- valid_udunits_prefixes()
  # change micro symbol from "µ" to "u"
  i <- match("micro", df$name)
  df$symbol[i] <- "u"
  return(df)
}

# TODO: make a worked example
# check with lm()
# n_gas <- 4
# x <- c(10, 20, 30, 40)
# y <- c(1.8, 2.2, 2.8, 4.1)
# m <- lm(y ~ x)
# summary(m)
# sigma_t <- sd(x)
# sigma_chi <- 0.3279
# se_b <- sqrt(sigma_chi^2 / ((n_gas - 1) * sigma_t^2))
# IS CORRECT!!

# TODO: should replace with sum of squares?
# = var(x) * (n_gas-1)
# Var(B) = sigma^2 / SSX
# or does n_gas just cancel out

get_v_t <- function(x_max, n_gas) {
  v_t <- seq(from = x_max / n_gas, to = x_max, length.out = n_gas)
  return(v_t)
}

get_sigma_t <- function(v_t) {
  sigma_t <- sd(v_t)
  units(sigma_t) <- units(v_t)
  return(sigma_t)
}


# calculates the confidence interval
# the lmits are flux ± ci_flux
  # sigma_t <- sapply(n_gas, get_sigma_t, x_max = t_max)
get_ci_flux <- function(sigma_chi = set_units(25, nmol_n2o / mol),
                     height   = set_units(0.23, m),
                     sigma_t  = set_units(90, s),
                     n_gas = 10) {  # no. samples
  ci_b <- sqrt(sigma_chi^2 / ((n_gas - 1) * sigma_t^2)) * 1.96 
  # convert ci_b from mol/mol/s to flux units mol/m2/s - check correct
  ci_flux <- ci_b * rho * height
  return(ci_flux)
}

get_dt_ci <- function(sigma_chi = set_units(25, nmol_n2o / mol),
                     height   = set_units(0.23, m),
                     sigma_t  = set_units(90, s),
                     n_gas = 10,        # no. samples
                     meanlog = 0,
                     sdlog = 1,
                     max_flux = set_units(5, nmol_n2o / m^2 /s)) {
  zero_flux <- max_flux - max_flux # get zero in flux units
  v_flux <- seq(zero_flux, max_flux, by = max_flux / 100)
  v_prob <- drop_units(dlnorm(v_flux, meanlog, sdlog))
  ci_flux <- get_ci_flux(
    sigma_chi = sigma_chi, height = height,
    sigma_t = sigma_t, n_gas = n_gas
  )
  snr <- v_flux^2 / ci_flux^2
  dt <- data.table(
    prob = v_prob, flux = v_flux, ci_flux = ci_flux, snr = snr,
    sigma_chi = sigma_chi, sigma_t = sigma_t, height = height
  )
  return(dt)
}

get_percent_detectable <- function(sigma_chi = set_units(25, nmol_n2o / mol),
                                   height   = set_units(0.23, m),
                                   sigma_t  = set_units(90, s),
                                   n_gas = 10,        # no. samples
                                   meanlog = 0,
                                   sdlog = 1
                                   ) {
  ci_flux <- get_ci_flux(
    sigma_chi = sigma_chi, height = height,
    sigma_t = sigma_t, n_gas = n_gas)
  percent_detectable <- (1 - drop_units(plnorm(ci_flux, meanlog, sdlog))) * 100
  return(percent_detectable)
}

get_sigma_spatial <- function(n_samples = 10, n_sims = 1e5, location = 0,
                              scale = 1.5, dist_lognormal = TRUE) {
  n <- rep(n_samples, n_sims)
  # assume same mean for either distribution (cannot have CV with mu = 0)
  mu_true <- exp(location + scale^2 / 2)
  if (dist_lognormal) {
    df <- as.data.frame(sapply(n, function(x) {
      rlnorm(x, location, scale)
    }))
  } else { # assume normal
    df <- as.data.frame(sapply(n, function(x) {
      rnorm(x, mu_true, exp(scale))
    }))
  }
  x <- colMeans(df)
  error <- x - mu_true
  sigma <- sd(error)
  cv <- sigma / mu_true * 100
  return(c(
    mu = mean(x), sigma = sigma, bias = mean(error),
    cv = cv, n_samples = n_samples, dist_lognormal = dist_lognormal
  ))
}

get_omega_sd <- function(
                         sigma_chi = set_units(5, nmol_n2o / mol),
                         height   = set_units(0.23, m),
                         t_max = set_units(5*60, s), # flux measurement time length, s
                         n_gas = 10,        # no. gas samples per mmnt
                         N_appl = set_units(0.1489069 * 1e9, nmol_n2o / m^2), # enter as mol/m2 * 1e9 = nmol/m2
                         omega = 0.01,
                         n_days = 18,
                         n_mmnt_per_day = 4,
                         n_sims = 3,
                         d_max = set_units(21 * secsPerDay, s), # time length over which measurements are taken, ~ 25 days but in secs
                         delta = 11.8, 
                         k = 0.86,
                         sigma_s = 1.5,
                         plot_graph = TRUE) {

df_params <- data.frame(omega   = omega,
                         N_appl = N_appl,
                         sigma_chi = sigma_chi,
                         n_days = n_days,
                         n_mmnt_per_day = n_mmnt_per_day,
                         n_sims = n_sims,
                         d_max = d_max, # time length over which measurements are taken, ~ 25 days but in secs
                         delta = delta, 
                         k = k,
                         sigma_s = sigma_s,
                         height = height,
                         t_max = t_max, # flux measurement time length, s
                         n_gas = n_gas)
                                                 
  n_mmnts <- n_days * n_mmnt_per_day
  # vector of measurement times
  v_times <- rep(seq(set_units(1, s), d_max, length.out = n_days), times = n_mmnt_per_day)
  # vector of true instantaneous mean flux
  # / set_units(1, s) is shortcut to make units correct
  v_F_mean <- drop_units(dlnorm(v_times, delta, k)) * omega * N_appl / set_units(1, s)
  # plot(v_times, v_F_mean)
  # mean(v_F_mean)
  
  # true cumulative flux
  F_cum <- drop_units(plnorm(d_max, delta, k)) * omega * N_appl
  
  F_mean <- plnorm(d_max, delta, k) * omega * N_appl  / d_max

  # add error from chamber flux mmnt
  v_t <- get_v_t(t_max, n_gas)
  sigma_t <- get_sigma_t(v_t)
  ci_flux <- get_ci_flux(
    sigma_chi = sigma_chi,
    height = height,
    sigma_t = sigma_t,
    n_gas = n_gas
  )
  sd_flux <- ci_flux / 1.96

  # add spatial variation
  # calculate mu_log from arithmetic mean mu
  v_F_meanlog <- drop_units(log(v_F_mean)) - 0.5 * sigma_s^2
  m_F_sample <- sapply(rep(n_mmnts, n_sims), rlnorm, v_F_meanlog, sigma_s)
  
  # add error from chamber flux mmnt
  m_F_obs <- sapply(1:n_sims, function(x) {
    rnorm(n_mmnts, mean = m_F_sample[, x], sd = sd_flux)
  })
  
  # get arithmetic mean of simulated flux at each time point
  df_F_obs_mean <- aggregate(m_F_obs, by = list(v_times), FUN = mean)
  # F_cum_obs <- trapz(df_F_obs_mean$Group.1, df_F_obs_mean[, 2:3])
  v_F_cum_obs <- sapply(1 + 1:n_sims, function(x) {
    trapz(df_F_obs_mean$Group.1, df_F_obs_mean[, x])
  })

  v_F_cum_obs <- set_units(v_F_cum_obs, nmol_n2o / m^2)

  v_omega_obs <- drop_units(v_F_cum_obs / N_appl)
  
  v_F_cum_error <- v_F_cum_obs - F_cum
  v_omega_error <- v_omega_obs - omega

# quantile(v_omega_error, probs = 0.66666666) - quantile(v_omega_error, probs = 0.6)
# rmse(v_omega_error)
# quantile(v_omega_error, probs = 0.975) - quantile(v_omega_error, probs = 0.5)
# rmse(v_omega_error) * 1.96
# hist(v_omega_error)

  F_cum_rmsd <- rmse(v_F_cum_error)
  omega_rmsd <- rmse(v_omega_error)
  # F_cum_error_sd <- sd(v_F_cum_error)
  # omega_error_sd <- sd(v_omega_error)

  if (plot_graph) {
    df <- data.frame(
      time = v_times / secsPerDay, F_mean = v_F_mean,
      F_sample = m_F_sample[, 1], F_obs = m_F_obs[, 1]
    )
    cols <- c("True mean" = "black", "With spatial varn" = "blue", 
      "With spatial varn + mmnt noise" = "red")
    p <- ggplot(df, aes(time, F_mean))
    p <- p + scale_colour_manual(name="", values = cols)
    p <- p + geom_line(aes(colour = "True mean"))
    p <- p + geom_point(aes(y = F_sample, colour = "With spatial varn"))
    p <- p + geom_point(aes(y = F_obs, colour = "With spatial varn + mmnt noise"))
    p <- p + geom_line(data = df_F_obs_mean, aes(Group.1 / secsPerDay, V1, 
      colour = "With spatial varn + mmnt noise"))
    p <- p + scale_y_unit(unit = "nmol_n2o / m^2 / s", limits = c(-2.5, 15))
    print(p)
  }

  df <- data.frame(
      F_cum_true = F_cum, omega_true = omega,
      F_cum_obs_mean = mean(v_F_cum_obs), omega_obs_mean = mean(v_omega_obs),
      F_cum_bias    = mean(v_F_cum_error), omega_bias = mean(v_omega_error),
      # F_cum_error_sd    = F_cum_error_sd, omega_error_sd = omega_error_sd,
      F_cum_rmsd    = F_cum_rmsd, 
      omega_rmsd = omega_rmsd,
      F_cum_error_ci    = F_cum_rmsd*1.96, 
      F_cum_error_pc = F_cum_error_sd / F_cum * 100,
      # omega_error_ci = omega_rmsd*1.96,
      omega_ci_pc = omega_rmsd * 1.96 * 100,
      F_cum_bias_pc = mean(v_F_cum_error) / F_cum * 100,
      omega_bias_pc = mean(v_omega_error) / omega * 100,
      ci_flux = ci_flux,
      df_params)
  # return(list(df = data.frame(
  # return(df = cbind.data.frame(df, df_params) #,
  return(df) #,
    # v_omega_obs = v_omega_obs,
    # v_F_cum_obs = v_F_cum_obs)
}

fit_logn <- function(v_F_obs) {
  # calculate offset
  if (min(v_F_obs) <= 0) {
    offset <- min(v_F_obs) * -1
    # need to add small frac of offset to be non-zero
    offset <- offset + offset * 0.01
  } else {
    offset <- 0
  }
  v_F_obs <- v_F_obs + offset
  fit <- fitdist(v_F_obs, "lnorm")
  v_F_obs <- v_F_obs - offset

  meanlog <- as.numeric(fit$estimate["meanlog"])
  sdlog <- as.numeric(fit$estimate["sdlog"])
  # arithmetic mean of logn dist
  mu <- exp(meanlog + sdlog^2 / 2)
  stdev <- exp(sdlog)
  return(c(meanlog = meanlog, sdlog = sdlog, mu = mu, stdev = stdev))
}



##* WIP - this is a bit different
# this gives a bootstrapped CI for a given sample
# we want to show the change in CI with n and sdlog
get_boot_ci <- function(x) {

  # Arithmetic mean in natural units
  mu <- function(x) {
    m <- mean(log(x))
    s <- mean((log(x) - m)^2)
    return(exp(m + s / 2))
  }

  # Simulated data
  data0 <- exp(rnorm(100))
  x <- rlnorm(10)
  mean(data0)

  # Bootstrap
  boot_out <- boot(data = x, statistic = function(d, ind) {
    mu(d[ind])
  }, R = 10000)
  plot(density(boot_out$t))

  # 4 types of Bootstrap confidence intervals
  return(boot.ci(boot_out, conf = 0.95, type = "all"))
}
