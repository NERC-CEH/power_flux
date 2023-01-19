library(data.table)
library(ggplot2)
library(units)
library(caTools)
library(fitdistrplus)
library(boot)
# library(MBESS)

# constants
rho <- 41.742 # air density at STP, mol m-3

# df_si <- get_SI_prefix_table()
get_SI_prefix_table <- function() {
  df <- valid_udunits_prefixes()
  # change micro symbol from "µ" to "u"
  i <- match("micro", df$name)
  df$symbol[i] <- "u"
  return(df)
}

get_sigma_x <- function(x_max, n) {
  x <- seq(from = 0, to = x_max, length.out = n)
  sigma_x <- sd(x) # set_units(sd(x), s)
  return(sigma_x)
}

# calculates the confidence interval
# the lmits are flux ± ci_flux
get_ci_flux <- function(noise = 20,
                        SI_prefix = "nano",
                        height = 0.23,
                        t_max = 300,
                        n = 10) {
  sigma_x <- get_sigma_x(x_max = t_max, n)
  ci_b <- sqrt(noise^2 / ((n - 1) * sigma_x^2)) * 1.96
  # convert ci_b from mol/mol/s to flux units mol/m2/s - check correct
  ci_flux <- ci_b * rho * height
  return(ci_flux)
}

get_ci_dist <- function(meanlog = 0,
                        sdlog = 1,
                        max_flux = 5,
                        noise = 20,
                        SI_prefix = "nano",
                        height = 0.23,
                        t_max = 300,
                        n = 10) {
  v_flux <- seq(0, max_flux, by = max_flux / 100)
  v_prob <- dlnorm(v_flux, meanlog, sdlog)
  snr <- v_flux^2 / noise^2
  ci_flux <- get_ci_flux(
    noise = noise, SI_prefix = SI_prefix, height = height,
    t_max = t_max, n = n
  )
  df <- data.table(
    prob = v_prob, flux = v_flux, ci_flux = ci_flux, snr = snr,
    noise = noise, t_max = t_max, height = height
  )
  return(df)
}

get_percent_detectable <- function(meanlog = 0,
                                   sdlog = 1,
                                   noise = 20,
                                   SI_prefix = "nano",
                                   height = 0.23,
                                   t_max = 300,
                                   n = 10) {
  ci_flux <- get_ci_flux(
    noise = noise, SI_prefix = SI_prefix, height = height,
    t_max = t_max, n = n
  )
  percent_detectable <- (1 - plnorm(ci_flux, meanlog, sdlog)) * 100
  return(percent_detectable)
}

get_sigma_spatial <- function(n_samples = 10, n_sims = 1e5, location = 1,
                              scale = 1, dist_lognormal = TRUE) {
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

get_omega_sd <- function(noise = 0, # 20,
                         n_times = 5,
                         n_samples_per_time = 3,
                         n_sims = 4,
                         time_max = 25,
                         delta = 1, k = 0.5,
                         sigma_s = 0, # 0.5,
                         omega = 0.02,
                         N_appl = 20014,
                         SI_prefix = "nano",
                         height = 0.23,
                         t_max = 300,
                         n_samples_per_mmnt = 10,
                         plot_graph = FALSE) {
  n <- n_times * n_samples_per_time
  # vector of measurement times
  v_times <- rep(seq(0, time_max, length.out = n_times), times = n_samples_per_time)
  # vector of true mean flux
  v_F_mean <- dlnorm(v_times, delta, k) * omega * N_appl

  # add error from chamber flux mmnt
  ci_flux <- get_ci_flux(
    noise = noise,
    SI_prefix = SI_prefix,
    height = height,
    t_max = t_max,
    n = n_samples_per_mmnt
  )
  sd_flux <- ci_flux / 1.96
  # true cumulative flux
  F_cum <- plnorm(time_max, delta, k) * omega * N_appl

  for (i in 1:n_sims) {
    # add spatial variation
    # v_F_sample <- rlnorm(n, log(v_F_mean), sigma_s)
    m_F_sample <- sapply(rep(n, n_sims), rlnorm, log(v_F_mean), sigma_s)
    # add error from chamber flux mmnt
    # v_F_obs <- rnorm(n, v_F_sample, sd_flux)
    m_F_obs <- sapply(rep(n, n_sims), rnorm, v_F_sample, sd_flux)
    # get arithmetic mean of simulated flux at each time point
    df_F_obs_mean <- aggregate(m_F_obs, by = list(v_times), FUN = mean)
    # F_cum_obs <- trapz(df_F_obs_mean$Group.1, df_F_obs_mean[, 2:3])
    v_F_cum_obs <- sapply(1 + 1:n_sims, function(x) {
      trapz(df_F_obs_mean$Group.1, df_F_obs_mean[, x])
    })
    # str(v_F_cum_obs)

    v_F_cum_error <- v_F_cum_obs - F_cum
    v_omega_obs <- v_F_cum_obs / N_appl
    v_omega_error <- v_omega_obs - omega
  }

  if (plot_graph) {
    df <- data.frame(
      time = v_times, F_mean = v_F_mean,
      F_sample = v_F_sample, F_obs = m_F_obs[, 1]
    )
    p <- ggplot(df, aes(time, F_mean))
    p <- p + geom_line()
    p <- p + geom_point(aes(y = F_sample), colour = "blue")
    p <- p + geom_point(aes(y = F_obs), colour = "red")
    p <- p + geom_point(data = df_F_obs_mean, aes(Group.1, x), colour = "green")
    p
  }
  return(data.frame(
    F_cum_obs = v_F_cum_obs, omega_obs = v_omega_obs,
    F_cum_true = F_cum, omega_true = omega,
    F_cum_error = v_F_cum_error, omega_error = v_omega_error
  ))
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
