library(data.table)
library(ggplot2)
library(units)
library(MBESS)

# constants
rho <- 41.742  # air density at STP, mol m-3

# df_si <- get_SI_prefix_table()
get_SI_prefix_table <- function() {
  df <- valid_udunits_prefixes()
  # change micro symbol from "µ" to "u"
  i <- match("micro", df$name)
  df$symbol[i] <- "u"
  return(df)
}

get_sigma_x <- function(x_max, n) {
  x <- seq(from = 0, to = x_max, length.out =  n)
  sigma_x <- sd(x) # set_units(sd(x), s)
  return(sigma_x)
}
# sigma_x <- get_sigma_x(x_max = 1800, n = 4)

# calculates the confidence interval
# the lmits are flux ± ci_flux
get_ci <- function(
  noise = 20, 
  SI_prefix = "nano", 
  height = 0.23, 
  t_max = 300, 
  n = 10) {
  
  sigma_x <- get_sigma_x(x_max = t_max, n)
  ci_b <- sqrt(noise^2 / ((n-1) * sigma_x^2)) * 1.96
  # convert ci_b from mol/mol/s to flux units mol/m2/s - check correct
  ci_flux <- ci_b * rho * height
  return(ci_flux)
}

get_ci_dist <- function(
  meanlog = 0, 
  sdlog = 1, 
  max_flux = 5, 
  noise = 20, 
  SI_prefix = "nano", 
  height = 0.23, 
  t_max = 300, 
  n = 10) {
  
  v_flux <- seq(0, max_flux, by = max_flux/100)
  v_prob <- dlnorm(v_flux, meanlog, sdlog)
  snr <- v_flux^2 / noise^2
  ci_flux <- get_ci(noise = noise, SI_prefix = SI_prefix, height = height, 
    t_max = t_max, n = n)
  df <- data.table(prob = v_prob, flux = v_flux, ci_flux = ci_flux, snr = snr, 
    noise = noise, t_max = t_max, height = height)
  return(df)
}

get_percent_detectable <- function(
  meanlog = 0, 
  sdlog = 1, 
  max_flux = 5, 
  noise = 20, 
  SI_prefix = "nano", 
  height = 0.23, 
  t_max = 300, 
  n = 10) {
  
  v_flux <- seq(0, max_flux, by = max_flux/100)
  v_prob <- dlnorm(v_flux, meanlog, sdlog)
  ci_flux <- get_ci(noise = noise, SI_prefix = SI_prefix, height = height, 
    t_max = t_max, n = n)
  percent_detectable <- 1 - plnorm(ci_flux, meanlog, sdlog)  * 100
  return(percent_detectable)
}
