---
title: 'Guidelines for statistical inference from chamber measurements of N$_2$O and CH$_4$ fluxes'
author: 
- Peter Levy and Nick Cowan
date: Centre for Ecology and Hydrology, Bush Estate, Penicuik, EH26 0QB, U.K.
output:
  html_document: 
    toc: no
    keep_md: yes
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!--- { rendering -->

<!--- } -->

<!--- { startup -->

<!--- } -->

## Abstract
Chamber-based measurements of N2O and CH4 fluxes have particular characteristics which make statistical inference problematic:  fluxes are typically not normally distributed, but are heavily right-skewed; the precision of analysers is low in relation to the measured signal; the signal to noise ratio is small, often below the detection limits as conventionally defined; measurements are variable in time and space, and show peaks and hot spots which are unpredictable. To estimate the effect of treatments on the cumulative flux requires us to quantify and propagate the uncertainty properly.
We present guidelines for how to do this, and analysis of the necessary sample size (power analysis). We consider the uncertainty in individual flux measurements, spatial means, cumulative fluxes at the field scale, and between-treatment differences in emission factors. We show how these depend on: chamber height; analyser precision; the number of chamber measurements; the length of time and number of gas samples per measurement; and the magnitude and skewness of the spatial and temporal distributions.  The number of samples required to detect differences is larger than is commonly assumed.

## Individual chamber flux measurements
The flux of a GHG from the soil surface within a chamber is calculated from the rate of change in the mixing ratio $d \chi / d t$ (in mol GHG/mol air/ s), adjusted for the density of air per unit surface area:

$$F = \frac{d \chi} {d t} \times \rho \frac{V}{A}$$

where $F$ is the surface flux in mol GHG/m^2^/s, $\rho$ is the density of air in mol/m^3^, and $V$ and $A$ are the volume and area of the chamber (which simplifies to the height $h = V/A$ for most common chamber shapes). The term $\rho h$ can be considered a constant for a given chamber under typical conditions.

The term $d \chi / d t$ is usually estimated as the slope $\beta$ of a linear regression between $\chi$ and time $t$ during the chamber closure (although various methods accounting for non-linearity are also used).
The uncertainty in this term is estimated by the standard algebra of linear regression, where the 95% confidence interval (CI) in the slope $\beta$ is given:

$$CI_{\beta}^{95} = \sqrt{ \frac{\sigma^2} {n-1 \sum (t_i - \bar{t})^2} } \times \mathbb{T}$$

where $\sigma^2$ is the residual variance, $n$ is the number of data points, $\sum (t_i - \bar{t})^2$ is the variance in the $x$ independent variable, time, and $\mathbb{T}$ represents the t statistic, with a value of 1.96 for a sample size greater than 20. We can define a function which calculates this 95% confidence interval in chamber flux measurements, given the characteristic of the chamber and analyser.


```r
get_ci_flux
```

```
## function(noise = 20,
##                         SI_prefix = "nano",
##                         height = 0.23,
##                         t_max = 300,
##                         n = 10) {
##   sigma_x <- get_sigma_x(x_max = t_max, n)
##   ci_b <- sqrt(noise^2 / ((n - 1) * sigma_x^2)) * 1.96
##   # convert ci_b from mol/mol/s to flux units mol/m2/s - check correct
##   ci_flux <- ci_b * rho * height
##   return(ci_flux)
## }
```

To illustrate, we can calculate the CI in flux measurements from a chamber of height 0.23 m and an analyser with noise of 20 nmol/mol (standard deviation) and 10 samples recorded over 5 minutes. 
 

```r
get_ci_dist(height = 0.23, meanlog = 0, sdlog = 1, max_flux = 5, 
  SI_prefix = "nano", noise = 20, t_max = 5*60, n = 10)
```

```
##            prob flux ci_flux        snr noise t_max height
##   1: 0.00000000 0.00 1.24303 0.00000000    20   300   0.23
##   2: 0.08977828 0.05 1.24303 0.00000625    20   300   0.23
##   3: 0.28159019 0.10 1.24303 0.00002500    20   300   0.23
##   4: 0.43983718 0.15 1.24303 0.00005625    20   300   0.23
##   5: 0.54626787 0.20 1.24303 0.00010000    20   300   0.23
##  ---                                                      
##  97: 0.02428655 4.80 1.24303 0.05760000    20   300   0.23
##  98: 0.02364735 4.85 1.24303 0.05880625    20   300   0.23
##  99: 0.02302884 4.90 1.24303 0.06002500    20   300   0.23
## 100: 0.02243021 4.95 1.24303 0.06125625    20   300   0.23
## 101: 0.02185071 5.00 1.24303 0.06250000    20   300   0.23
```

The CI is a constant value given these characteristics. However, we might also consider this in relative terms compared to the flux itself, rather like a signal-to-noise ratio (SNR).
Therefore, we express this in relation to a range of typical or possible fluxes, up to a maximum value `max_flux`, with their probability of occuring specified by a lognormal distribution, with the parameter `meanlog` and `sdlog`.


```r
dt <- get_ci_dist(height = 0.23, meanlog = 0, sdlog = 1, max_flux = 5, 
  SI_prefix = "nano", noise = 20, t_max = 5*60, n = 10)
p <- ggplot(dt, aes(flux, prob))
p <- p + geom_area(aes(fill = flux > ci_flux))
p <- p + geom_line()
p
```

![](README_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The figure above shows the specified distribution of flux values in relation to the calculated CI of 1.24 nmol/m^2^/s. This shows that 
31% of fluxes (shaded red) will be lower than the CI, which can be considered as a limit of detection for the chamber-analyser system. 

### Variation in CI with chamber height

Using this function, we can examine how CI varies with various properties of the chamber-analyser system. Firstly, we see how CI is influenced by chamber height.


```r
v_height <- seq(0.1, 0.9, by = 0.1)
l_dt <- lapply(v_height, function(i){get_ci_dist(height = i)})
dt_h <- rbindlist(l_dt)
p <- ggplot(dt_h, aes(flux, prob))
p <- p + geom_area(aes(fill = flux > ci_flux))
p <- p + geom_line()
p <- p + facet_wrap(~ height) + ggtitle("Chamber height in m")
p
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### Variation in CI with analyser noise

Secondly, we see how CI is influenced by analyser noise.


```r
v_noise <- seq(1, 27, by = 3)
l_dt <- lapply(v_noise, function(i){get_ci_dist(noise = i)})
dt_noise <- rbindlist(l_dt)
p <- ggplot(dt_noise, aes(flux, prob))
p <- p + geom_area(aes(fill = flux > ci_flux))
p <- p + geom_line()
p <- p + facet_wrap(~ noise) + ggtitle("Analyser noise in nmol/mol")
p
```

![](README_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### Combine sources of variation in CI

We can plot the effect of both of these on CI as an image, showing how the percentage of fluxes which will be detectable (greater than the upper CI) varies with chamber height and analyser noise. 


```r
dt <- data.table(expand.grid(height = v_height, noise = v_noise))
dt[, ci_flux := get_ci_flux(height = height, noise = noise)]
dt[, percent_detectable := get_percent_detectable(height = height, noise = noise)]

p <- ggplot(dt, aes(height, noise, fill = percent_detectable, 
  text = paste("CI:", round(ci_flux, 3)))) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum()
p
```

![](README_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
# ggplotly(p, tooltip = "text")
```

We can add further dimensions for measurement length and number of data points.

## Uncertainty in Spatial Means
The lognormal distribution of fluxes makes it more difficult to estimate the mean flux over an area at a given time. We can illustrate this by comparing the sampling error as a function of sample size for lognormal versus normally-distributed data. It is hard to compare like-with-like here, but we compare a
lognormal distribution with both mean $\mu_{log}$ and standard deviation $\sigma_{log}$ equal to 1 on the log scale against a normal distribution
with the same mean and standard deviation in original units, given by
$\mu = \mathrm{exp}(\mu_{log} + \sigma_{log}^2/2)$
  and $\sigma = \mathrm{exp}(\sigma_{log})$.
  
Given these two distributions, we sample from these with a range of sample sizes and examine how close the sample mean is to the true mean. By repeating the simulation 10000 times, we can quantify how the typical error changes with sample size in normal and lognormally-distributed data.
We express this as the coefficient of variation ($CV  = \sigma / \mu \times 100$).
The figure below shows that the CV is twice as large for lognormally-distributed data at the equivalent sample size.


```r
df_logn <- data.frame(t(sapply(c(3:20, 25, 30, 40), get_sigma_spatial, location = 1)))
df_norm <- data.frame(t(sapply(c(3:20, 25, 30, 40), get_sigma_spatial, location = 1, dist_lognormal = FALSE)))
df <- rbind(df_norm, df_logn)
df$lognormal <- as.logical(df$dist_lognormal)
p <- ggplot(df, aes(n_samples, cv, colour = lognormal))
p <- p + geom_line()
p
```

## Uncertainty in Cumulative Fluxes at the Field Scale
tbc



## Uncertainty in  Between-Treatment Differences in Emission Factors
tbc



