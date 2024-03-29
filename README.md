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
Chamber-based measurements of N2O and CH4 fluxes have particular characteristics which make statistical inference problematic:  fluxes are typically not normally distributed, but are heavily right-skewed; the precision of analysers is low in relation to the measured signal; the signal to sigma_chi ratio is small, often below the detection limits as conventionally defined; measurements are variable in time and space, and show peaks and hot spots which are unpredictable. To estimate the effect of treatments on the cumulative flux requires us to quantify and propagate the uncertainty properly.
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
# sigma_chi <- set_units(0.3279, nmol_n2o/mol)
sigma_chi <- set_units(20, nmol_n2o/mol)
height <- set_units(0.23, m)
t_max <- set_units(5 * 60, s)
n_gas <- 10
v_t <- get_v_t(t_max, n_gas)
sigma_t <- get_sigma_t(v_t)

get_ci_flux
```

```
## function(sigma_chi = set_units(25, nmol_n2o / mol),
##                      height   = set_units(0.23, m),
##                      sigma_t  = set_units(90, s),
##                      n_gas = 10) {  # no. samples
##   ci_b <- sqrt(sigma_chi^2 / ((n_gas - 1) * sigma_t^2)) * 1.96 
##   # convert ci_b from mol/mol/s to flux units mol/m2/s - check correct
##   ci_flux <- ci_b * rho * height
##   return(ci_flux)
## }
```

```r
get_ci_flux(sigma_chi = sigma_chi,
            height = height,
            sigma_t = sigma_t,
            n_gas = n_gas)
```

```
## 1.381144 [nmol_n2o/m^2/s]
```

To illustrate, we can calculate the CI in flux measurements from a chamber of height 0.23 m and an analyser with sigma_chi of 20 nmol/mol (standard deviation) and 10 samples recorded over 5 minutes. 
 

```r
meanlog <- 0
sdlog <- 1.5
max_flux <- set_units(5, nmol_n2o / m^2 /s)
get_dt_ci(height = height, sigma_chi = sigma_chi, sigma_t = sigma_t, n_gas = n_gas,
            meanlog = 0, sdlog = 1, max_flux = max_flux)
```

```
##            prob                  flux                   ci_flux
##   1: 0.00000000 0.00 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##   2: 0.08977828 0.05 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##   3: 0.28159019 0.10 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##   4: 0.43983718 0.15 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##   5: 0.54626787 0.20 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##  ---                                                           
##  97: 0.02428655 4.80 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##  98: 0.02364735 4.85 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##  99: 0.02302884 4.90 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
## 100: 0.02243021 4.95 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
## 101: 0.02185071 5.00 [nmol_n2o/m^2/s] 1.381144 [nmol_n2o/m^2/s]
##                   snr         sigma_chi      sigma_t   height
##   1:  0.000000000 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##   2:  0.001310576 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##   3:  0.005242303 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##   4:  0.011795182 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##   5:  0.020969213 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##  ---                                                         
##  97: 12.078266816 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##  98: 12.331207950 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
##  99: 12.586770236 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
## 100: 12.844953674 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
## 101: 13.105758263 [1] 20 [nmol_n2o/mol] 90.82951 [s] 0.23 [m]
```

The CI is a constant value given these characteristics. However, we might also consider this in relative terms compared to the flux itself, rather like a signal-to-sigma_chi ratio (SNR).
Therefore, we express this in relation to a range of typical or possible fluxes, up to a maximum value `max_flux`, with their probability of occuring specified by a lognormal distribution, with the parameter `meanlog` and `sdlog`.


```r
dt <- get_dt_ci(height = height, sigma_chi = sigma_chi, sigma_t = sigma_t, n_gas = n_gas,
            meanlog = 0, sdlog = 1, max_flux = max_flux)
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
v_height <- set_units(seq(0.1, 0.9, by = 0.1), m)
l_dt <- lapply(v_height, function(i){get_dt_ci(height = i)})
dt_h <- rbindlist(l_dt)
p <- ggplot(dt_h, aes(flux, prob))
p <- p + geom_area(aes(fill = flux > ci_flux))
p <- p + geom_line()
p <- p + facet_wrap(~ height) + ggtitle("Chamber height in m")
p
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### Variation in CI with analyser sigma_chi

Secondly, we see how CI is influenced by analyser sigma_chi.


```r
v_noise <- set_units(seq(1, 27, by = 3), nmol_n2o / mol)
l_dt <- lapply(v_noise, function(i){get_dt_ci(sigma_chi = i)})
dt_noise <- rbindlist(l_dt)
p <- ggplot(dt_noise, aes(flux, prob))
p <- p + geom_area(aes(fill = flux > ci_flux))
p <- p + geom_line()
p <- p + facet_wrap(~ sigma_chi) + ggtitle("Analyser sigma_chi in nmol/mol")
p
```

![](README_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### Combine sources of variation in CI

We can plot the effect of both of these on CI as an image, showing how the percentage of fluxes which will be detectable (greater than the upper CI) varies with chamber height and analyser sigma_chi. 


```r
dt <- data.table(expand.grid(height = v_height, sigma_chi = v_noise))
dt[, ci_flux := get_ci_flux(height = height, sigma_chi = sigma_chi)]
dt[, percent_detectable := get_percent_detectable(height = height, sigma_chi = sigma_chi)]

p <- ggplot(dt, aes(height, sigma_chi, fill = percent_detectable, 
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

To examine the effect of the above uncertainties on cumulative fluxes at the field scale, we use output from the DNDC model to represent the true time course of fluxes following a fertiliser application. The figure below shows daily values, with an initial peak and declining to near zero by day 10. This is typical of output from DNDC and similar process-based models.

![Fig x. DNDC model output.](README_files/figure-html/unnamed-chunk-8-1.png)

Using this output as the time course of the true mean, we can simulate sets of chamber samples which might be realised, given specified spatial variability and sigma_chi in the measurements system. We do this using the model of Levy et al (2018):

$$ \mu_{t} = \mathtt{dlnorm}(t, \Delta, k) N_\mathrm{in} \Omega $$

which gives the (arithmetic) mean flux at time $t$ following fertilisation with a quantity of nitrogen $N_\mathrm{in}$ and emission factor $\Omega$, and a time course given by the lognormal function with location $\Delta$ and scale $k$. At any time $t$, the N$_2$O flux has a lognormal distribution in space:

$$ F \sim \ln\mathcal{N}(\mu_{\mathrm{log,}t}, \sigma_s^2)$$
where
$$ \mu_{\mathrm{log,}t} = \mathrm{log}(\mu_{t}) - 0.5 \sigma_s^2 $$

where $\sigma_s$ is the spatial standard deviation of the log-transformed flux. As examples, the figure below shows five realisations in turn. 

<video controls loop><source src="README_files/figure-html/unnamed-chunk-9.webm" /></video>


The time course of the true mean is shown in black, and is the same in each case. Blue points show the true flux at each of the chamber locations, assuming a spatial variability given by $\sigma_s$ = 0.6. The red points show simulated measurements, with measurement sigma_chi in the system specified by $\sigma$ = 10 nmol/mol  (the residual term in the regression of $d \chi$ versus $d t$ in Equation 2 above) The red line shows the daily mean of measurements; the measured cumulative flux is commonly taken as the trapezoidal integrgation of the area under this curve. This can be compared with the true cumulative flux, the area under the black curve. By iterating many times with different values representing spatial variability, measurement sigma_chi and time length, sampling intensity, and chamber height, we can examine the variation in the uncertainty in estimates of the cumulative flux.

### Spatial variability
Assume a range of spatial variability with no measurement sigma_chi.
The uncertainty in the resulting cumulative flux is expressed as the coefficient of variation (CV): the standard deviation in the simulated estimates divided by the true mean $\times$ 100. We also show the bias: the mean difference between the simulated estimates and the true mean.

![](README_files/figure-html/sigma_s-1.png)<!-- -->

### Measurement sigma_chi
Assume a range of measurement sigma_chi with no spatial variability, and a chamber height of 23 cm.
![](README_files/figure-html/sigma_chi-1.png)<!-- -->

### Chamber height
Assume a range of chamber heights with measurement sigma_chi of 10 nmol/mol.

Chamber height interacts linearly with measurement sigma_chi.
![](README_files/figure-html/height-1.png)<!-- -->

### Sampling intensity: number of measurement days
Assume a range of measurement days with other values at their defaults.
![](README_files/figure-html/n_days-1.png)<!-- -->

### Sampling intensity: number of samples per flux measurement
Assume a range of number of samples per flux measurement with other values at their defaults.
![](README_files/figure-html/n_gas-1.png)<!-- -->

### Sampling intensity: number of flux measurements per day
Assume a range of number of flux measurements per day (at each sampling interval) with other values at their defaults.
![](README_files/figure-html/n_mmnt_per_day-1.png)<!-- -->

### Sampling intensity: length of chamber closure for flux measurement
Assume a range of chamber closure times for flux measurement, with other values at their defaults.
![](README_files/figure-html/t_max-1.png)<!-- -->

### Magnitude of fluxes: amount of N applied
Assume a range of chamber closure times for flux measurement, with measurement sigma_chi of 20 nmol/mol.
![](README_files/figure-html/N_appl-1.png)<!-- -->

### Summary
Tabulate values for specific system set-up.

![](README_files/figure-html/summary-1.png)<!-- -->

|system           |          F_cum_error_cv|          F_cum_error_ci|
|:----------------|-----------------------:|-----------------------:|
|GC               | 52.26282 [m^2/nmol_n2o]| 102.4351 [m^2/nmol_n2o]|
|Skyline          | 52.26282 [m^2/nmol_n2o]| 102.4351 [m^2/nmol_n2o]|
|QCL-Fast chamber | 52.26282 [m^2/nmol_n2o]| 102.4351 [m^2/nmol_n2o]|




## Uncertainty in  Between-Treatment Differences in Emission Factors
The table above shows the typical uncertainty (95 % CI) in cumulative fluxes to be of the order of 102 % of the true value. If we want to detect differences in cumulative fluxes between two treatments, the variance is additive, so this becomes 145 %. In practical terms, this means that the difference between cumulative fluxes (and therefore emission factors) in two treatments would need to be greater than around 145% of the mean to be statistically detectable.




