# =============================================================================
# calibration_prior_distribution_6.r - Molecular Clock Calibration Prior Distributions
# =============================================================================
#
# DESCRIPTION:
#   This R script generates and visualizes prior probability distributions for fossil
#   calibration points in molecular clock dating analysis. It implements multiple
#   strategies for setting calibration priors using truncated Cauchy, skew-t, and
#   uniform distributions to accommodate uncertainty in fossil age estimates.
#
# PURPOSE:
#   Molecular clock dating requires specification of prior distributions on node ages
#   based on fossil calibrations. This script provides a systematic approach to
#   translate fossil age ranges into appropriate statistical distributions that
#   account for both hard and soft bounds on calibration ages.
#
# KEY FEATURES:
#   - Truncated Cauchy distributions for minimum age constraints
#   - Skew-t distributions for flexible calibration priors
#   - Uniform distributions for bounded age ranges
#   - Highest posterior density interval (HDI) calculations
#   - Automated plotting and table generation for multiple calibration strategies
#
# CALIBRATION STRATEGIES:
#   Strategy 1: Truncated Cauchy (conservative minimum age)
#   Strategy 2: Skew-t (flexible minimum age with upper bound)
#   Strategy 3: Truncated Cauchy (relaxed minimum age)
#   Strategy 4: Skew-t (most flexible calibration)
#
# DEPENDENCIES:
#   - sn: Skew-normal and skew-t distributions
#   - tfprobability: TensorFlow probability distributions
#   - purrr: Functional programming tools
#   - doParallel: Parallel computing support
#   - iterators: Iterator constructs
#
# OUTPUT:
#   - SVG plots of prior density distributions
#   - CSV tables with distribution parameters and statistics
#   - HDI95: 95% highest density intervals
#   - PTI: Prior tail intervals (expected quantiles)
#
# AUTHOR: [Original author]
# MODIFIED: [Date] - Added comprehensive documentation
# =============================================================================

## This code is for setting the skew-t and gamma distributions of fossil calibration
setwd('F:/BristolProjects/Nematoda/timing/20250929/PriorDistributionPDF')
# save.image('calibration_distribution.RData')

# Load required libraries for statistical distributions and parallel computing
library(sn)           # Skew-normal and skew-t distributions
library(tfprobability) # TensorFlow probability distributions
library(purrr)        # Functional programming tools
library(doParallel)   # Parallel computing support
library(iterators)    # Iterator constructs

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

# HPrDI95 - Calculate Highest Posterior Density Interval (95%)
# ------------------------------------------------------------------------------
# This function calculates the 95% highest posterior density interval for any
# quantile function. The HDI represents the narrowest interval containing 95%
# of the probability density.
#
# Parameters:
#   qfunc: Quantile function that takes a probability p and returns the quantile
#   st, ed: Start and end probabilities for interval search (default: 0.0001 to 0.9999)
#   by: Step size for probability grid (default: 0.0001)
#
# Returns:
#   Vector of length 2: [lower_bound, upper_bound] of the 95% HDI
#
# Algorithm:
#   1. Evaluate quantile function across probability grid
#   2. Calculate width of all possible 95% intervals
#   3. Return the narrowest interval (minimum width)
HPrDI95 <- function(qfunc, st = 0.0001, ed = 0.9999, by = 0.0001){
    # Define function to calculate interval width for each starting probability
    wqq <- function(p0){
        q0 = qfunc(p0)                    # Lower quantile
        q1 = qfunc(p0 + 0.95)            # Upper quantile (95% interval)
        x <- data.frame(q0, q1, w = q1 - q0)  # Width of interval
        return(x)
    }
    # Evaluate interval widths across probability grid and find minimum
    seq(st, ed - 0.95, by) |> map_dfr(\(x) wqq(x)) -> t
    return(unlist(t[which.min(t[,3]), -3]))  # Return narrowest interval bounds
}

# =============================================================================
# TRUNCATED CAUCHY DISTRIBUTION FUNCTIONS
# =============================================================================

# plot cauchy curve----
# The following two functions of truncated Cauchy distribution are from:
# Saralees Nadarajah & Samuel Kotz (2007) Programs in R for Computing Truncated
# Cauchy Distributions, Quality Technology & Quantitative Management, 4:3, 407-412

# dtcauchy - Density function for truncated Cauchy distribution
# ------------------------------------------------------------------------------
# Computes the probability density function for a truncated Cauchy distribution.
# The truncation allows for minimum age constraints in calibration priors.
#
# Parameters:
#   x: Vector of quantiles
#   tl: Truncation point (minimum age bound)
#   p: Location parameter offset (default: 0.1)
#   c: Scale parameter multiplier (default: 0.2)
#   pl: Lower tail probability (default: 0.025)
#
# Returns:
#   Vector of density values
dtcauchy <- function(x, tl = 1, p = 0.1, c = 0.2, pl = 0.025){
    loc = tl * (1 + p)    # Location parameter (mean of untruncated distribution)
    scl = tl * c          # Scale parameter (spread of distribution)

    # Calculate normalization constant for truncation
    r0 = pl / pcauchy(tl, loc, scl)

    d <- NULL
    for (t in x){
        if (t > tl) {
            # Above truncation point: standard Cauchy density normalized by survival probability
            d <- c(d, dcauchy(t, location = loc, scale = scl) /
                      (1 - pcauchy(tl, location = loc, scale = scl) - pl))
        } else {
            # Below truncation point: constant density for soft bounds
            d <- c(d, r0 * dcauchy(t, location = loc, scale = scl) /
                      (1 - pcauchy(tl, location = loc, scale = scl) - pl))
        }
    }
    return(d)
}

# dtcauchy2 - Alternative density function for truncated Cauchy distribution
# ------------------------------------------------------------------------------
# An alternative implementation of truncated Cauchy density with different
# parameterization, allowing for both upper and lower tail probabilities.
#
# Parameters:
#   x: Vector of quantiles
#   tl: Truncation point (minimum age bound)
#   p: Location parameter offset (default: 0.1)
#   c: Scale parameter multiplier (default: 0.2)
#   pl: Lower tail probability (default: 0.025)
#   pu: Upper tail probability (default: 0.975)
#
# Returns:
#   Vector of density values
dtcauchy2 <- function(x, tl = 1, p = 0.1, c = 0.2, pl = 0.025, pu = 0.975) {
    # Calculate distribution parameters
    A <- 0.5 + (1/pi) * atan(p/c)  # Normalization constant
    om <- (pu/pl) * (1 / (pi * A * c * (1 + (p/c)^2)))  # Omega parameter

    d <- NULL
    for (t in x){
        if (t > tl){
            # Above truncation: Cauchy-like density
            d <- c(d, pu * (1 / (A * pi * c * tl * (1 + ((t - tl * (1 + p)) / (c * tl)) ^ 2))))
        } else {
            # Below truncation: exponential-like density
            d <- c(d, ifelse(pl == 0, 0, pl * (om/tl) * (t/tl) ^ (om - 1)))
        }
    }
    return(d)
}

# qtcauchy - Quantile function for truncated Cauchy distribution
# ------------------------------------------------------------------------------
# Computes quantiles for the truncated Cauchy distribution.
#
# Parameters:
#   pt: Vector of probabilities
#   p: Location parameter offset (default: 0.1)
#   c: Scale parameter multiplier (default: 0.2)
#   tl: Truncation point (default: 1)
#   pl: Lower tail probability (default: 0.025)
#
# Returns:
#   Vector of quantile values
qtcauchy <- function(pt, p = 0.1, c = 0.2, tl = 1, pl = 0.025){
    loc = tl * (1 + p)    # Location parameter
    scl = tl * c          # Scale parameter

    # Calculate truncation adjustment
    L <- pcauchy(tl, location = loc, scale = scl) - pl

    # Transform probabilities to account for truncation
    qcauchy(L + pt * (1 - L), location = loc, scale = scl)
}

# pdf_tcauchy - Create truncated Cauchy distribution object
# ------------------------------------------------------------------------------
# Creates a comprehensive object representing a truncated Cauchy distribution
# for calibration priors, including density function, quantiles, and statistics.
#
# Parameters:
#   tl: Truncation point (minimum age bound)
#   p: Location parameter offset (default: 0.1)
#   c: Scale parameter multiplier (default: 0.2)
#   pl: Lower tail probability (default: 0, hard bound)
#
# Returns:
#   List containing distribution parameters, quantiles, HDI, and statistics
pdf_tcauchy <- function(tl, p, c, pl = 0){
    loc = tl * (1 + p)    # Location parameter
    scl = tl * c          # Scale parameter

    # Create density function
    func <- function(x){dtcauchy2(x, tl = tl, p = p, c = c, pl = pl)}

    # Calculate key quantiles
    q01 <- qtcauchy(0.01, p, c, tl, pl)   # 1st percentile
    q99 <- qtcauchy(0.99, p, c, tl, pl)   # 99th percentile
    E.ql <- qtcauchy(pl, p, c, tl, pl)    # Expected lower bound
    E.qu <- qtcauchy(0.975, p, c, tl, pl) # Expected upper bound (97.5th percentile)

    # Calculate 95% highest density interval
    hpdi <- HPrDI95(qfunc = function(pt){qtcauchy(pt, p, c, tl, pl)})

    # Evaluate density at HDI points to find maximum
    xxx <- seq(hpdi[1], hpdi[2], by = 0.0001)
    xxx |> map_dbl(\(xxx) func(xxx)) -> dens
    vmax <- max(dens)     # Maximum density
    pmax <- xxx[which.max(dens)]  # Position of maximum density

    # Create model string for output
    mod <- paste('L(', tl, ',', p, ',', c, ',', ifelse(pl == 0, '1e-300', pl), ')', sep = '')

    return(list(tl = tl, p = p, c = c, pl = pl, loc = loc, scl = scl,
                q01 = q01, q99 = q99, E.ql = E.ql, E.qu = E.qu, hpdi = hpdi,
                vmax = vmax, pmax = pmax, mod = mod))
}

# plot.tcauchy - Plot truncated Cauchy distribution
# ------------------------------------------------------------------------------
# Creates a visualization of the truncated Cauchy distribution with key statistics.
#
# Parameters:
#   pdf: Distribution object from pdf_tcauchy()
#   hmin, hmax: Plot axis limits
#   TitlePrefix: Prefix for plot title
#   last: Boolean indicating if this is the last plot in a multi-plot layout
#
# Features plotted:
#   - Probability density curve (red)
#   - Location parameter (blue dashed line)
#   - 1st and 99th percentiles (dark cyan lines)
#   - 95% HDI (deep pink shaded rectangle)
#   - Prior tail interval (orange shaded rectangle)
#   - Distribution parameters and statistics as text labels
plot.tcauchy <- function(pdf, hmin, hmax, TitlePrefix = '', last = F){
    # Create density function for plotting
    func <- function(x){dtcauchy2(x, tl = pdf$tl, p = pdf$p, c = pdf$c, pl = pdf$pl)}

    # Position text legend based on mode location
    hleg <- pdf$tl + 0.66 * (hmax - pdf$tl)

    # Set up plot margins
    par(mar = c(ifelse(last, 5, 2), 5, 4, 2) + 0.1)

    # Plot density curve
    curve(func, hmin, hmax, n = 1001, type = 'l', col = 'red', lwd = 4,
          xlim = c(hmin, hmax), xaxt = ifelse(last, 's', 'n'),
          xlab = ifelse(last, '100 Ma', ''), ylab = ifelse(last, 'Probability Density', ''),
          cex = 2, cex.axis = 2, cex.lab = 2, mar = c(ifelse(last, 5, 2), 10, 4, 2))

    # Add title
    title(main = paste(TitlePrefix, 'Truncated Cauchy Distribution: q.01 - q.99\n',
                       '(tl=', pdf$tl, ', p=', pdf$p, ', c=', pdf$c,
                       ', location=', pdf$loc, ', scale=', pdf$scl, ')', sep = ''),
          cex = 2)

    # Add reference lines and shaded regions
    abline(v = pdf$loc, col = 'blue', lty = 2, lwd = 3)  # Location parameter
    axis(3, at = c(pdf$q01, pdf$q99), labels = FALSE, lty = 1, lwd = 5, col = 'darkcyan')  # Quantiles
    rect(pdf$hpdi[1], pdf$vmax * 0.25, pdf$hpdi[2], pdf$vmax * 0.45,
         density = 20, angle = 45, col = 'deeppink')  # 95% HDI
    rect(pdf$E.ql, pdf$vmax * 0, pdf$E.qu, pdf$vmax * 0.2,
         density = 20, angle = 45, col = 'darkorange')  # Prior tail interval

    # Add text labels with statistics
    text(x = hleg, y = pdf$vmax * 0.90, adj = c(0.5, 0.5),
         labels = pdf$mod, col = 'gray25', cex = 2)
    text(x = hleg, y = pdf$vmax * 0.80, adj = c(0.5, 0.5),
         labels = paste('95% HDI = [', round(pdf$hpdi[1], 4), ', ', round(pdf$hpdi[2], 4), ']', sep = ''),
         col = 'deeppink', cex = 2)
    text(x = hleg, y = pdf$vmax * 0.70, adj = c(0.5, 0.5),
         labels = paste('PTI = [', round(pdf$E.ql, 4), ', ', round(pdf$E.qu, 4), ']', sep = ''),
         col = 'darkorange', cex = 2)
}

# =============================================================================
# SKEW-T DISTRIBUTION FUNCTIONS
# =============================================================================

# plot skew-t curve----
# pdf_skewt - Create skew-t distribution object for calibration priors
# ------------------------------------------------------------------------------
# Creates a comprehensive skew-t distribution object optimized for fossil
# calibration priors. The skew-t distribution provides flexible modeling of
# asymmetric calibration constraints with heavy tails.
#
# Parameters:
#   xi: Location parameter (central tendency)
#   nu: Degrees of freedom (controls tail heaviness, default: inferred)
#   ql: Lower quantile bound (minimum age)
#   qu: Upper quantile bound (maximum age, default: 10)
#   pl: Lower tail probability (default: 0.025)
#   pu: Upper tail probability (default: 0.025)
#   single.end: Boolean for one-sided vs two-sided constraints
#
# Returns:
#   List containing distribution parameters, quantiles, HDI, and statistics
#
# Algorithm:
#   Uses parallel optimization to find skew-t parameters that match desired
#   quantile constraints, allowing for flexible calibration prior specification.
pdf_skewt <- function(xi, nu, ql, qu = 10, pl = 0.025, pu = 0.025, single.end = F){
    # Nested optimization functions for parameter estimation
    QUx <- function(omx = 10){  # Find omega for upper bound only
        require(doParallel)
        qq <- function(x) {
            q <- data.frame(om = x,
                            qu = qst(omega = x, p = 1 - pu, xi = xi, alpha = alpha, nu = nu, tol = 1e-08, method = 0))
            return(q)
        }
        # Parallel computation setup
        cl <- makeCluster(ifelse(Sys.info()['sysname'] == 'Windows', detectCores(), detectCores() - 4))
        registerDoParallel(cl)

        # Iterative search for optimal omega
        omi <- 4
        om <- seq(0.001, omi, by = 0.0001)
        while (omi < omx){
            q <- foreach(o = om, .export = c('alpha', 'pu', 'xi', 'nu'), .combine = rbind, .packages = c('sn', 'purrr')) %dopar% {qq(o)}
            d <- abs(qu - q[,2]) # L1 norm minimization
            res <- as.numeric(q[which.min(d),])
            if (res[1] == omi){
                om <- seq(omi, omi + 2, by = 0.0001); omi = omi + 2
            } else {
                break
            }
        }
        stopCluster(cl)
        return(res)
    }
    QLx <- function(omx=10, alx=100){
        require(doParallel)
        qq <- function(x,y){
            q <- data.frame(om=x, al=y, 
                            ql = qst(p=pl, xi=xi, omega=x, alpha=y, nu=nu, tol=1e-08, method = 0))
            return(q)
        }
        cl <- makeCluster(ifelse(Sys.info()['sysname']=='Windows',detectCores(),detectCores()-4))
        registerDoParallel(cl)
        omi <- 4; alj <- 10
        om <- seq(0.001, omi, by=0.0001)
        al <- seq(0, alj, by=1)
        while (omi < omx && alj < alx){
            q <- foreach(o=om, .export = c('pl','pu','xi','nu'), .combine = rbind, .packages = c('sn','purrr')) %dopar% {map2_dfr(o, al, qq)}
            d <- abs(ql - q[,3]) # L1 norm minimum
            dm <- min(d)
            res <- as.numeric(q[d==dm,])
            if (res[1] == omi && res[2] == alj) {
                om <- seq(omi, omi+2, by=0.0001); omi = omi+2
                al <- seq(alj, alj+10, by=1); alj = alj+10
            }else if (res[1] == omi){
                om <- seq(omi, omi+2, by=0.0001); omi = omi+2
            }else if (res[2] == alj){
                al <- seq(alj, alj+10, by=1); alj = alj+10
            }else {
                break
            }
        }
        stopCluster(cl)
        return(res)
    }
    QLUx1 <- function(omx=10, alx=100){
        require(doParallel)
        qq <- function(x,y){
            q <- data.frame(om=x, al=y, 
                            ql = qst(p=pl, xi=xi, omega=x, alpha=y, nu=nu, tol=1e-08, method = 0), 
                            qu = qst(p=1-pu, xi=xi, omega=x, alpha=y, nu=nu, tol=1e-08, method = 0))
            return(q)
        }
        cl <- makeCluster(ifelse(Sys.info()['sysname']=='Windows',detectCores(),detectCores()-4))
        registerDoParallel(cl)
        omi <- 4; alj <- 10
        om <- seq(0.001, omi, by=0.0001)
        al <- seq(0, alj, by=1)
        while (omi < omx && alj < alx){
            qlu <- foreach(o=om, .export = c('pl','pu','xi','nu'), .combine = rbind, .packages = c('sn','purrr')) %dopar% {map2_dfr(o, al, qq)}
            ss <- (qlu[,3] - ql) ^ 2 / pl + (qlu[,4] - qu) ^ 2 / pu # weighted L2 norm minimum
            res <- as.numeric(qlu[which.min(ss),])
            if (res[1] == omi && res[2] == alj) {
                om <- seq(omi, omi+2, by=0.0001); omi = omi+2
                al <- seq(alj, alj+10, by=1); alj = alj+10
            }else if (res[1] == omi){
                om <- seq(omi, omi+2, by=0.0001); omi = omi+2
            }else if (res[2] == alj){
                al <- seq(alj, alj+10, by=1); alj = alj+10
            }else {
                break
            }
        }
        stopCluster(cl)
        return(res)
    }
    # Main parameter estimation logic for skew-t distribution
    if (pl == 0){
        # Hard lower bound case
        if (missing(nu)) nu <- ifelse(single.end, 1, 100)
        alpha <- 10000  # Extreme skewness for hard bound
        qux <- QUx(omx = 10)
        omega <- qux[1]
        E.ql <- ql  # Expected lower bound = specified lower bound
        E.qu <- qux[2]  # Expected upper bound from optimization
    } else {
        # Soft bounds case - optimize both bounds
        if (single.end) {
            # One-sided constraint (minimum age only)
            if (missing(nu)) nu <- 2
            qlx <- QLx(omx = 10, alx = 100)  # Optimize for lower bound
            omega <- qlx[1]
            alpha <- ceiling(qlx[2])
            E.ql <- qlx[3]
            E.qu <- qst(p = 1 - pu, xi = xi, omega = omega, alpha = alpha, nu = nu, tol = 1e-08, method = 0)
        } else if (pl >= 0.00001){
            # Two-sided constraint with moderate lower tail
            if (missing(nu)) nu <- 100
            qlux <- QLUx1(omx = 10, alx = 100)  # Joint optimization
            omega <- qlux[1]
            alpha <- ceiling(qlux[2])
            E.ql <- qlux[3]
            E.qu <- qlux[4]
        } else {
            # Two-sided constraint with very small lower tail
            if (missing(nu)) nu <- 5
            qlux <- QLUx2(omx = 10, alx = 4000)  # Alternative optimization
            omega <- qlux[1]
            alpha <- ceiling(qlux[2])
            E.ql <- qlux[3]
            E.qu <- qlux[4]
        }
    }

    # Create density function and calculate statistics
    func <- function(x){dst(x, xi = xi, omega = omega, alpha = alpha, nu = nu)}
    q01 <- qst(0.01, xi = xi, omega = omega, alpha = alpha, nu = nu, tol = 1e-08, method = 0)
    q99 <- qst(0.99, xi = xi, omega = omega, alpha = alpha, nu = nu, tol = 1e-08, method = 0)
    hpdi <- HPrDI95(qfunc = function(p) qst(p, xi = xi, omega = omega, alpha = alpha, nu = nu, tol = 1e-08, method = 0))

    # Find maximum density and its location
    xxx <- seq(hpdi[1], hpdi[2], by = 0.0001)
    xxx |> map_dbl(\(xxx) func(xxx)) -> dens
    vmax <- max(dens)
    pmax <- xxx[which.max(dens)]

    # Create model string for output
    mod <- paste('ST(', round(xi, 4), ',', round(omega, 4), ',', alpha, ',', nu, ')', sep = '')

    return(list(xi = xi, omega = omega, alpha = alpha, nu = nu, ql = ql, qu = qu, pl = pl, pu = pu,
                E.ql = E.ql, E.qu = E.qu, q01 = q01, q99 = q99, hpdi = hpdi, vmax = vmax, pmax = pmax, mod = mod))
}

# plot.skewt - Visualize skew-t distribution calibration prior
# ------------------------------------------------------------------------------
# Creates a comprehensive visualization of a skew-t distribution calibration
# prior, showing the probability density curve with key statistical annotations.
# Used for displaying flexible asymmetric calibration priors in molecular dating.
#
# Parameters:
#   pdf: Skew-t distribution object created by pdf_skewt()
#   hmin: Minimum x-axis value for plotting
#   hmax: Maximum x-axis value for plotting
#   TitlePrefix: Text prefix for plot title (e.g., "Strategy 1: ")
#   last: Boolean indicating if this is the last plot in a multi-panel figure
#
# Visual Elements:
#   - Red density curve showing probability distribution
#   - Blue dashed line: Location parameter (xi)
#   - Dark cyan lines: 1st and 99th percentiles
#   - Deep pink shaded rectangle: 95% highest density interval (HDI)
#   - Orange shaded rectangle: Prior tail interval (PTI)
#   - Text labels showing distribution parameters and key statistics
#
# Layout:
#   - Adjusts margins based on position in multi-panel figure
#   - X-axis labels only on bottom panel (when last=TRUE)
#   - Y-axis labels only on bottom panel (when last=TRUE)
#   - Legend positioned dynamically based on mode location
plot.skewt <- function(pdf, hmin, hmax, TitlePrefix='', last=F){
    # Create density function for the skew-t distribution
    func <- function(x){dst(x, xi = pdf$xi, omega = pdf$omega, alpha = pdf$alpha, nu = pdf$nu)}

    # Position legend dynamically based on mode location to avoid overlap
    if (pdf$pmax < hmin + 0.66 * diff(c(hmin, hmax))){
        hleg <- pdf$pmax + 0.66 * (hmax - pdf$pmax)  # Mode near left, place legend right
    }else{
        hleg <- pdf$pmax - 0.66 * (pdf$pmax - hmin)  # Mode near right, place legend left
    }

    # Set plot margins (larger bottom margin for last plot in series)
    par(mar = c(ifelse(last,5,2),5,4,2)+0.1)

    # Plot the density curve
    curve(func, hmin, hmax, n=1001, type = 'l', col = 'red', lwd = 4,
          xlim=c(hmin, hmax), xaxt = ifelse(last, 's', 'n'),
          xlab = ifelse(last,'100 Ma',''), ylab = ifelse(last, 'Probability Density',''),
          cex = 2, cex.axis = 2, cex.lab = 2, mar = c(ifelse(last,5,2),10,4,2))

    # Add title with distribution parameters
    title(main = paste(TitlePrefix, 'Skew-t Distribution: q.01 - q.99\n',
                       '(xi=', pdf$xi, ', omega=', round(pdf$omega,4),
                       ', alpha=', pdf$alpha, ', nu=', pdf$nu, ')', sep = ''),
          cex = 2)

    # Add reference lines and shaded regions
    abline(v=pdf$xi, col='blue', lty = 2, lwd = 3)  # Location parameter
    axis(3, at=c(pdf$q01,pdf$q99), labels=FALSE, lty = 1, lwd = 5, col = 'darkcyan')  # 1st/99th percentiles
    rect(pdf$hpdi[1], pdf$vmax * 0.25, pdf$hpdi[2], pdf$vmax * 0.45,
         density = 20, angle = 45, col = 'deeppink')  # 95% HDI
    rect(pdf$E.ql, pdf$vmax * 0, pdf$E.qu, pdf$vmax * 0.2,
         density = 20, angle = 45, col = 'darkorange')  # Prior tail interval

    # Add text labels with key statistics
    text(x = hleg, y = pdf$vmax * 0.90, adj = c(0.5,0.5),
         labels = pdf$mod, col = 'gray25', cex = 2)  # Distribution formula
    text(x = hleg, y = pdf$vmax * 0.80, adj = c(0.5,0.5),
         labels = paste('95% HDI = [', round(pdf$hpdi[1], 4), ', ', round(pdf$hpdi[2], 4), ']', sep = ''),
         col = 'deeppink', cex = 2)  # Highest density interval
    text(x = hleg, y = pdf$vmax * 0.70, adj = c(0.5,0.5),
         labels = paste('PTI = [', round(pdf$E.ql, 4), ', ', round(pdf$E.qu, 4), ']', sep = ''),
         col = 'darkorange', cex = 2)  # Prior tail interval
}

# =============================================================================
# UNIFORM DISTRIBUTION FUNCTIONS
# =============================================================================

# plot uniform curve----
# pdf_between - Create uniform distribution object for bounded age ranges
# ------------------------------------------------------------------------------
# Creates a uniform distribution object for calibration priors with hard or
# soft bounds on both minimum and maximum ages.
#
# Parameters:
#   lower: Lower bound (minimum age)
#   upper: Upper bound (maximum age)
#   p.lower: Lower tail probability (default: 0, hard bound)
#   p.upper: Upper tail probability (default: 0.025, soft upper bound)
#
# Returns:
#   List containing distribution parameters and metadata
pdf_between <- function(lower, upper, p.lower = 0, p.upper = 0.025){
    d = abs(upper - lower)  # Range width
    h.poss <- (1 - p.lower - p.upper) / d  # Height of uniform density

    # Bound type indicators
    l.bound <- ifelse(p.lower == 0, '(hard)', paste('(', p.lower, ')', sep = ''))
    u.bound <- ifelse(p.upper == 0, '(hard)', paste('(', p.upper, ')', sep = ''))

    # Model string for output
    mod <- paste('B(', lower, ',', upper, ',', ifelse(p.lower == 0, '1e-300', p.lower), ',', ifelse(p.upper == 0, '1e-300', p.upper), ')', sep = '')

    return(list(lower = lower, upper = upper, p.lower = p.lower, p.upper = p.upper,
                h.poss = h.poss, l.bound = l.bound, u.bound = u.bound, mod = mod))
}

# plot.between - Plot uniform distribution
# ------------------------------------------------------------------------------
# Creates a visualization of the uniform distribution for calibration priors.
#
# Parameters:
#   pdf: Distribution object from pdf_between()
#   hmin, hmax: Plot axis limits
#   TitlePrefix: Prefix for plot title
#   last: Boolean indicating if this is the last plot in a multi-plot layout
plot.between <- function(pdf, hmin, hmax, TitlePrefix = '', last = F){
    par(mar = c(ifelse(last, 5, 2), 5, 4, 2) + 0.1)

    # Create empty plot
    plot(NULL, xlim = c(hmin, hmax), ylim = c(0, pdf$h.poss * 1.2),
         xaxt = ifelse(last, 's', 'n'), xlab = ifelse(last, '100 Ma', ''),
         ylab = ifelse(last, 'Probability Density', ''), cex = 2, cex.axis = 2, cex.lab = 2)

    # Add title
    title(main = paste(TitlePrefix, 'Uniform Distribution:\n(', pdf$lower, pdf$l.bound, '-', pdf$upper, pdf$u.bound, ')'), cex = 2)

    # Plot uniform density as horizontal line
    segments(x0 = pdf$lower, y0 = pdf$h.poss, x1 = pdf$upper, col = 'red', lwd = 4)

    # Add vertical lines for bounds (solid for hard bounds, dashed for soft bounds)
    segments(x0 = pdf$lower, y0 = 0, y1 = pdf$h.poss, col = 'red', lwd = 4, lty = ifelse(pdf$p.lower == 0, 1, 4))
    segments(x0 = pdf$upper, y0 = pdf$h.poss, y1 = 0, col = 'red', lwd = 4, lty = ifelse(pdf$p.upper == 0, 1, 4))

    # Add model label
    hleg <- mean(c(pdf$lower, pdf$upper))
    text(x = hleg, y = pdf$h.poss * 0.90, adj = c(0.5, 0.5), labels = pdf$mod, col = 'gray25', cex = 2)
}

# =============================================================================
# MAIN CALIBRATION PRIOR GENERATION FUNCTION
# =============================================================================

# FinalPlot - Generate calibration prior distributions and plots
# ------------------------------------------------------------------------------
# Main function for creating molecular clock calibration priors. Implements
# four different strategies for translating fossil age constraints into
# statistical distributions suitable for Bayesian phylogenetic dating.
#
# Parameters:
#   lower.cali: Minimum age constraint (fossil age)
#   upper.cali: Maximum age constraint (NA for minimum-only constraints)
#   relax.lower: Relaxation probability for lower bound (default: NA)
#   output.table: Whether to write results to CSV (default: TRUE)
#   append.table: Whether to append to existing CSV (default: TRUE)
#   ord: Calibration number/order for file naming
#
# Strategies:
#   1. Truncated Cauchy (conservative minimum age prior)
#   2. Skew-t with upper bound (flexible minimum with maximum constraint)
#   3. Relaxed truncated Cauchy (less conservative minimum age)
#   4. Flexible skew-t (most accommodating prior)
#
# Outputs:
#   - SVG plot file with four strategy comparisons
#   - CSV table with distribution parameters and statistics
FinalPlot <- function(lower.cali = NA, upper.cali = NA, relax.lower = NA, output.table = T, append.table = T, ord = 1){
    pl0 <- ifelse(is.na(relax.lower), 0, relax.lower)  # Lower tail probability

    if (is.na(upper.cali)){
        # MINIMUM AGE CONSTRAINT ONLY (upper.cali = NA)
        # Create four different prior strategies for minimum age constraints

        # Strategy 1: Conservative truncated Cauchy
        pdf1 <- pdf_tcauchy(tl = lower.cali, p = 0.1, c = 0.2, pl = pl0)

        # Strategy 2: Skew-t with inferred upper bound from Strategy 1
        pdf2 <- pdf_skewt(xi = lower.cali, ql = lower.cali, qu = pdf1$E.qu, pl = pl0, pu = 0.025, single.end = T)

        # Strategy 3: More relaxed truncated Cauchy
        pdf3 <- pdf_tcauchy(tl = lower.cali, p = 0.5, c = 0.2, pl = pl0)

        # Strategy 4: Flexible skew-t with inferred upper bound
        pdf4 <- pdf_skewt(xi = lower.cali * 1.5, ql = lower.cali, qu = pdf3$E.qu, pl = ifelse(is.na(relax.lower), 0.001, max(relax.lower, 0.001)), pu = 0.025, single.end = T)

        # Calculate plot axis limits encompassing all distributions
        hrange <- range(c(lower.cali, pdf1$q01, pdf1$q99, pdf1$E.ql, pdf1$E.qu, pdf1$hpdi,
                          pdf2$q01, pdf2$q99, pdf2$E.ql, pdf2$E.qu, pdf2$hpdi,
                          pdf3$q01, pdf3$q99, pdf3$E.ql, pdf3$E.qu, pdf3$hpdi,
                          pdf4$q01, pdf4$q99, pdf4$E.ql, pdf4$E.qu, pdf4$hpdi))
        ehmin <- hrange[1] - 0.1 * diff(hrange)
        ehmax <- hrange[2] + 0.1 * diff(hrange)

        # Generate SVG plot with four panels
        svg(filename = paste('./', ord, '.PriorDensity.CST.1-4.svg', sep = ''), width = 12, height = 16)
        layout(matrix(1:4, 4, 1))  # 4-row, 1-column layout
        plot.tcauchy(pdf1, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 1: ', last = F)
        plot.skewt(pdf2, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 2: ', last = F)
        plot.tcauchy(pdf3, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 3: ', last = F)
        plot.skewt(pdf4, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 4: ', last = T)
        dev.off()

        # Generate CSV table with distribution parameters
        if (output.table) {
            TAB <- data.frame(CaliNum = rep(ord, 4),
                              Strategy = paste('CST', 1:4),
                              Minima = rep(lower.cali, 4),
                              Maxima = rep(upper.cali, 4),
                              Quantile.l = c(pdf1$ql, pdf2$ql, pdf3$ql, pdf4$ql),
                              Viol.Prob.l = c(pdf1$pl, pdf2$pl, pdf3$pl, pdf4$pl),
                              Quantile.u = c(pdf1$qu, pdf2$qu, pdf3$qu, pdf4$qu),
                              Viol.Prob.u = c(pdf1$pu, pdf2$pu, pdf3$pu, pdf4$pu),
                              Formula = c(pdf1$mod, pdf2$mod, pdf3$mod, pdf4$mod),
                              PTI.l = round(c(pdf1$E.ql, pdf2$E.ql, pdf3$E.ql, pdf4$E.ql), 4),
                              PTI.u = round(c(pdf1$E.qu, pdf2$E.qu, pdf3$E.qu, pdf4$E.qu), 4),
                              HDI95.l = round(c(pdf1$hpdi[1], pdf2$hpdi[1], pdf3$hpdi[1], pdf4$hpdi[1]), 4),
                              HDI95.u = round(c(pdf1$hpdi[2], pdf2$hpdi[2], pdf3$hpdi[2], pdf4$hpdi[2]), 4),
                              PeakPos = round(c(pdf1$pmax, pdf2$pmax, pdf3$pmax, pdf4$pmax), 4))

            # Write or append to CSV file
            if (append.table){
                write.table(TAB, file = paste('./AllPriorDensityDistribution.CST.1-4.csv', sep = ''),
                           sep = ",", na = "", col.names = ifelse(ord == 1, TRUE, FALSE),
                           row.names = FALSE, append = ifelse(ord == 1, FALSE, TRUE))
            } else {
                write.table(TAB, file = paste('./', ord, '.PriorDensityDistribution.CST.1-4.csv', sep = ''),
                           sep = ",", na = "", col.names = TRUE, row.names = FALSE, append = F)
            }
        }
    } else {
        # BOUNDED AGE CONSTRAINT (both minimum and maximum specified)
        # Create four different prior strategies for bounded age constraints

        # Strategy 1: Uniform distribution between bounds
        pdf1 <- pdf_between(lower = lower.cali, upper = upper.cali, p.lower = pl0, p.upper = 0.025)

        # Strategy 2: Skew-t centered on lower bound
        pdf2 <- pdf_skewt(xi = lower.cali, ql = lower.cali, qu = upper.cali, pl = pl0, pu = 0.025, single.end = F)

        # Strategy 3: Skew-t centered on midpoint
        pdf3 <- pdf_skewt(xi = mean(c(lower.cali, upper.cali)), ql = lower.cali, qu = upper.cali, pl = ifelse(is.na(relax.lower), 0.00001, max(relax.lower, 0.001)), pu = 0.025, single.end = F)

        # Strategy 4: Skew-t with relaxed lower bound
        pdf4 <- pdf_skewt(xi = mean(c(lower.cali, upper.cali)), ql = lower.cali, qu = upper.cali, pl = ifelse(is.na(relax.lower), 0.001, max(relax.lower, 0.001)), pu = 0.025, single.end = F)

        # Calculate plot axis limits
        hrange <- range(c(lower.cali, upper.cali, pdf2$q01, pdf2$q99, pdf2$E.ql, pdf2$E.qu, pdf2$hpdi,
                          pdf3$q01, pdf3$q99, pdf3$E.ql, pdf3$E.qu, pdf3$hpdi,
                          pdf4$q01, pdf4$q99, pdf4$E.ql, pdf4$E.qu, pdf4$hpdi))
        ehmin <- hrange[1] - 0.1 * diff(hrange)
        ehmax <- hrange[2] + 0.1 * diff(hrange)

        # Generate SVG plot
        svg(filename = paste('./', ord, '.PriorDensity.CST.1-4.svg', sep = ''), width = 12, height = 16)
        layout(matrix(1:4, 4, 1))
        plot.between(pdf1, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 1: ', last = F)
        plot.skewt(pdf2, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 2: ', last = F)
        plot.skewt(pdf3, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 3: ', last = F)
        plot.skewt(pdf4, hmin = ehmin, hmax = ehmax, TitlePrefix = 'Strategy 4: ', last = T)
        dev.off()

        # Generate CSV table
        if (output.table) {
            TAB <- data.frame(CaliNum = rep(ord, 4),
                              Strategy = paste('CST', 1:4),
                              Minima = rep(lower.cali, 4),
                              Maxima = rep(upper.cali, 4),
                              Quantile.l = c(lower.cali, pdf2$ql, pdf3$ql, pdf4$ql),
                              Viol.Prob.l = c(pdf1$p.lower, pdf2$pl, pdf3$pl, pdf4$pl),
                              Quantile.u = c(upper.cali, pdf2$qu, pdf3$qu, pdf4$qu),
                              Viol.Prob.u = c(pdf1$p.upper, pdf2$pu, pdf3$pu, pdf4$pu),
                              Formula = c(pdf1$mod, pdf2$mod, pdf3$mod, pdf4$mod),
                              PTI.l = round(c(pdf1$lower, pdf2$E.ql, pdf3$E.ql, pdf4$E.ql), 4),
                              PTI.u = round(c(pdf1$upper, pdf2$E.qu, pdf3$E.qu, pdf4$E.qu), 4),
                              HDI95.l = round(c(NA, pdf2$hpdi[1], pdf3$hpdi[1], pdf4$hpdi[1]), 4),
                              HDI95.u = round(c(NA, pdf2$hpdi[2], pdf3$hpdi[2], pdf4$hpdi[2]), 4),
                              PeakPos = round(c(NA, pdf2$pmax, pdf3$pmax, pdf4$pmax), 4))

            if (append.table){
                write.table(TAB, file = paste('./AllPriorDensityDistribution.CST.1-4.csv', sep = ''),
                           sep = ",", na = "", col.names = ifelse(ord == 1, TRUE, FALSE),
                           row.names = FALSE, append = ifelse(ord == 1, FALSE, TRUE))
            } else {
                write.table(TAB, file = paste('./', ord, '.PriorDensityDistribution.CST.1-4.csv', sep = ''),
                           sep = ",", na = "", col.names = TRUE, row.names = FALSE, append = F)
            }
        }
    }
}

# =============================================================================
# EXAMPLE CALIBRATION PRIOR GENERATION CALLS
# =============================================================================

# These examples demonstrate how to generate calibration priors for different
# fossil age constraints using the FinalPlot function. Each call creates:
# - An SVG plot comparing four calibration strategies
# - A CSV table with distribution parameters and statistics
#
# Notation:
# - B(min, max, p_lower, p_upper): Bounded constraint with min/max ages
# - L(min, p, c, p_lower): Lower bound constraint with truncated Cauchy parameters
#
# The 'ord' parameter numbers the calibrations sequentially for file naming.

# BOUNDED AGE CONSTRAINTS (minimum and maximum ages known)
# ------------------------------------------------------------

# 1. B(5.7351,6.09,1e-300,0.025) - Wide bounded range
FinalPlot(lower.cali=5.7351, upper.cali=6.09, relax.lower=NA, output.table=T, append.table = T, ord=1)

# 2. B(5.7351,6.3497,1e-300,0.025) - Even wider upper bound
FinalPlot(lower.cali=5.7351, upper.cali=6.3497, relax.lower=NA, output.table=T, append.table = T, ord=2)

# 3. B(5.7351,5.908,1e-300,0.025) - Narrow bounded range
FinalPlot(lower.cali=5.7351, upper.cali=5.908, relax.lower=NA, output.table=T, append.table = T, ord=3)

# 4. B(5.5025,6.09,1e-300,0.025) - Asymmetric bounded range
FinalPlot(lower.cali=5.5025, upper.cali=6.09, relax.lower=NA, output.table=T, append.table = T, ord=4)

# 6. B(5.2882,5.908,1e-300,0.025) - Another narrow range
FinalPlot(lower.cali=5.2882, upper.cali=5.908, relax.lower=NA, output.table=T, append.table = T, ord=6)

# 7. B(5.14,5.908,1e-300,0.025) - Narrow range with lower minimum
FinalPlot(lower.cali=5.14, upper.cali=5.908, relax.lower=NA, output.table=T, append.table = T, ord=7)

# 8. B(4.443,5.4209,1e-300,0.025) - Wide range with lower bounds
FinalPlot(lower.cali=4.443, upper.cali=5.4209, relax.lower=NA, output.table=T, append.table = T, ord=8)

# MINIMUM AGE CONSTRAINTS ONLY (no maximum age known)
# ---------------------------------------------------

# 5. L(2.591,0.1,0.2,1e-300) - Moderate minimum age
FinalPlot(lower.cali=2.591, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=5)

# 9. L(0.9817,0.1,0.2,1e-300) - Young minimum age
FinalPlot(lower.cali=0.9817, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=9)

# 10. L(4.054,0.1,0.2,1e-300) - Older minimum age
FinalPlot(lower.cali=4.054, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=10)

# 11. L(0.9415,0.1,0.2,1e-300) - Very young minimum age
FinalPlot(lower.cali=0.9415, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=11)

# 12. L(1.202,0.1,0.2,1e-300) - Young minimum age
FinalPlot(lower.cali=1.202, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=12)

# 13. L(1.202,0.1,0.2,0.025) - Relaxed lower bound probability
FinalPlot(lower.cali=1.202, upper.cali=NA, relax.lower=0.025, output.table=T, append.table = T, ord=13)

# 14. L(1.983,0.1,0.2,1e-300) - Moderate minimum age
FinalPlot(lower.cali=1.983, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=14)

# 15. L(0.347,0.1,0.2,1e-300) - Very young minimum age
FinalPlot(lower.cali=0.347, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=15)

# 16. L(2.355,0.1,0.2,1e-300) - Moderate minimum age
FinalPlot(lower.cali=2.355, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=16)

# 17. L(0.1762,0.1,0.2,1e-300) - Extremely young minimum age
FinalPlot(lower.cali=0.1762, upper.cali=NA, relax.lower=NA, output.table=T, append.table = T, ord=17)

# dev.off()  # Close any open graphics devices
