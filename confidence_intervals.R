library(ggplot2)
library(kdensity)
library(MASS)
source("bivariate_normal.R")
source("distributions.R")
source("util.R")

get_credible_interval <- function(posteriorFunc, statistic=NULL, numData=NULL, fiducial=FALSE) {
  if (fiducial) {
    posterior <- posteriorFunc
  }
  else {
    posterior <- function(x) posteriorFunc(x, statistic[1], statistic[2], numData) / normalise_posterior(posteriorFunc, statistic, numData)
  }
  epsilon <- 0.001
  lower <- -1
  step <- 0.1
  right <- TRUE
  while (step > epsilon) {
    if (right) {
      lower <- lower + step
      prob_temp <- integrate(posterior, lower=-1, upper=lower)$value
      if (abs(prob_temp - 0.025) < epsilon) {
        break
      }
      if (prob_temp > 0.025) {
        right <- FALSE
        step <- step / 2
      }
    }
    else {
      if (abs(lower - step) - 1 < 0.0001) {
        step <- step/2
      }
      lower <- lower - step
      prob_temp <- integrate(posterior, lower=-1, upper=lower)$value
      if (abs(prob_temp - 0.025) < epsilon) {
        break
      }
      if (prob_temp < 0.025) {
        right <- TRUE
        step <- step / 2
      }
    }
  }
  upper <- 1
  step <- 0.1
  left <- TRUE
  while (step > epsilon) {
    if (left) {
      upper <- upper - step
      prob_temp <- integrate(posterior, lower=upper, upper=1)$value
      if (abs(prob_temp - 0.025) < epsilon) {
        break
      }
      if (prob_temp > 0.025) {
        left <- FALSE
        step <- step / 2
      }
    }
    else {
      if (abs(upper + step) - 1 < 0.0001) {
        step <- step / 2
      }
      upper <- upper + step
      prob_temp <- integrate(posterior, lower=upper, upper=1)$value
      if (abs(prob_temp - 0.025) < epsilon) {
        break
      }
      if (prob_temp < 0.025) {
        left <- TRUE
        step <- step / 2
      }
    }
  }
  return(c(lower, upper))
}

is_in_credible_interval <- function(rhoValue, interval) {
  if (rhoValue < interval[1] | rhoValue > interval[2]) {
    return(FALSE)
  }
  return(TRUE)
}

get_coverage_probability <- function(posteriorFunc, rho, numIter, numData) {
  counter <- 0
  iter <- 0
  while (iter < numIter) {
    statistic <- simulate_bivariate_normal(rho, numData)
    skip_to_next <- FALSE
    tryCatch(
      {
        interval <- get_credible_interval(posteriorFunc, statistic, numData)
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    if(skip_to_next) {
      next
    }
    counter <- counter + is_in_credible_interval(rho, interval)
    iter <- iter + 1
  }
  return(counter / numIter)
}

get_coverage_probability_fiducial <- function(rho, numIter, numData, numSamples) {
  counter <- 0
  iter <- 0
  while (iter < numIter) {
    if (iter %% 100 == 0) print(iter)
    samples <- get_fiducial_distribution(rho, numData, numSamples, knownMean=TRUE)
    density_approximation <- get_density_approximation(samples)
    skip_to_next <- FALSE
    tryCatch(
      {
        interval <- get_credible_interval(density_approximation, fiducial=TRUE)
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    if(skip_to_next) {
      next
    }
    counter <- counter + is_in_credible_interval(rho, interval)
    iter <- iter + 1
  }
  return(counter / numIter)
}

get_gg_object_coverage_probability <- function(posteriorFunc, rho, numIter, numData, numSamples=NULL, fiducial=FALSE) {
  numRhos <- length(rho)
  coverage_list <- rep(NA, numRhos)
  for (i in 1:numRhos) {
    print(rho[i])
    if (fiducial) {
      coverage_list[i] <- get_coverage_probability_fiducial(rho[i], numIter, numData, numSamples)
    }
    else {
      coverage_list[i] <- get_coverage_probability(posteriorFunc, rho[i], numIter, numData)
    }
    print(coverage_list[i])
  }
  df <- data.frame(rho, coverage_list)
  gg <- ggplot(df, aes(x=rho, coverage_list)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2) +
    geom_hline(yintercept=0.95, linetype="dashed", size=0.8) +
    xlab("Correlation") +
    ylab("Coverage probability") +
    ylim(0.6, 1)
  return(gg)
}

plot_coverage_probabilities_together <- function(gg1, gg2, gg3) {
  df1 <- gg1$data
  df2 <- gg2$data
  df3 <- gg3$data
  df1$num_data <- rep(3, 19)
  df2$num_data <- rep(10, 19)
  df3$num_data <- rep(100, 19)
  df <- rbind(df1, df2, df3)
  gg <- ggplot(data=df, aes(rho, coverage_list)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2) +
    geom_hline(yintercept=0.95, linetype="dashed", size=0.8) +
    labs(
      y="Coverage probability",
      x="Correlation"
    ) +
    #ylim(c(0.94, 0.96)) +
    facet_wrap(~ num_data)
  return(gg)
}

plot_coverage_probabilities_together_PC <- function(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9) {
  df1 <- gg1$data
  df2 <- gg2$data
  df3 <- gg3$data
  df4 <- gg4$data
  df5 <- gg5$data
  df6 <- gg6$data
  df7 <- gg7$data
  df8 <- gg8$data
  df9 <- gg9$data
  df1$num_data <- rep(3, 19)
  df2$num_data <- rep(10, 19)
  df3$num_data <- rep(100, 19)
  df4$num_data <- rep(3, 19)
  df5$num_data <- rep(10, 19)
  df6$num_data <- rep(100, 19)
  df7$num_data <- rep(3, 19)
  df8$num_data <- rep(10, 19)
  df9$num_data <- rep(100, 19)
  df1$lambda <- rep("0.0818", 19)
  df2$lambda <- rep("0.0818", 19)
  df3$lambda <- rep("0.0818", 19)
  df4$lambda <- rep("1.249", 19)
  df5$lambda <- rep("1.249", 19)
  df6$lambda <- rep("1.249", 19)
  df7$lambda <- rep("5.360", 19)
  df8$lambda <- rep("5.360", 19)
  df9$lambda <- rep("5.360", 19)
  df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  gg <- ggplot(data=df, aes(rho, coverage_list)) +
    theme_grey(base_size = 15) +
    geom_line(size=1) +
    geom_hline(yintercept=0.95, linetype="dashed", size=0.6) +
    labs(
      y="Coverage probability",
      x="Correlation"
    ) +
    facet_grid(lambda ~ num_data)
  return(gg)
}

plot_coverage_probabilities_fiducial <- function() {
  df <- as.data.frame(
    rbind(
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(0.972, 0.987, 0.995, 0.994, 0.997, 0.999, 0.996, 0.998, 0.998, 1.000, 0.999, 0.998, 0.999, 0.997, 0.994, 0.998, 0.996, 0.983, 0.967),
        n = 3,
        knownMean = "Unknown mean"
      ),
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(0.997, 0.988, 0.970, 0.938, 0.930, 0.875, 0.815, 0.805, 0.765, 0.755, 0.747, 0.788, 0.840, 0.862, 0.900, 0.928, 0.969, 0.980, 0.996),
        n = 10,
        knownMean = "Unknown mean"
      ),
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(1.000, 0.999, 0.994, 0.973, 0.932, 0.897, 0.800, 0.589, 0.403, 0.171, 0.374, 0.610, 0.780, 0.887, 0.944, 0.970, 0.990, 1.000, 0.999),
        n = 100,
        knownMean = "Unknown mean"
      ),
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(0.993, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.998, 0.995),
        n = 3,
        knownMean = "Known mean"
      ),
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(0.997, 0.996, 0.971, 0.948, 0.907, 0.884, 0.825, 0.794, 0.802, 0.775, 0.789, 0.816, 0.827, 0.891, 0.916, 0.958, 0.977, 0.988, 0.999),
        n = 10,
        knownMean = "Known mean"
      ),
      cbind(
        rho = seq(-0.9, 0.9, by=0.1),
        prob = c(1.000, 0.997, 0.996, 0.980, 0.940, 0.878, 0.773, 0.618, 0.359, 0.173, 0.345, 0.602, 0.779, 0.882, 0.961, 0.978, 0.988, 0.996, 0.999),
        n = 100,
        knownMean = "Known mean"
      )
    )
  )
  df$rho <- sapply(df$rho, as.numeric)
  df$prob <- sapply(df$prob, as.numeric)
  df$n <- sapply(df$n, as.numeric)
  gg <- ggplot(data=df, aes(rho, prob)) +
    theme_grey(base_size = 25) +
    geom_line(size=1.2, aes(group=1)) +
    geom_hline(yintercept=0.95, linetype="dashed", size=0.6) +
    labs(
      y="Coverage probability",
      x="Correlation"
    ) +
    facet_grid(knownMean ~ n)
  return(gg)
}



rho <- seq(-0.9, 0.9, by=0.1)
numIter <- 1000
numData <- 100
#prior <- "flat"
#lambda <- 5.360
#posteriorFunc <- function(x1, x2, x3, x4) posterior_PC(x1, x2, x3, x4, lambda=lambda)
#posteriorFunc <- posterior_flat
posteriorFunc <- NULL
numSamples <- 1000
fiducial <- TRUE

cov_prob <- rep(NA, length(rho))
for (i in (length(cov_prob)-7):length(cov_prob)) {
  print(rho[i])
  cov_prob[i] <- get_coverage_probability_fiducial(rho[i], numIter, numData, numSamples)
}

#gg <- get_gg_object_coverage_probability(posteriorFunc, rho, numIter, numData, numSamples, fiducial)
#gg
#filename <- sprintf("figures/cov_prob/%s_iter%d_data%d_lambda%.3f.RDS", prior, numIter, numData, lambda)
#filename <- sprintf("figures/cov_prob/%s_iter%d_data%d.RDS", prior, numIter, numData)
#filename <- sprintf("figures/cov_prob/fiducial_samples%d_iter%d_data%d.RDS", numSamples, numIter, numData)
#saveRDS(gg, filename)



# For making PC prior plot
#gg1 <- readRDS("figures/cov_prob/PC_iter10000_data3_lambda0.0818.RDS")
#gg2 <- readRDS("figures/cov_prob/PC_iter10000_data10_lambda0.0818.RDS")
#gg3 <- readRDS("figures/cov_prob/PC_iter10000_data100_lambda0.0818.RDS")
#gg4 <- readRDS("figures/cov_prob/PC_iter10000_data3_lambda1.249.RDS")
#gg5 <- readRDS("figures/cov_prob/PC_iter10000_data10_lambda1.249.RDS")
#gg6 <- readRDS("figures/cov_prob/PC_iter10000_data100_lambda1.249.RDS")
#gg7 <- readRDS("figures/cov_prob/PC_iter10000_data3_lambda5.360.RDS")
#gg8 <- readRDS("figures/cov_prob/PC_iter10000_data10_lambda5.360.RDS")
#gg9 <- readRDS("figures/cov_prob/PC_iter10000_data100_lambda5.360.RDS")
#plot_coverage_probabilities_together_PC(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9)
