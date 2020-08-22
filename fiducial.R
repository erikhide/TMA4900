library(ggplot2)
library(ggthemes)
library(kdensity)
library(MASS)
source("util.R")

sample_from_fiducial <- function(r, numData) {
  x <- r / sqrt(1 - r^2)
  u1 <- rchisq(1, df=numData-1)
  u2 <- rchisq(1, df=numData-2)
  u3 <- mvrnorm(1, 0, 1)
  theta <- (x * u2 - u3) / u1
  return(theta / sqrt(1 + theta^2))
}

get_fiducial_distribution <- function(rho, numData, numSamples, knownMean=FALSE) {
  if (knownMean) {
    r <- get_empirical_with_known_means(rho, numData)
  }
  else {
    r <- get_empirical(rho, numData)
  }
  list <- rep(NA, numSamples)
  for (i in 1:numSamples) {
    list[i] <- sample_from_fiducial(r, numData)
  }
  return(list)
}

get_fiducial_approximate_density <- function(rho, numData, numSamples=1000, knownMean=FALSE) {
  samples <- get_fiducial_distribution(rho, numData, numSamples, knownMean)
  density_approximation <- get_density_approximation(samples)
}

plot_fiducial_distribution <- function(rho, numData, numSamples) {
  estimate_list <- data.frame(get_fiducial_distribution(rho, numData, numSamples))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}

plot_fiducial_density <- function(density_approximation) {
  rho <- seq(-1, 1, by=0.001)
  df <- data.frame(rho, density=density_approximation(rho))
  gg <- ggplot(df, aes(rho, density)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2) +
    xlab("Correlation") +
    ylab("Density")
  gg
}

plot_fiducials_together <- function(rho, numSamples) {
  df <- data.frame()
  rho_values <- seq(-1, 1, by=0.001)
  density <- rep(0, length(rho_values))
  n <- c(3, 10, 100)
  for (k in 1:2) {
    if (k == 1) {
      knownMean = FALSE
      empiricalType <- "Unknown mean"
    }
    if (k == 2) {
      knownMean = TRUE
      empiricalType <- "Known mean"
    }
    for (i in 1:3) {
      num_average <- 100
      for (j in 1:num_average) {
        if (j %% 10 == 0) {
          print(i)
          print(j)
        }
        density_approximation <- get_fiducial_approximate_density(rho, n[i], numSamples, knownMean)
        density <- density + density_approximation(rho_values)
      }
      density <- density / num_average
      df <- rbind(df, data.frame(rho=rho_values, density=density, n=n[i], empiricalType=empiricalType))
    }
  }
  gg <- ggplot(df, aes(rho, density)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    ylim(c(0, 3)) +
    facet_grid(empiricalType ~ n)
  return(gg)
}


get_MAP_estimate_fiducial <- function(rho, numData, numSamples) {
  samples <- get_fiducial_distribution(rho, numData, numSamples)
  density_approximation <- get_density_approximation(samples)
  result <- optimize(density_approximation, interval=c(-1, 1), maximum=TRUE)
  return(result$maximum)
}

get_MAP_estimate_fiducial_distribution <- function(rho, numData, numIter, numSamples=1000) {
  MAP_list <- rep(NA, numIter)
  for (i in 1:numIter) {
    MAP_list[i] <- get_MAP_estimate_fiducial(rho, numData, numSamples)
  }
  return(MAP_list)
}


get_Bayes_estimate_fiducial <- function(density_approximation) {
  result <- integrate(function(x) x * density_approximation(x), lower=-1, upper=1)
  return(result$value)
}

get_Bayes_estimate_fiducial_distribution <- function(rho, numData, numIter, numSamples=1000, knownMean=FALSE) {
  density_approximation <- get_fiducial_approximate_density(rho, numData, numSamples, knownMean)
  list <- rep(NA, numIter)
  for (i in 1:numIter) {
    list[i] <- get_Bayes_estimate_fiducial(density_approximation)
  }
  return(list)
}


minimise_Bayes_estimate_with_KL_divergence_fiducial <- function(density_approximation) {
  expected_value_func <- function(x) x * density_approximation(x)
  expected_value <- integrate(expected_value_func, lower=-1, upper=1)$value
  func_to_minimise <- function(x) 1/2 * log(1-x^2) + 1 / (1-x^2) - x / (1-x^2) * expected_value
  return(optimize(func_to_minimise, interval=c(-1, 1), maximum=FALSE)$minimum)
}

get_Bayes_estimate_fiducial_with_KL_divergence_distribution <- function(rho, numData, numIter, numSamples=1000, knownMean=FALSE) {
  list <- rep(NA, numIter)
  iter <- 0
  while (iter < numIter) {
    density_approximation <- get_fiducial_approximate_density(rho, numData, numSamples, knownMean)
    skip_to_next <- FALSE
    tryCatch(
      {
        list[iter+1] <- minimise_Bayes_estimate_with_KL_divergence_fiducial(density_approximation)
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    if(skip_to_next) {
      next
    }
    iter <- iter + 1
  }
  return(list)
}


minimise_Bayes_estimate_fiducial_with_Fisher_information <- function(density_approximation) {
  func_to_minimise <- function(rho_hat) {
    FIM_function_help <- function(x) abs(sqrt(2)*atanh(sqrt(2)*x/sqrt(1+x^2)) - asinh(x) - sqrt(2)*atanh(sqrt(2)*rho_hat/sqrt(1+rho_hat^2)) + asinh(rho_hat))
    FIM_function <- function(x) FIM_function_help(x) * density_approximation(x)
    return(integrate(FIM_function, lower=-1, upper=1)$value)
  }
  return(optimize(func_to_minimise, interval=c(-1, 1), maximum=FALSE)$minimum)
}

get_Bayes_estimate_fiducial_with_Fisher_information_distribution <- function(rho, numData, numIter, numSamples=1000, knownMean=FALSE) {
  list <- rep(NA, numIter)
  iter <- 0
  while (iter < numIter) {
    density_approximation <- get_fiducial_approximate_density(rho, numData, numSamples, knownMean)
    skip_to_next <- FALSE
    tryCatch(
      {
        list[iter+1] <- minimise_Bayes_estimate_fiducial_with_Fisher_information(density_approximation)
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    if(skip_to_next) {
      next
    }
    iter <- iter + 1
  }
  return(list)
}


rho <- 0.5
numData <- 3
numSamples <- 10000
#x <- get_fiducial_distribution(rho, numData, numSamples)
#density_approximation <- get_density_approximation(x)
#plot_fiducial_density(density_approximation)
#gg <- plot_fiducials_together(rho, numSamples)
