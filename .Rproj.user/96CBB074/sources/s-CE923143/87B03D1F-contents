library(ggplot2)
library(kdensity)
library(MASS)

get_MAP_estimate <- function(posteriorFunc, statistic, numData) {
  result <- optimize(function(x) posteriorFunc(x, statistic[1], statistic[2], numData), interval=c(-1, 1), maximum=TRUE)
  return(result$maximum)
}

get_MAP_distribution <- function(posteriorFunc, rho, numData, numIter) {
  MAP_list <- rep(NA, numIter)
  for (i in 1:numIter) {
    statistic <- simulate_bivariate_normal(rho, numData)
    MAP_list[i] <- get_MAP_estimate(posteriorFunc, statistic, numData)
  }
  return(MAP_list)
}

plot_MAP_distribution <- function(posteriorFunc, rho, numData, numIter) {
  MAP_list <- data.frame(get_MAP_distribution(posteriorFunc, rho, numData, numIter))
  colnames(MAP_list) <- c("value")
  gg <- ggplot(MAP_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


get_Bayes_estimate <- function(posteriorFunc, statistic, numData) {
  result <- integrate(function(x) x * posteriorFunc(x, statistic[1], statistic[2], numData) / normalise_posterior(posteriorFunc, statistic, numData), lower=-1, upper=1)
  return(result$value)
}

get_Bayes_estimate_distribution <- function(posteriorFunc, rho, numData, numIter) {
  value_list <- rep(NA, numIter)
  iter <- 0
  while (iter < numIter) {
    skip_to_next <- FALSE
    statistic <- simulate_bivariate_normal(rho, numData)
    tryCatch(
      {
        value_list[iter+1] <- get_Bayes_estimate(posteriorFunc, statistic, numData)
      },
      error = function(e) {
        skip_to_next <<- TRUE
      }
    )
    if(skip_to_next) next
    if (value_list[iter+1] >= 1 | value_list[iter+1] <= -1) next
    iter <- iter + 1
  }
  return(value_list)
}

plot_Bayes_estimate_distribution <- function(posteriorFunc, rho, numData, numIter) {
  estimate_list <- data.frame(get_Bayes_estimate_distribution(posteriorFunc, rho, numData, numIter))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


minimise_Bayes_estimate_with_KL_divergence <- function(posteriorFunc, statistic, numData) {
  expected_value_func <- function(x) x * posteriorFunc(x, statistic[1], statistic[2], numData) / normalise_posterior(posteriorFunc, statistic, numData)
  expected_value <- integrate(expected_value_func, lower=-1, upper=1)$value
  func_to_minimise <- function(x) 1/2 * log(1-x^2) + 1 / (1-x^2) - x / (1-x^2) * expected_value
  return(optimize(func_to_minimise, interval=c(-1, 1), maximum=FALSE)$minimum)
}

get_Bayes_estimate_with_KL_divergence_distribution <- function(posteriorFunc, rho, numData, numIter) {
  list <- rep(NA, numIter)
  iter <- 0
  while (iter < numIter) {
    statistic <- simulate_bivariate_normal(rho, numData)
    skip_to_next <- FALSE
    tryCatch(
      {
        list[iter+1] <- minimise_Bayes_estimate_with_KL_divergence(posteriorFunc, statistic, numData)
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


minimise_Bayes_estimate_with_Fisher_information <- function(posteriorFunc, statistic, numData) {
  func_to_minimise <- function(rho_hat) {
    FIM_function_help <- function(x) abs(sqrt(2)*atanh(sqrt(2)*x/sqrt(1+x^2)) - asinh(x) - sqrt(2)*atanh(sqrt(2)*rho_hat/sqrt(1+rho_hat^2)) + asinh(rho_hat))
    FIM_function <- function(x) FIM_function_help(x) * posteriorFunc(x, statistic[1], statistic[2], numData)
    return(integrate(FIM_function, lower=-1, upper=1)$value)
  }
  return(optimize(func_to_minimise, interval=c(-1, 1), maximum=FALSE)$minimum)
}

minimise_Bayes_estimate_with_Fisher_information_squared <- function(posteriorFunc, statistic, numData) {
  expected_value_func_1 <- function(x) atanh(sqrt(2)*x / sqrt(1+x^2)) * posteriorFunc(x, statistic[1], statistic[2], numData) / normalise_posterior(posteriorFunc, statistic, numData)
  expected_value_func_2 <- function(x) asinh(x) * posteriorFunc(x, statistic[1], statistic[2], numData) / normalise_posterior(posteriorFunc, statistic, numData)
  expected_value_1 <- integrate(expected_value_func_1, lower=-1, upper=1)$value
  expected_value_2 <- integrate(expected_value_func_2, lower=-1, upper=1)$value
  tanh_func <- function(x) atanh(sqrt(2)*x / sqrt(1+x^2))
  sinh_func <- function(x) asinh(x)
  func_to_minimise <- function(x) 2*tanh_func(x)^2 - 2*sqrt(2)*sinh_func(x)*tanh_func(x) + sinh_func(x)^2 - (4*tanh_func(x) - 2*sqrt(2)*sinh_func(x)) * 
                                  expected_value_1 + (2*sqrt(2)*tanh_func(x) - 2*sinh_func(x)) * expected_value_2
  result <- optimize(func_to_minimise, interval=c(-1, 1), maximum=FALSE)
  return(result$minimum)
}

get_Bayes_estimate_with_Fisher_information_distribution <- function(posteriorFunc, rho, numData, numIter) {
  list <- rep(NA, numIter)
  iter <- 0
  while (iter < numIter) {
    statistic <- simulate_bivariate_normal(rho, numData)
    skip_to_next <- FALSE
    tryCatch(
      {
        list[iter+1] <- minimise_Bayes_estimate_with_Fisher_information(posteriorFunc, statistic, numData)
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



#posteriorFunc <- function(x1, x2, x3, x4) posterior_PC(x1, x2, x3, x4, lambda=5.360)
#posteriorFunc <- posterior_J
rho <- seq(0.1, 0.9, by=0.2)
numData <- 3
numIter <- 10000
numSamples <- 1000
knownMean <- TRUE

for (i in 1:length(rho)) {
  #list <- get_Bayes_estimate_distribution(posteriorFunc, rho[i], numData, numIter)
  #list <- get_Bayes_estimate_with_KL_divergence_distribution(posteriorFunc, rho[i], numData, numIter)
  #list <- get_Bayes_estimate_with_Fisher_information_distribution(posteriorFunc, rho[i], numData, numIter)
  #list <- get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples=numSamples, knownMean=knownMean)
  #list <- get_Bayes_estimate_fiducial_with_KL_divergence_distribution(rho[i], numData, numIter, numSamples=numSamples, knownMean=knownMean)
  #list <- get_Bayes_estimate_fiducial_with_Fisher_information_distribution(rho[i], numData, numIter, numSamples=numSamples, knownMean=knownMean)
  #list <- get_MLE_distribution(rho[i], numData, numIter)
  list <- get_empirical_with_known_means_distribution(rho[i], numData, numIter)
  
  density_approximation <- get_density_approximation(list)
  #plot_density_approximation(density_approximation)
  expected_value <- get_expected_value_of_density_approximation(density_approximation)
  variance <- get_variance_of_density_approximation(density_approximation)
  print(rho[i])
  print(expected_value)
  print(variance)
}


