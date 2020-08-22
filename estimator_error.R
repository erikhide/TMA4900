library(kdensity)
source("bivariate_normal.R")

get_KL_divergence <- function(rho_1, rho_2) {
  return(-1/2 * log((1-rho_1^2)/(1-rho_2^2)) + (1-rho_1*rho_2)/(1-rho_2^2) - 1)
}

get_Fisher_information_metric <- function(rho_1, rho_2) {
  return(abs(sqrt(2)*atanh(sqrt(2)*rho_1/sqrt(1+rho_1^2)) - asinh(rho_1) - sqrt(2)*atanh(sqrt(2)*rho_2/sqrt(1+rho_2^2)) + asinh(rho_2)))
}

calculate_estimates <- function(estimateFunc, rho, numData, numIter) {
  n_rho <- length(rho)
  list <- matrix(NA, nrow=n_rho, ncol=numIter)
  for (i in 1:n_rho) {
    print(rho[i])
    list[i, ] <- estimateFunc(rho[i], numData, numIter)
  }
  return(list)
}


get_mean_squared_error <- function(rho, list) {
  n_rho <- length(rho)
  numIter <- length(list[1, ])
  MSE <- rep(NA, n_rho)
  for (i in 1:n_rho) {
    MSE[i] <- as.numeric(sum((list[i, ] - rho[i])^2) / numIter)
  }
  return(MSE)
}

get_KL_divergence_error <- function(rho, list) {
  n_rho <- length(rho)
  numIter <- length(list[1, ])
  KL_div <- rep(NA, n_rho)
  for (i in 1:n_rho) {
    KL_div[i] <- 0
    for (j in 1:numIter) {
      KL_div[i] <- KL_div[i] + get_KL_divergence(rho[i], list[i, j])
    }
    KL_div[i] <- KL_div[i] / numIter
  }
  return(KL_div)
}

get_Fisher_information_error <- function(rho, list) {
  n_rho <- length(rho)
  numIter <- length(list[1, ])
  fisher_information <- rep(NA, n_rho)
  for (i in 1:n_rho) {
    fisher_information[i] <- 0
    for (j in 1:numIter) {
      fisher_information[i] <- fisher_information[i] + get_Fisher_information_metric(rho[i], list[i, j])
    }
    fisher_information[i] <- fisher_information[i] / numIter
  }
  return(fisher_information)
}


get_MSE_results <- function(rho, numData, numIter) {
  df <- rbind(
    as.data.frame(
      cbind(rho,
            mse=get_mean_squared_error(function(x1, x2, x3)
              get_Bayes_estimate_distribution(posterior_J, x1, x2, x3), rho, numData, numIter),
            name="MSE w/Jeffreys prior"
      )
    ),
    as.data.frame(
      cbind(rho,
            mse=get_mean_squared_error(function(x1, x2, x3)
              get_Bayes_estimate_distribution(posterior_flat, x1, x2, x3), rho, numData, numIter),
            name="MSE w/Flat prior"
      )
    ),
    as.data.frame(
      cbind(rho,
            mse=get_mean_squared_error(function(x1, x2, x3)
              get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818), x1, x2, x3),
              rho, numData, numIter),
            name="MSE w/PC prior (lambda = 0.0818)"
      )
    ),
    as.data.frame(
      cbind(rho,
            mse=get_mean_squared_error(function(x1, x2, x3)
              get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249), x1, x2, x3),
              rho, numData, numIter),
            name="MSE w/PC prior (lambda = 1.249)"
      )
    ),
    as.data.frame(
      cbind(rho,
            mse=get_mean_squared_error(function(x1, x2, x3)
              get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360), x1, x2, x3),
              rho, numData, numIter),
            name="MSE w/PC prior (lambda = 5.360)"
      )
    ),
    as.data.frame(
      cbind(
        rho,
        value=get_mean_squared_error()
      )
    )
  )

  df$rho <- sapply(df$rho, as.numeric)
  df$mse <- sapply(df$mse, as.numeric)
  
  return(df)
}

plot_results <- function(df, y_label="") {
  gg <- ggplot(data=df, aes(rho, value)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2, group=1) +
    labs(
      y=y_label,
      x="Correlation"
    ) +
    facet_wrap(~ name)
  return(gg)
}



#posteriorFunc <- function(x1, x2, x3, x4) posterior_PC(x1, x2, x3, x4, lambda=5.360)
#posteriorFunc <- posterior_J
#estimateFunc <- function(x1, x2, x3) get_Bayes_estimate_distribution(posteriorFunc, x1, x2, x3)
#estimateFunc <- function(x1, x2, x3) get_Bayes_estimate_with_KL_divergence_distribution(posteriorFunc, x1, x2, x3)
#estimateFunc <- function(x1, x2, x3) get_Bayes_estimate_with_Fisher_information_distribution(posteriorFunc, x1, x2, x3)
#estimateFunc <- get_MLE_distribution
#estimateFunc <- get_empirical_distribution
#estimateFunc <- get_empirical_with_known_means_distribution
knownMean <- TRUE
numSamples <- 1000
estimateFunc <- function(x1, x2, x3) get_Bayes_estimate_fiducial_with_Fisher_information_distribution(x1, x2, x3, numSamples=numSamples, knownMean=knownMean)
rho <- seq(0.1, 0.9, by=0.2)
numData <- 3
numIter <- 10000

estimates <- calculate_estimates(estimateFunc, rho, numData, numIter)
#get_mean_squared_error(rho, estimates)
get_KL_divergence_error(rho, estimates)
get_Fisher_information_error(rho, estimates)
