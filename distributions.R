source("util.R")

likelihood <- function(rho, T1, T2, n) {
  1 / ((2*pi)^n * (1-rho**2)^(n/2+1)) * exp((-1/2 * T1 + rho * T2) / (1-rho^2))
}

prior_PC <- function(rho, lambda=1) {
  prior <- lambda * rho * sign(rho) / (2 * (1-rho^2) * sqrt(-log(1-rho^2))) * exp(-lambda * sqrt(-log(1-rho^2)))
  for (i in 1:length(rho)) {
    if (rho[i] == 0) {
      prior[i] <- lambda / 2
    }
  }
  return(prior)
}

posterior_PC <- function(rho, T1, T2, n, lambda=1) {
  likelihood(rho, T1, T2, n) * prior_PC(rho, lambda)
}

get_user_defined_scaling_PC <- function(alpha, rho_0) {
  -log(alpha) / sqrt(-log(1-rho_0^2))
}

prior_J <- function(rho) {
  sqrt(1+rho^2) / (1-rho^2)
}

posterior_J <- function(rho, T1, T2, n) {
  likelihood(rho, T1, T2, n) * prior_J(rho)
}

prior_flat <- function(rho) {
  1/2 * (rho > -2)
}

posterior_flat <- function(rho, T1, T2, n) {
  likelihood(rho, T1, T2, n) * prior_flat(rho)
}

prior_extra <- function(rho) {
  1 / (1-rho^2)
}

posterior_extra <- function(rho, T1, T2, n) {
  likelihood(rho, T1, T2, n) * prior_extra(rho)
}

find_reference_prior_numerically <- function(numPoints, numReps, numSamples) {
  epsilon <- 0.001
  lengthBetweenRhos <- (2 - 2*epsilon) / numPoints
  rho <- seq(-1+epsilon, 1-epsilon, by=lengthBetweenRhos)
  density <- rep(NA, numPoints+1)
  for (i in 1:numPoints + 1) {
    density_point <- 0
    for (j in 1:numReps) {
      statistic <- simulate_bivariate_normal(rho[i], numSamples)
      temp_function <- function(x) 1 / ((2*pi)^numSamples * (1-x^2)^(numSamples/2)) * 
        exp(-statistic[1] / (2*(1-x^2)) + x * statistic[2] / (1-x^2))
      temp_const <- integrate(temp_function, lower=-1, upper=1)$value
      density_point <- density_point + log(temp_function(rho[i]) / temp_const)
    }
    density[i] <- exp(density_point / numReps)
  }
  return(as.data.frame(cbind(rho, density)))
}

