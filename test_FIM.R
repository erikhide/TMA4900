library(ggplot2)

get_BE_FIM_value <- function(rho_hat, posteriorFunc, statistic, numData) {
  FIM <- function(x) abs(sqrt(2)*atanh(sqrt(2)*x/sqrt(1+x^2)) - asinh(x) - sqrt(2)*atanh(sqrt(2)*rho_hat/sqrt(1+rho_hat^2)) + asinh(rho_hat))
  FIM_function <- function(x) FIM(x) * posteriorFunc(x, statistic[1], statistic[2], numData)
  return(integrate(FIM_function, lower=-1, upper=1)$value)
}

plot_BE_FIM <- function(posteriorFunc, statistic, numData) {
  rho_hat <- seq(-0.999, 0.999, by=0.001)
  numPoints <- length(rho_hat)
  value <- rep(NA, numPoints)
  for (i in 1:numPoints) {
    value[i] <- get_BE_FIM_value(rho_hat[i], posteriorFunc, statistic, numData)
  }
  df <- data.frame(rho_hat=rho_hat, value=value)
  gg <- ggplot(df, aes(x=rho_hat, y=value)) +
    geom_line()
  return(gg)
}


posteriorFunc <- posterior_J
rho <- 0
numData <- 3
statistic <- simulate_bivariate_normal(rho, numData)

plot_BE_FIM(posteriorFunc, statistic, numData)

test <- function(x) get_BE_FIM_value(x, posteriorFunc, statistic, numData)
print(optimize(test, interval=c(-1, 1), maximum=FALSE))
