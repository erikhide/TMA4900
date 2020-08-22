simulate_bivariate_normal <- function(rho, numData) {
  data <- mvrnorm(numData, mu=c(0, 0), Sigma=matrix(c(1, rho, rho, 1), 2))
  T1 <- sum(data[,1]^2) + sum(data[,2]^2)
  T2 <- sum(data[,1]*data[,2])
  return(c(T1, T2))
}

normalise_prior <- function(priorFunc) {
  const <- integrate(priorFunc, lower=-1, upper=1)
  return(const$value)
}

normalise_posterior <- function(posteriorFunc, statistic, numData) {
  const <- integrate(function(x) posteriorFunc(x, statistic[1], statistic[2], numData), lower=-1, upper=1)
  return(const$value)
}

normalise_function <- function(density_approximation) {
  const <- integrate(density_approximation, lower=-1, upper=1)
  return(const$value)
}

get_density_approximation <- function(samples) {
  density_approximation <- kdensity(samples, normalized=FALSE)
  normalised_density <- function(x) density_approximation(x) / normalise_function(density_approximation)
  return(density_approximation)
}

get_expected_value_of_density_approximation <- function(density_approximation) {
  expected_value_func <- function(x) x * density_approximation(x)
  return(integrate(expected_value_func, lower=-1, upper=1)$value)
}

get_variance_of_density_approximation <- function(density_approximation) {
  expected_value <- get_expected_value_of_density_approximation(density_approximation)
  variance_func <- function(x) (x - expected_value)^2 * density_approximation(x)
  return(integrate(variance_func, lower=-1, upper=1)$value)
}

plot_histogram <- function(samples) {
  samples <- data.frame(samples)
  colnames(samples) <- c("value")
  gg <- ggplot(samples, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  return(gg)
}

plot_density_approximation <- function(density_approximation) {
  epsilon <- 0.001
  x <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  df <- data.frame(rho=x, density=density_approximation(x))
  df$rho <- sapply(df$rho, as.numeric)
  df$density <- sapply(df$density, as.numeric)
  gg <- ggplot(data=df, aes(x=rho, y=density)) +
    theme_grey(base_size = 22) +
    geom_line(size=1, group=1) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    xlim(c(-1, 1))
  return(gg)
}
