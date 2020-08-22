library(plotly)
library(raster)
source("util.R")

get_LRT_function <- function(rho_0, statistic, numData) {
  rho_MLE <- get_MLE(statistic, numData)
  T1 <- statistic[1]
  T2 <- statistic[2]
  return(((1-rho_MLE^2) / (1-rho_0^2))^(numData/2) * 
           exp(-T1/2 * (1/(1-rho_0^2) - 1/(1-rho_MLE^2)) + T2 * (rho_0/(1-rho_0^2) - rho_MLE/(1-rho_MLE^2))))
}

plot_LRT_function <- function(rho_0, numData) {
  T1 <- seq(0, 1000, by=1)
  T2 <- seq(-200, 200, by=1)
  mat <- array(NA, c(length(T1), length(T2)))
  for (i in 1:length(T1)) {
    for (j in 1:length(T2)) {
      mat[i, j] <- get_LRT_function(rho_0, c(T1[i], T2[j]), numData)
    }
  }
  return(mat)
}

mat <- plot_LRT_function(0, 3)
plot(raster(mat))
