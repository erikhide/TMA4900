library(ggplot2)
library(ggthemes)
library(gridExtra)
library(MASS)
source("distributions.R")
source("util.R")

plot_prior <- function(priorFunc, scaling=1) {
  epsilon <- 0.001
  rho <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  density <- scaling * priorFunc(rho)
  df <- as.data.frame(cbind(rho, density))
  gg <- ggplot(df, aes(x=rho, y=density)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2) +
    xlab("Correlation") +
    ylab("Density") +
    ylim(0, 2)
  gg
}

plot_PC_prior <- function(lambda_values=c(0.0818, 1.249, 5.360)) {
  epsilon <- 0.001
  rho <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  if (length(lambda_values) != 3) {
    return(FALSE)
  }
  
  density1 <- prior_PC(rho, lambda=lambda_values[1])
  density2 <- prior_PC(rho, lambda=lambda_values[2])
  density3 <- prior_PC(rho, lambda=lambda_values[3])
  index <- rep(c(1, 2, 3), each=length(rho))
  df <- as.data.frame(cbind(rho, c(density1, density2, density3), index))
  colnames(df)[colnames(df) == "V2"] <- "densities"
  legend_title1 <- expression(paste(" ", lambda, " = 0.0818"))
  legend_title2 <- expression(paste(" ", lambda, " = 1.249"))
  legend_title3 <- expression(paste(" ", lambda, " = 5.360"))
  gg <- ggplot(df, aes(x=rho, y=densities, group=index)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.1, aes(linetype = factor(index))) +
    labs(
      x = "Correlation",
      y = "Density",
      colour = "Priors"
    ) +
    scale_linetype_manual(name="", labels = c(legend_title1, legend_title2, legend_title3), 
                          values = c("solid", "twodash", "dotted")) +
    theme(legend.position = "top") +
    theme(legend.background = element_rect(fill="grey90", 
                                           size=0.5, linetype="solid")) +
    theme(legend.key.size = unit(3,"line")) +
    ylim(0, 3)
  gg
}

plot_PC_prior_limit <- function(lambda_values=c(0.01, 0.001, 0.0001)) {
  epsilon <- 0.0000001
  rho <- seq(0.999, 1 - epsilon, by=epsilon)
  if (length(lambda_values) != 3) {
    return(FALSE)
  }
  
  density1 <- prior_PC(rho, lambda=lambda_values[1])
  density2 <- prior_PC(rho, lambda=lambda_values[2])
  density3 <- prior_PC(rho, lambda=lambda_values[3])
  index <- rep(c(1, 2, 3), each=length(rho))
  df <- as.data.frame(cbind(rho, c(density1, density2, density3), index))
  colnames(df)[colnames(df) == "V2"] <- "densities"
  legend_title1 <- expression(paste(" ", lambda, " = 0.01"))
  legend_title2 <- expression(paste(" ", lambda, " = 0.001"))
  legend_title3 <- expression(paste(" ", lambda, " = 0.0001"))
  gg <- ggplot(df, aes(x=rho, y=densities, group=index)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.1, aes(linetype = factor(index))) +
    labs(
      x = "Correlation",
      y = "Density",
      colour = "Priors"
    ) +
    scale_linetype_manual(name="", labels = c(legend_title1, legend_title2, legend_title3), 
                          values = c("solid", "twodash", "dotted")) +
    theme(legend.position = "top") +
    theme(legend.background = element_rect(fill="grey90", 
                                           size=0.5, linetype="solid")) +
    theme(legend.key.size = unit(3,"line")) +
    ylim(0, 6)
  gg
}

plot_prior_reference <- function(numPoints, numReps, numSamples) {
  df <- find_reference_prior_numerically(numPoints, numReps, numSamples)
  df$J <- 2.45 * prior_J(df$rho)
  gg <- ggplot(df) +
    theme_grey(base_size = 30) +
    geom_point(aes(x=rho, y=density), size=1.2) +
    geom_line(aes(x=rho, y=J), size=1.2) +
    xlab("Correlation") +
    ylab("Density") +
    ylim(0, 15)
  gg
}

plot_priors <- function() {
  epsilon <- 0.001
  rho <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  density_flat <- prior_flat(rho)
  density_J <- prior_J(rho) / 2
  density_PC1 <- prior_PC(rho, lambda=0.0818)
  density_PC2 <- prior_PC(rho, lambda=1.249)
  density_PC3 <- prior_PC(rho, lambda=5.360)
  index <- rep(c(1, 2, 3, 4, 5), each=length(rho))
  df <- as.data.frame(cbind(rho, c(density_flat, density_J, density_PC1, density_PC2, density_PC3), index))
  colnames(df)[colnames(df) == "V2"] <- "densities"
  
  legend_title1 <- expression(paste(" PC, ", lambda, " = 0.0818"))
  legend_title2 <- expression(paste(" PC, ", lambda, " = 1.249"))
  legend_title3 <- expression(paste(" PC, ", lambda, " = 5.360"))
  gg <- ggplot(df, aes(x=rho, y=densities, group=index)) +
    theme_grey(base_size = 26) +
    geom_line(size=1.1, aes(linetype = factor(index))) +
    labs(
      x = "Correlation",
      y = "Density",
      colour = "Priors"
    ) +
    scale_linetype_manual(name="", labels = c("Flat", "Jeffreys", legend_title1, legend_title2, legend_title3), 
                          values = c("solid", "dashed", "dotdash", "twodash", "dotted")) +
    theme(legend.position = "top") +
    theme(legend.key.size = unit(4,"line")) +
    ylim(0, 3)
  gg
}


plot_posterior <- function(posteriorFunc, statistic, numData) {
  rho <- seq(-1, 1, by=0.001)
  if (length(statistic) > 2) {
    density <- posteriorFunc(rho, statistic[1,1], statistic[1,2], numData) / 
      normalise_posterior(posteriorFunc, statistic[1, ], numData)
    for (i in 2:length(statistic[, 1])) {
      density <- density + posteriorFunc(rho, statistic[i, 1], statistic[i, 2], numData) / 
        normalise_posterior(posteriorFunc, statistic[i, ], numData)
    }
    density <- density / length(statistic[, 1])
  }
  else {
    density <- posteriorFunc(rho, statistic[1], statistic[2], numData) / 
      normalise_posterior(posteriorFunc, statistic, numData)
  }
  
  df <- as.data.frame(cbind(rho, density))
  gg <- ggplot(df, aes(x=rho, y=density)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2) +
    xlab("Correlation") +
    ylab("Density")
  return(gg)
}

plot_posterior_PC <- function(statistic, numData, lambda_list) {
  rho <- seq(-1, 1, by=0.001)
  density <- array(NA, c(length(rho), 3))
  for (j in 1:3) {
    posteriorFunc <- function(x1, x2, x3, x4) posterior_PC(x1, x2, x3, x4, lambda=lambda_list[j])
    density[, j] <- posteriorFunc(rho, statistic[1,1], statistic[1,2], numData) / 
      normalise_posterior(posteriorFunc, statistic[1, ], numData)
    for (i in 2:length(statistic[, 1])) {
      density[, j] <- density[, j] + posteriorFunc(rho, statistic[i, 1], statistic[i, 2], numData) / 
        normalise_posterior(posteriorFunc, statistic[i, ], numData)
    }
  }
  density <- density / length(statistic[, 1])

  index <- rep(c(0.0818, 1.249, 5.360), each=length(rho))
  df <- as.data.frame(cbind(rho, c(density[, 1], density[, 2], density[, 3]), index))
  colnames(df)[colnames(df) == "V2"] <- "densities"
  gg <- ggplot(df, aes(x=rho, y=densities, group=index)) +
    theme_grey(base_size = 22) +
    geom_line(size=1.2, aes(linetype = factor(index))) +
    xlab("Correlation") +
    ylab("Density") +
    scale_linetype_manual(name="PC priors", labels = c("legend_title1", "legend_title2", "legend_title3"), 
                          values = c("solid", "twodash", "dotted")) +
    theme(legend.position = c(.75, .8))
  return(gg)
}

plot_posteriors_together <- function(gg_list) {
  df1 <- gg_list[[1]]$data
  df2 <- gg_list[[2]]$data
  df3 <- gg_list[[3]]$data
  df1$num_data <- rep(3, length(df1$rho))
  df2$num_data <- rep(10, length(df1$rho))
  df3$num_data <- rep(100, length(df1$rho))
  df <- rbind(df1, df2, df3)
  gg <- ggplot(data=df, aes(rho, density)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    #ylim(c(0.94, 0.96)) +
    facet_wrap(~ num_data)
  return(gg)
}

plot_posteriors_together_PC <- function(gg_list) {
  df1 <- gg_list[[1]]$data
  df2 <- gg_list[[2]]$data
  df3 <- gg_list[[3]]$data
  df1$num_data <- rep(3, length(df1$rho))
  df2$num_data <- rep(10, length(df1$rho))
  df3$num_data <- rep(100, length(df1$rho))
  df <- rbind(df1, df2, df3)
  legend_title1 <- expression(paste(" ", lambda, " = 0.0818"))
  legend_title2 <- expression(paste(" ", lambda, " = 1.249"))
  legend_title3 <- expression(paste(" ", lambda, " = 5.360"))
  gg <- ggplot(data=df, aes(x=rho, y=densities, group=index)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2, aes(linetype = factor(index))) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    scale_linetype_manual(name="", labels = c(legend_title1, legend_title2, legend_title3), 
                          values = c("solid", "twodash", "dotted")) +
    theme(legend.position = "top") +
    theme(legend.background = element_rect(fill="grey90", 
                                           size=0.5, linetype="solid")) +
    theme(legend.key.size = unit(3,"line")) +
    facet_grid(. ~ num_data)
  return(gg)
}


get_MLE <- function(statistic, numData) {
  MLE <- optimize(function(x) likelihood(x, statistic[1], statistic[2], numData), interval=c(-1, 1), maximum=TRUE)
  return(MLE$maximum)
}

get_MLE_distribution <- function(rho, numData, numIter) {
  MLE_list <- rep(NA, numIter)
  for (i in 1:numIter) {
    statistic <- simulate_bivariate_normal(rho, numData)
    MLE_list[i] <- get_MLE(statistic, numData)
  }
  return(MLE_list)
}

plot_MLE_distribution <- function(rho, numData, numIter) {
  estimate_list <- data.frame(get_MLE_distribution(rho, numData, numIter))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


get_empirical_with_known_variances <- function(rho, numData) {
  data <- mvrnorm(numData, mu=c(0, 0), Sigma=matrix(c(1, rho, rho, 1), 2))
  empirical <- sum(data[, 1] * data[, 2]) / numData
  if (empirical > 1) return(1)
  if (empirical < -1) return(-1)
  return(empirical)
}

get_empirical_with_known_variances_distribution <- function(rho, numData, numIter) {
  list <- rep(NA, numIter)
  for (i in 1:numIter) {
    list[i] <- get_empirical_with_known_variances(rho, numData)
  }
  return(list)
}

plot_empirical_with_known_variances_distribution <- function(rho, numData, numIter) {
  estimate_list <- data.frame(get_empirical_with_known_variances_distribution(rho, numData, numIter))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


get_empirical_with_known_means <- function(rho, numData) {
  data <- mvrnorm(numData, mu=c(0, 0), Sigma=matrix(c(1, rho, rho, 1), 2))
  empirical <- sum(data[, 1] * data[, 2]) / (sqrt(sum(data[, 1]^2) * sum(data[, 2]^2)))
  return(empirical)
}

get_empirical_with_known_means_distribution <- function(rho, numData, numIter) {
  list <- rep(NA, numIter)
  for (i in 1:numIter) {
    list[i] <- get_empirical_with_known_means(rho, numData)
  }
  return(list)
}


get_empirical <- function(rho, numData) {
  data <- mvrnorm(numData, mu=c(0, 0), Sigma=matrix(c(1, rho, rho, 1), 2))
  empirical <- sum((data[, 1]-mean(data[, 1])) * (data[, 2]-mean(data[, 2]))) / 
                (sqrt(sum((data[, 1]-mean(data[, 1]))^2) * sum((data[, 2]-mean(data[, 2]))^2)))
  return(empirical)
}

get_empirical_distribution <- function(rho, numData, numIter) {
  list <- rep(NA, numIter)
  for (i in 1:numIter) {
    list[i] <- get_empirical(rho, numData)
  }
  return(list)
}

plot_empirical_distribution <- function(rho, numData, numIter) {
  estimate_list <- data.frame(get_empirical_distribution(rho, numData, numIter))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


get_SLR_estimate <- function(rho, numData) {
  data <- mvrnorm(numData, mu=c(0, 0), Sigma=matrix(c(1, rho, rho, 1), 2))
  y <- data[, 1]
  x <- data[, 2]
  SLR_model <- lm(y ~ 0 + x)
  return(as.numeric(SLR_model$coefficients))
}

get_SLR_estimate_distribution <- function(rho, numData, numIter) {
  SLR_list <- rep(NA, numIter)
  for (i in 1:numIter) {
    SLR_list[i] <- get_SLR_estimate(rho, numData)
  }
  return(SLR_list)
}

plot_SLR_estimate_distribution <- function(rho, numData, numIter) {
  estimate_list <- data.frame(get_SLR_estimate_distribution(rho, numData, numIter))
  colnames(estimate_list) <- c("value")
  gg <- ggplot(estimate_list, aes(x=value)) +
    geom_histogram(bins=50, aes(y=..density..)) +
    xlab("Correlation") +
    ylab("Density")
  gg
}


rho <- 0.5
numData <- c(3, 10, 100)
posteriorFunc <- posterior_flat
#posteriorFunc <- function(x1, x2, x3, x4) posterior_PC(x1, x2, x3, x4, lambda=0.0818)

gg <- list()

for (i in 1:length(numData)) {
  #statistic <- simulate_bivariate_normal(rho, numData[i])
  numPosteriors <- 1000
  statistic <- array(numeric(), c(numPosteriors, 2))
  for (j in 1:numPosteriors) {
    statistic[j,] <- simulate_bivariate_normal(rho, numData[i])
  }
  #gg[[i]] <- eval(substitute(plot_posterior_PC(statistic, numData[i], c(0.0818, 1.249, 5.360)), list(i = i)))
  gg[[i]] <- eval(substitute(plot_posterior(posteriorFunc, statistic, numData[i]), list(i = i)))
}

plot_posteriors_together(gg)
