library(ggplot2)
library(ggthemes)
library(kdensity)
source("bivariate_normal.R")
source("util.R")

get_estimates <- function(rho, numData, numIter) {
  n_rho <- length(rho)
  df <- rbind(
    as.data.frame(
      cbind(rho,
            estimate=get_SLR_estimate_distribution(rho, numData, numIter),
            name="SLR"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_MLE_distribution(rho, numData, numIter),
            name="MLE"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_empirical_with_known_variances_distribution(rho, numData, numIter),
            name="Empirical w/known variance"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_empirical_distribution(rho, numData, numIter),
            name="Empirical w/uknown variance"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_Bayes_estimate_distribution(posterior_J, rho, numData, numIter),
            name="MMSE w/Jeffreys prior"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_Bayes_estimate_distribution(posterior_flat, rho, numData, numIter),
            name="MMSE w/Flat prior"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818),
                                                     rho, numData, numIter),
            name="MMSE w/PC prior (lambda = 0.0818)"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249),
                                                     rho, numData, numIter),
            name="MMSE w/PC prior (lambda = 1.249)"
      )
    ),
    as.data.frame(
      cbind(rho,
            estimate=get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360),
                                                     rho, numData, numIter),
            name="MMSE w/PC prior (lambda = 5.360)"
      )
    )
  )
  
  df$rho <- sapply(df$rho, as.numeric)
  df$estimate <- sapply(df$estimate, as.numeric)
  
  return(df)
}

get_estimates_densities <- function(rho, numData, numIter) {
  epsilon <- 0.001
  x <- seq(-1 + epsilon, 1 - epsilon, by=0.01)
  df <- rbind(
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_SLR_estimate_distribution(rho, numData, numIter))(x),
            name="SLR"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_MLE_distribution(rho, numData, numIter))(x),
            name="MLE"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_empirical_with_known_variances_distribution(rho, numData, numIter))(x),
            name="Empirical w/known variance"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_empirical_distribution(rho, numData, numIter))(x),
            name="Empirical w/uknown variance"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_Bayes_estimate_distribution(posterior_J, rho, numData, numIter))(x),
            name="MMSE w/Jeffreys prior"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_Bayes_estimate_distribution(posterior_flat, rho, numData, numIter))(x),
            name="MMSE w/Flat prior"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818),
                                                     rho, numData, numIter))(x),
            name="MMSE w/PC prior (lambda = 0.0818)"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249),
                                                     rho, numData, numIter))(x),
            name="MMSE w/PC prior (lambda = 1.249)"
      )
    ),
    as.data.frame(
      cbind(rho=x,
            estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360),
                                                     rho, numData, numIter))(x),
            name="MMSE w/PC prior (lambda = 5.360)"
      )
    )
  )
  
  df$rho <- sapply(df$rho, as.numeric)
  df$estimate <- sapply(df$estimate, as.numeric)
  
  return(df)
}

plot_estimate_histograms <- function(rho, numData, numIter) {
  df <- get_estimates(rho, numData, numIter)
  gg <- ggplot(data=df, aes(x=estimate)) +
    theme_grey(base_size = 22) +
    geom_histogram(bins=50, aes(y=..density..), group=1) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    xlim(c(-1, 1)) +
    facet_wrap(~ name)
  return(gg)
}

plot_estimate_densities <- function(rho, numData, numIter) {
  df <- get_estimates_densities(rho, numData, numIter)
  gg <- ggplot(data=df, aes(x=rho, y=estimate)) +
    theme_grey(base_size = 22) +
    geom_line(size=1, group=1) +
    labs(
      y="Density",
      x="Correlation"
    ) +
    xlim(c(-1, 1)) +
    facet_wrap(~ name)
  return(gg)
}


get_estimate_densities_Bayes_estimators <- function(rho, numData, numIter) {
  epsilon <- 0.001
  x <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  df <- data.frame()
  for (i in 1:length(rho)) {
    df <- rbind(df,
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_distribution(posterior_J, rho[i], numData, numIter))(x),
                    name="Jeffreys prior",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_distribution(posterior_flat, rho[i], numData, numIter))(x),
                    name="Flat prior",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818),
                                                                                       rho[i], numData, numIter))(x),
                    name="PC prior (0.0818)",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249),
                                                                                       rho[i], numData, numIter))(x),
                    name="PC prior (1.249)",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360),
                                                                                       rho[i], numData, numIter))(x),
                    name="PC prior (5.360)",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_KL_divergence_distribution(posterior_J, rho[i], numData, numIter))(x),
                    name="Jeffreys prior",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_KL_divergence_distribution(posterior_flat, rho[i], numData, numIter))(x),
                    name="Flat prior",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_KL_divergence_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818),
                                                                                                          rho[i], numData, numIter))(x),
                    name="PC prior (0.0818)",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_KL_divergence_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249),
                                                                                                          rho[i], numData, numIter))(x),
                    name="PC prior (1.249)",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_KL_divergence_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360),
                                                                                                          rho[i], numData, numIter))(x),
                    name="PC prior (5.360)",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_Fisher_information_distribution(posterior_J, rho[i], numData, numIter))(x),
                    name="Jeffreys prior",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_Fisher_information_distribution(posterior_flat, rho[i], numData, numIter))(x),
                    name="Flat prior",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_Fisher_information_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=0.0818),
                                                                                                               rho[i], numData, numIter))(x),
                    name="PC prior (0.0818)",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_Fisher_information_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=1.249),
                                                                                                               rho[i], numData, numIter))(x),
                    name="PC prior (1.249)",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_with_Fisher_information_distribution(function(t1, t2, t3, t4) posterior_PC(t1, t2, t3, t4, lambda=5.360),
                                                                                                               rho[i], numData, numIter))(x),
                    name="PC prior (5.360)",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                )
    )
  }
  df$rho <- sapply(df$rho, as.numeric)
  df$estimate <- sapply(df$estimate, as.numeric)
  return(df)
}

plot_estimate_densities_together <- function(df, title="") {
  gg <- ggplot(data=df, aes(x=rho, y=estimate, group=name)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2, aes(linetype = factor(name))) +
    labs(
      y="Density",
      x="Correlation",
      title=title
    ) +
    scale_linetype_manual(name="", values = c("solid", "dashed", "dotdash", "twodash", "dotted")) +
    theme(legend.key.size = unit(4,"line")) +
    theme(legend.position = "top") +
    xlim(c(-1, 1)) +
    ylim(c(0,16)) +
    facet_grid(true_rho ~ loss)
  return(gg)
}


get_estimate_densities_fiducial <- function(rho, numData, numIter, numSamples=1000) {
  epsilon <- 0.001
  x <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  df <- data.frame()
  for (i in 1:length(rho)) {
    df <- rbind(df,
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=FALSE))(x),
                    name="Unknown mean",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=TRUE))(x),
                    name="Known mean",
                    loss="Mean square error",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=FALSE))(x),
                    name="Unknown mean",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=TRUE))(x),
                    name="Known mean",
                    loss="Kullback-Leibler divergence",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=FALSE))(x),
                    name="Known mean",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_Bayes_estimate_fiducial_distribution(rho[i], numData, numIter, numSamples, knownMean=TRUE))(x),
                    name="Unknown mean",
                    loss="Fisher information metric",
                    true_rho=rho[i]
                  )
                )
    )
  }
  df$rho <- sapply(df$rho, as.numeric)
  df$estimate <- sapply(df$estimate, as.numeric)
  return(df)
}

get_estimate_densities_freq <- function(rho, numData, numIter) {
  epsilon <- 0.001
  x <- seq(-1 + epsilon, 1 - epsilon, by=0.001)
  df <- data.frame()
  for (i in 1:length(rho)) {
    df <- rbind(df,
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_MLE_distribution(rho[i], numData, numIter))(x),
                    name="MLE",
                    true_rho=rho[i]
                  )
                ),
                as.data.frame(
                  cbind(
                    rho=x,
                    estimate=get_density_approximation(get_empirical_with_known_means_distribution(rho[i], numData, numIter))(x),
                    name="Empirical correlation w/known mean",
                    true_rho=rho[i]
                  )
                )
    )
  }
  df$rho <- sapply(df$rho, as.numeric)
  df$estimate <- sapply(df$estimate, as.numeric)
  return(df)
}

plot_estimate_densities_together_freq <- function(df, title="") {
  gg <- ggplot(data=df, aes(x=rho, y=estimate, group=name)) +
    theme_grey(base_size = 30) +
    geom_line(size=1.2, aes(linetype = factor(name))) +
    labs(
      y="Density",
      x="Correlation",
      title=title
    ) +
    scale_linetype_manual(name="", values = c("solid", "dashed", "dotdash", "twodash", "dotted")) +
    theme(legend.key.size = unit(4,"line")) +
    theme(legend.position = "top") +
    xlim(c(-1, 1)) +
    #ylim(c(0,16)) +
    facet_wrap(~ true_rho)
  return(gg)
}


rho <- c(0, 0.4, 0.8)
numData <- 3
numIter <- 10000
#df <- get_estimate_densities_fiducial(rho, numData, numIter, numSamples=1000)
df <- get_estimate_densities_freq(rho, numData, numIter)
plot_estimate_densities_together_freq(df, title="")
