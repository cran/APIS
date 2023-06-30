# THIS FILE CONTAINS ALL THE FUNCTIONS FOR PLOTING
#
# =====================================================================================

#' Plot probabilities
#'
#' @param log_file log file from the APIS_2n() or APIS_3n function
#' @param threshold threshold
#' @param simulated_individuals names of the simulated individuals
#'
#' @return plot of the distribution of probabilities
#'
#' @import ggplot2
#'
#' @export

plot_probabilities <- function(log_file, threshold = NULL, simulated_individuals = NULL) {
  # DESCRIPTION
  # plot probabilities histogram
  #
  # INPUTS
  # log_file : the log file
  # threshold : the threshold value
  # simulated_individuals : vector of simulated individual names
  #
  # OUPUTS
  # p : the graph

  if (!is.null(simulated_individuals)) {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$probability_1, log_file$probability_2),
                           P = rep(c('P1', 'P2'), each = nrow(log_file)),
                           simulated = log_file$simulated)
  } else {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$probability_1, log_file$probability_2),
                           P = rep(c('P1', 'P2'), each = nrow(log_file)),
                           simulated = 0)
  }

  p = ggplot2::ggplot(data = data_plot, aes(x = .data$values)) +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "P2"),
                            alpha = 0.25, bins = 30, fill = "blue") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "P1"),
                            alpha = 0.25, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "P1"),
                            alpha = 0.3, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "P2"),
                            alpha = 0.3, bins = 30, fill = "blue") +
    ggplot2::xlab('mendelian transmission probability') +
    ggplot2::theme(axis.title.y = element_blank())

  if (!is.null(threshold)) {
    p = p + ggplot2::geom_vline(xintercept = threshold)
  }

  return(p)
}

#' Plot deltas
#'
#' @param log_file log file from the APIS_2n() or APIS_3n function
#' @param threshold threshold
#' @param simulated_individuals names of the simulated individuals
#'
#' @return plot of the distribution of delta
#'
#' @import ggplot2
#'
#' @export

plot_delta <- function(log_file, threshold = NULL, simulated_individuals = NULL) {
  # DESCRIPTION
  # plot delta histogram
  #
  # INPUTS
  # log_file : the log file
  # threshold : the threshold value
  # simulated_individuals : vector of simulated individual names
  #
  # OUPUTS
  # p : the graph

  if (!is.null(simulated_individuals)) {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$delta_1_2, log_file$delta_2_3),
                           P = rep(c('delta_1_2', 'delta_2_3'), each = nrow(log_file)),
                           simulated = log_file$simulated)
  } else {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$delta_1_2, log_file$delta_2_3),
                           P = rep(c('delta_1_2', 'delta_2_3'), each = nrow(log_file)),
                           simulated = 0)
  }

  p = ggplot2::ggplot(data = data_plot, aes(x = .data$values)) +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "delta_2_3"),
                            alpha = 0.25, bins = 30, fill = "blue") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "delta_1_2"),
                            alpha = 0.25, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "delta_1_2"),
                            alpha = 0.3, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "delta_2_3"),
                            alpha = 0.3, bins = 30, fill = "blue") +
    ggplot2::xlab('delta') +
    ggplot2::theme(axis.title.y = element_blank())

  if (!is.null(threshold)) {
    p = p + ggplot2::geom_vline(xintercept = threshold)
  }

  return(p)
}

#' Plot mismatches
#'
#' @param log_file log file from the APIS_2n() or APIS_3n function
#' @param threshold threshold
#' @param simulated_individuals names of the simulated individuals
#'
#' @return plot of the distribution of mismatches
#'
#' @import ggplot2
#'
#' @export

plot_mismatches <- function(log_file, threshold = NULL, simulated_individuals = NULL) {
  # DESCRIPTION
  # plot mistmaches histogram
  #
  # INPUTS
  # log_file : the log file
  # threshold : the threshold value
  # simulated_individuals : vector of simulated individual names
  #
  # OUPUTS
  # p : the graph

  if (!is.null(simulated_individuals)) {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$mismatch_1, log_file$mismatch_2),
                           P = rep(c('mismatch_1', 'mismatch_2'), each = nrow(log_file)),
                           simulated = log_file$simulated)
  } else {
    log_file$simulated = ifelse(test = log_file$offspring %in% simulated_individuals, yes = 1, no = 0)
    data_plot = data.frame(values = c(log_file$mismatch_1, log_file$mismatch_2),
                           P = rep(c('mismatch_1', 'mismatch_2'), each = nrow(log_file)),
                           simulated = 0)
  }

  p = ggplot2::ggplot(data = data_plot, aes(x = .data$values)) +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "mismatch_2"),
                            alpha = 0.25, bins = 30, fill = "blue") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 1 & data_plot$P == "mismatch_1"),
                            alpha = 0.25, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "mismatch_1"),
                            alpha = 0.3, bins = 30, fill = "red") +
    ggplot2::geom_histogram(data = subset(data_plot, data_plot$simulated == 0 & data_plot$P == "mismatch_2"),
                            alpha = 0.3, bins = 30, fill = "blue") +
    ggplot2::xlab('number of mismatch') +
    ggplot2::theme(axis.title.y = element_blank())

  if (!is.null(threshold)) {
    p = p + ggplot2::geom_vline(xintercept = threshold)
  }

  return(p)
}
