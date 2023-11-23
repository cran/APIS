# THIS FILE CONTAINS ALL THE FUNCTIONS FOR ASSIGNMENT
#
# =====================================================================================

#' Estimate threshold
#'
#' @param log_file log file
#' @param error error
#'
#' @return threshold for delta
#'
#' @importFrom stats quantile
#'
#' @keywords internal
#' @noRd

estimate_mendel_threshold <- function(log_file, error) {
  # DESCRIPTION
  # estimate the delta mendel threshold for assignment
  #
  # INPUTS
  # log_file : the log file
  # error : percent of error allowed in "mendel" method
  #
  # OUTPUTS
  # threshold_delta : value of the delta threshold

  missing_range <- estimate_missing_parents(log_file = log_file)

  # If the number of individuals with missing parents is lower than the error
  if ((missing_range[1]) <= round(error * nrow(log_file))) {
    # All the individuals are assigned to their most likely parent pair

    threshold_delta = min(log_file$delta_1_2)
  } else {
    sorted_delta_2_3 <- sort(log_file$delta_2_3, decreasing = TRUE)

    threshold_delta <- quantile(sorted_delta_2_3[1:(nrow(log_file) - (missing_range[2]))], probs = (1 - error), type = 5, na.rm = TRUE)
  }

  return(threshold_delta)

}

#' Estimate missing parents
#'
#' @param log_file log file
#'
#' @return range of individuals with missing parents
#'
#' @importFrom stats median
#'
#' @keywords internal
#' @noRd

estimate_missing_parents <- function(log_file) {
  # DESCRIPTION
  # estimate the number of individuals with missing parents
  #
  # INPUTS
  # log_file : the log file
  #
  # OUTPUTS
  # vector of 2 elements : estimated number of offspring with missing parents with strict and relaxed method

  # median_p2 <- median(log_file$probability_2)
  #
  # probability_2_sort <- sort(log_file$probability_2)
  #
  # missing_0         <- 2 * length(log_file$probability_1[which(log_file$probability_1 <= median_p2)])
  # number_individual <- length(log_file$probability_1)
  #
  # threshold <- min(log_file$probability_2)
  #
  # N2_l <- round(length(log_file$probability_2[which(log_file$probability_2 <= threshold)]) -
  #                 (length(log_file$probability_3[which(log_file$probability_3 <= threshold)]) * ((missing_0) / number_individual)))
  #
  # cpt <- 1
  # while(N2_l < ((number_individual - (missing_0)) / 2) & cpt <= nrow(log_file)) {
  #   threshold <- probability_2_sort[cpt]
  #
  #   N2_l <- round(length(log_file$probability_2[which(log_file$probability_2 <= threshold)]) -
  #                   (length(log_file$probability_3[which(log_file$probability_3 <= threshold)]) * ((missing_0) / number_individual)))
  #
  #   cpt <- cpt+1
  # }
  #
  # median_1 <- threshold
  # missing_1 <- 2 * length(log_file$probability_1[which(log_file$probability_1 <= median_1)])

  N1_1min <- 2 * length(which(log_file$probability_1 <= median(log_file$probability_2)))
  # N1_1 <- ifelse(test = missing_1 > round(nrow(log_file)/2), yes = round(nrow(log_file)/2), no = missing_1)
  N1_1 <- ifelse(test = N1_1min > nrow(log_file), yes = nrow(log_file), no = N1_1min)

  return(c(N1_1min, N1_1))
}

# ggplot(data=log_file)+geom_histogram(aes(x=probability_2),fill='blue',alpha=0.5)+geom_histogram(aes(x=probability_1),fill='pink',alpha=0.5)+geom_histogram(aes(x=probability_3),fill='green',alpha=0.5)+
#   geom_vline(xintercept=threshold)+geom_vline(xintercept = median_p2,col='red')


#' Estimate threshold
#'
#' @param log_file log file
#' @param error error
#'
#' @return threshold for exclusion
#'
#' @importFrom stats quantile median
#'
#' @keywords internal
#' @noRd

estimate_exclusion_threshold <- function(log_file, error) {
  # DESCRIPTION
  # estimate the delta mendel threshold for assignment
  #
  # INPUTS
  # log_file : the log file
  # error : percent of error allowed in "mendel" method
  #
  # OUTPUTS
  # threshold_exclu : value of the delta threshold

  missing_0 = 2 * length(which(log_file$mismatch_1 >= median(log_file$mismatch_2)))

  missing_range = ifelse(test = missing_0 > nrow(log_file),yes = nrow(log_file),no = missing_0)

  # If the number of individuals with missing parents is lower than the error
  if (missing_range <= round(error * nrow(log_file))) {
    # All the individuals are assigned to their most likely parent pair

    threshold_exclu = max(log_file$mismatch_1)
  } else {
    sorted_miss_2 <- sort(log_file$mismatch_2, decreasing = FALSE) # different from proba_mendel because axis is inversed

    threshold_exclu <- quantile(sorted_miss_2[1:(nrow(log_file) - (missing_range))], probs = error, type = 5, na.rm = TRUE)

  }
  return(threshold_exclu)
}
