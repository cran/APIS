# THIS FILE CONTAINS ALL THE FUNCTIONS FOR CREATING PEDIGREE
#
# =====================================================================================

#' Create pedigree for delta method
#'
#' @param log_file log file
#' @param threshold_delta threshold for delta
#'
#' @return pedigree using the delta method
#'
#' @keywords internal
#' @noRd

create_pedigree_delta <- function(log_file, threshold_delta) {
  # DESCRIPTION
  # create the pedigree file for "mendel" method
  #
  # INPUTS
  # log_file : the log file
  # threshold_delta : the threshold value for deltas
  #
  # OUPUTS
  # pedigree : the pedigree

  pedigree = log_file[,1:3]
  colnames(pedigree) = c('offspring', 'sire', 'dam')

  individual_not_assigned = which(log_file$delta_1_2 < threshold_delta)

  if(length(individual_not_assigned) > 0) {
    pedigree[individual_not_assigned, 2:3] = c(NA, NA)
  }

  return(pedigree)
}

#' Create pedigree for exclusion method
#'
#' @param log_file log file
#' @param threshold_exclusion threshold for mismatches
#'
#' @return pedigree using the exclusion method
#'
#' @keywords internal
#' @noRd

create_pedigree_exclusion <- function(log_file, threshold_exclusion) {
  # DESCRIPTION
  # create the pedigree file for "exclusion" method
  #
  # INPUTS
  # log_file : the log file
  # threshold_exclusion : the threshold value for mismatches
  #
  # OUPUTS
  # pedigree : the pedigree

  pedigree = log_file[,1:3]
  colnames(pedigree) = c('offspring', 'sire', 'dam')

  individual_not_assigned = which((log_file$mismatch_1 > threshold_exclusion) |
                                    (log_file$mismatch_1 == log_file$mismatch_2))

  if(length(individual_not_assigned) > 0) {
    pedigree[individual_not_assigned, 2:3] = c(NA, NA)
  }

  return(pedigree)
}

#' Create pedigree for probability method
#'
#' @param log_file log file
#' @param threshold_probability threshold for probabilities
#'
#' @return pedigree using the probabilty method
#'
#' @keywords internal
#' @noRd

create_pedigree_probability <- function(log_file, threshold_probability) {
  # DESCRIPTION
  # create the pedigree file for "probability" method
  #
  # INPUTS
  # log_file : the log file
  # threshold_exclusion : the threshold value for probabilities
  #
  # OUPUTS
  # pedigree : the pedigree

  pedigree = log_file[,1:3]
  colnames(pedigree) = c('offspring', 'sire', 'dam')

  individual_not_assigned = which(log_file$probability_1 < threshold_probability)

  if(length(individual_not_assigned) > 0) {
    pedigree[individual_not_assigned, 2:3] = c(NA, NA)
  }

  return(pedigree)
}
