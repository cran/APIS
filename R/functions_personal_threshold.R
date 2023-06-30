# THIS FILE CONTAINS ALL THE FUNCTIONS FOR PERSONAL THRESHOLD
#
# =====================================================================================

#' Personal threshold
#'
#' @param log_file log file from the APIS_2n() or APIS_3n function
#' @param method method : "delta", "probability" or "exclusion"
#' @param threshold threshold
#' @param verbose verbose
#'
#' @return pedigree file
#'
#' @importFrom cowplot plot_grid
#'
#' @keywords internal
#' @noRd
#'

personal_threshold <- function(log_file, method, threshold,
                               verbose = FALSE) {
  # DESCRIPTION
  # process APIS outputs with new threshold value
  #
  # INPUTS
  # log_file : the log file
  # method : method of assignment
  # threshold : the new threshold value
  # verbose : display information on screen
  #
  # OUPUTS
  # pedigree : the pedigree

  if(method == 'delta') {
    pedigree = create_pedigree_delta(log_file = log_file, threshold_delta = threshold)

    # Plot graphs

    p1 = plot_delta(log_file = log_file, threshold = threshold)
    p2 = plot_probabilities(log_file = log_file, threshold = NULL)
    p3 = plot_mismatches(log_file = log_file, threshold = NULL)
  } else if(method == 'probability') {
    pedigree = create_pedigree_probability(log_file = log_file, threshold_probability = threshold)

    # Plot graphs

    p1 = plot_delta(log_file = log_file, threshold = NULL)
    p2 = plot_probabilities(log_file = log_file, threshold = threshold)
    p3 = plot_mismatches(log_file = log_file, threshold = NULL)
  } else if(method == 'exclusion') {
    pedigree = create_pedigree_exclusion(log_file = log_file, threshold_exclusion = threshold)

    # Plot graphs

    p1 = plot_delta(log_file = log_file, threshold = NULL)
    p2 = plot_probabilities(log_file = log_file, threshold = NULL)
    p3 = plot_mismatches(log_file = log_file, threshold = threshold)
  } else {
    stop("Invalid method : must be 'delta', 'probability' or 'exclusion'")
  }

  print(cowplot::plot_grid(p1, p2, p3, nrow = 3))

  if (verbose) {
    print_summary(theoretical_assignment_power = NULL,
                  assignment_rate = get_assignment_rate(pedigree))
  }

  return(pedigree)
}
