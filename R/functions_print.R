# THIS FILE CONTAINS ALL THE FUNCTIONS FOR PRINTING
#
# =====================================================================================


#' Print title of APIS
#'
#' @keywords internal
#' @noRd

print_title <- function() {
  # DESCRIPTION
  # Print the title
  #
  # INPUTS
  #
  # OUPUTS
  # display the title on the screen

  cat('===================================================', sep = '\n')
  cat('           ___   _____   _   _____  ', sep = '\n')
  cat('          /   | |  _  \\ | | /  ___/ ', sep = '\n')
  cat('         / /| | | |_| | | | | |___  ', sep = '\n')
  cat('        / / | | |  ___/ | | \\ __  \\ ', sep = '\n')
  cat('       / /  | | | |     | |  ___| | ', sep = '\n')
  cat('      /_/   |_| |_|     |_| /_____/ ', sep = '\n')
  cat('\n')
  cat('---------------------------------------------------', sep = '\n')
  cat('AUTO-ADAPTIVE PARENTAGE INFERENCE SOFTWARE', sep = '\n')
  cat('---------------------------------------------------', sep = '\n')
}

#' Print the summary of the assignment
#'
#' @keywords internal
#' @noRd

print_summary <- function(theoretical_assignment_power,
                          assignment_rate) {
  # DESCRIPTION
  # Print the assignment summary
  #
  # INPUTS
  # theoretical_assignment_power : theoretical assignment power
  # # assignment_rate : assignment rate

  # OUPUTS
  # display the summary on the screen

  cat('--------------------------------------', sep = '\n')
  cat('             APIS SUMMARY', sep = '\n')
  cat('--------------------------------------', sep = '\n')

  if (!is.null(theoretical_assignment_power)) {
    P = substr(theoretical_assignment_power, 1, 5)
    cat("theoretical assignment power : ", P)
    cat("\n")
  }

  AR = assignment_rate * 100
  cat("realized assignment rate : ", AR, "%")
  cat("\n")

}
