# THIS FILE CONTAINS ALL THE FUNCTIONS FOR ASSIGNMENT RATE
#
# =====================================================================================

#' Assignment rate
#'
#' @param pedigree pedigree file (from the APIS_2n() or APIS_3n() function)
#'
#' @return assignment rate
#' @keywords internal
#' @noRd
#'

get_assignment_rate <- function(pedigree) {
  # DESCRIPTION
  # compute the assignment rate (in %)
  #
  # INPUTS
  # pedigree : the pedigree file
  #
  # OUTPUTS
  # assignment_rate : the assignment rate
  individual_assigned = length(which(!is.na(pedigree$sire)))

  assignment_rate = individual_assigned / nrow(pedigree)
  return(assignment_rate)
}
