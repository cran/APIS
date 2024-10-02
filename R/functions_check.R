# THIS FILE CONTAINS ALL THE FUNCTIONS FOR CHECKING INPUTS
#
# =====================================================================================

# GENERAL PART
# ======================================================================================

#' Check matrix format
#'
#' @param genotype_matrix genotype matrix
#' @param verbose verbose
#'
#' @importFrom methods is
#'
#' @keywords internal
#' @noRd

check_matrices_format <- function(genotype_matrix, verbose = FALSE) {
  # DESCRIPTION
  # This function checks the matrices format
  #
  # INPUTS
  # genotype_matrix : genotype matrix in the APIS format
  # verbose : display information on the screen
  #
  # OUTPUTS
  # stop if an error is detected

  # Check if genotype_matrix is a matrix
  if (!is(genotype_matrix, 'matrix')) {
    stop("The genotype matrices should be a 'matrix' object : use the function as.matrix() ")
  }

  # Check if genotype_matrix is a character matrix
  if (!is.character(genotype_matrix)) {
    stop("The genotype matrices should be filled with 'character' values : check the genotype format and use the function storage.mode() ")
  }
  
  sep_verif=grepl(pattern = "/",x = genotype_matrix)
  if (!all(sep_verif)){
    warning(paste0("Unusual genotype found: ",paste(unique(genotype_matrix[!sep_verif]),collapse = ", ")))
    warning("The genotypes must be separated by '/': 'A/A' or 'NA/NA' for diploids and 'A/A/A' or 'NA/NA/NA' for triploids. At least one genotype was missing '/'. See the supplied datasets in APIS for reference.")
  }

  if (verbose) {
    cat(paste0(deparse(substitute(genotype_matrix)), " matrix format : OK"))
    cat("\n")
  }
}

#' Check genotypes
#'
#' @param offspring_genotype offspring genotype matrix
#' @param sire_genotype sire genotype matrix
#' @param dam_genotype dam genotype matrix
#' @param verbose verbose
#'
#' @return list of all genotypes matrices
#'
#' @keywords internal
#' @noRd

check_genotypes <- function(offspring_genotype, sire_genotype, dam_genotype,
                            offspring_ploidy_level = 2,
                            parental_ploidy_level = 2,
                            verbose = FALSE) {
  # DESCRIPTION
  # Check different genotype format errors
  #
  # INPUTS
  # offspring_genotype : genotype matrix of the offspring in the APIS format
  # sire_genotype : genotype matrix of the dams in the APIS format
  # dam_genotype : genotype matrix of the offspring in the APIS format
  # offspring_ploidy_level : ploidy level of the offspring
  # parental_ploidy_level : ploidy level of the parental population
  # verbose : display information on the screen
  #
  # OUTPUTS
  # list of the tree matrices after quality control

  # Check if all matrices have the same number of markers
  if (ncol(offspring_genotype) != ncol(sire_genotype) || ncol(offspring_genotype) != ncol(dam_genotype)) {
    stop("All genotype matrices (offspring_genotype, sire_genotype and dam_genotype) should have the same numbers of markers")
  }

  # Check for offspring genotype marker that is not genotyped
  offspring_ploidy_NA = paste0(rep("NA", times = offspring_ploidy_level), collapse = "/")
  check_NA_offspring <- function(X) {
    if (length(X[which(X == offspring_ploidy_NA)]) == length(X)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  parental_ploidy_NA = paste0(rep("NA", times = parental_ploidy_level), collapse = "/")
  check_NA_parent <- function(X) {
    if (length(X[which(X == parental_ploidy_NA)]) == length(X)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  offspring_marker_genotype <- apply(offspring_genotype, 2, check_NA_offspring)
  sire_marker_genotype <- apply(sire_genotype, 2, check_NA_parent)
  dam_marker_genotype <- apply(dam_genotype, 2, check_NA_parent)

  marker_non_genotyped <- which(offspring_marker_genotype | sire_marker_genotype | dam_marker_genotype)
  if (length(marker_non_genotyped) == 0) {
    # all markers are genotyped
  } else {
    id_marker <- colnames(offspring_genotype)[marker_non_genotyped]
    offspring_genotype  <- offspring_genotype[, -marker_non_genotyped]
    sire_genotype       <- sire_genotype[, -marker_non_genotyped]
    dam_genotype        <- dam_genotype[, -marker_non_genotyped]

    if (verbose) {
      message(paste0("marker(s) : ", paste0(id_marker, collapse = " / "), " with no genotypes were removed"))
    }
  }

  return(list(offspring_genotype, sire_genotype, dam_genotype))
}

#' Check parameters
#'
#' @param exclusion_threshold the exclusion threshold
#' @param number_marker number of markers genotyped
#' @param error maximum error rate tolerated in the pedigree
#' @param verbose verbose
#'
#' @importFrom methods is
#'
#' @keywords internal
#' @noRd

check_input_parameters <- function(exclusion_threshold, number_marker, error,
                                   verbose = FALSE) {
  # DESCRIPTION
  # Check the input parameters : exclusion threshold and maximum error accepted
  #
  # INPUTS
  # exclusion_threshold : exclusion threshold
  # number_marker : number of markers
  # error : percent of error allowed in "mendel" method
  # verbose : display information on the screen
  #
  # OUPUTS
  # stop if errors are detected

  # Check if exclusion_threshold is a numeric variable
  if (!is(exclusion_threshold, 'numeric')) {
    if (! is.null(exclusion_threshold)){
      stop("exclusion_threshold parameter should be a 'numeric variable' or NULL (default) : use the function as.numeric() ")
    }
  } else {
    #	Check if the number of mismatches allowed is lower than the number of markers and positive
    if ((exclusion_threshold >= 0) && (exclusion_threshold <= number_marker)) {
      # OK
    } else {
      stop("The exclusion threshold is not in the correct range : should be in the range [0, number of markers]")
    }
  }

  # Check if the error parameter is a numeric variable
  if (!is(error, 'numeric')) {
    stop("error parameter should be a 'numeric variable' : use the function as.numeric() ")
  }

  #	Check if the error in the range [0, 1]
  if ((error < 0) | (error > 1)) {
    stop("The error is not in the correct range : should be in the range [0, 1]")
  }

  if (verbose) {
    cat('exclusion_threshold parameter : OK', sep = "\n")
    cat('error parameter : OK', sep = "\n")
    cat('---------------------------------------------------', sep = '\n')
  }

}
