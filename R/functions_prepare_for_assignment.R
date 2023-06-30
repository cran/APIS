# THIS FILE CONTAINS ALL THE FUNCTIONS FOR RECODING FOR C FUNCTION
#
# =====================================================================================

#' Genotypes recoding
#'
#' @param genotype_matrix genotype matrix
#' @param variant_dataframe output of the get_variant_dataframe() function
#' @param number_cores number of cores for parallel functions
#'
#' @return recoded genotype matrix
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#'
#' @keywords internal
#' @noRd

recode_genotypes_for_assignment_2n <- function(genotype_matrix, variant_dataframe, number_cores) {
  # DESCRIPTION
  # Recode a genotype matrix according to the variant coding for diploid
  #
  # INPUTS
  # genotype_matrix : genotype matrix in the APIS format
  # variant_dataframe : the variant dataframe (output of get_variant_dataframe())
  # number_cores : number of cores for parallel programming
  #
  # OUPUTS
  # recode_genotype : recoded genotype matrix

  # Recode for C subroutine
  recode_genotype_for_subroutine <- function(marker_genotype, marker_list) {
    tmp <- unlist(strsplit(marker_genotype, '/'))
    all1 <- marker_list[which(marker_list[,1] == tmp[1]),2]
    all2 <- marker_list[which(marker_list[,1] == tmp[2]),2]

    return(sort(c(all1, all2), decreasing = FALSE))
  }

  # Initialize variables
  i <- NULL

  # Create parallel nodes
  cl = parallel::makeCluster(number_cores)
  doParallel::registerDoParallel(cl)

  # Recode offspring genotypes
  recode_genotype <- foreach(i = 1:nrow(genotype_matrix), .combine = rbind) %dopar% {
    tmp <- as.numeric(as.vector(sapply(genotype_matrix[i,, drop = FALSE], recode_genotype_for_subroutine, marker_list = variant_dataframe)))
  }

  parallel::stopCluster(cl)

  if (!is(recode_genotype, "matrix")) {
    recode_genotype <- t(as.matrix(recode_genotype))
  } else {
  }

  rownames(recode_genotype) <- rownames(genotype_matrix)
  return(recode_genotype)
}

#' Genotypes recoding
#'
#' @param genotype_matrix genotype matrix
#' @param variant_dataframe output of the get_variant_dataframe() function
#' @param number_cores number of cores for parallel functions
#'
#' @return recoded genotype matrix
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#'
#' @keywords internal
#' @noRd

recode_genotypes_for_assignment_3n <- function(genotype_matrix, variant_dataframe, number_cores) {
  # DESCRIPTION
  # Recode a genotype matrix according to the variant coding for triploid
  #
  # INPUTS
  # genotype_matrix : genotype matrix in the APIS format
  # variant_dataframe : the variant dataframe (output of get_variant_dataframe())
  # number_cores : number of cores for parallel programming
  #
  # OUPUTS
  # recode_genotype : recoded genotype matrix

  # Recode for C subroutine
  recode_genotype_3n_for_subroutine <- function(marker_genotype, marker_list) {
    tmp <- unlist(strsplit(marker_genotype, '/'))
    all1 <- marker_list[which(marker_list[,1] == tmp[1]),2]
    all2 <- marker_list[which(marker_list[,1] == tmp[2]),2]
    all3 <- marker_list[which(marker_list[,1] == tmp[3]),2]

    new_genotype = c(all1, all2, all3)

    return(sort(new_genotype, decreasing = FALSE))
  }

  # Initialize variables
  i <- NULL

  # Create parallel nodes
  cl = parallel::makeCluster(number_cores)
  doParallel::registerDoParallel(cl)

  # Recode offspring genotypes
  recode_genotype <- foreach(i = 1:nrow(genotype_matrix), .combine = rbind) %dopar% {
    tmp <- as.numeric(as.vector(sapply(genotype_matrix[i,, drop = FALSE], recode_genotype_3n_for_subroutine, marker_list = variant_dataframe)))
  }

  parallel::stopCluster(cl)

  if (!is(recode_genotype, "matrix")) {
    recode_genotype <- t(as.matrix(recode_genotype))
  } else {
  }

  rownames(recode_genotype) <- rownames(genotype_matrix)
  return(recode_genotype)
}

#' Recode allele frequencies
#'
#' @param allele_frequency allele frequencies from the "get_allele_frequencies_2n()" function
#' @param variant_dataframe output of the get_variant_dataframe() function
#'
#' @return recoded allele frequencies
#'
#' @keywords internal
#' @noRd

recode_allele_frequencies_for_assignment <- function(allele_frequency, variant_dataframe) {
  # DESCRIPTION
  # Recode allele frequencies to match the variant coding
  #
  # INPUTS
  # allele_frequency : dataframe of allele frequencies (output of the get_allele_frequencies())
  # variant_dataframe : the variant dataframe (output of get_variant_dataframe())
  #
  # OUPUTS
  # allele_frequency : recoded allele frequencies dataframe

  # Recode function for each marker
  recode_allele_frequencies <- function(column_name, marker_list) {
    tmp <- unlist(strsplit(column_name, '_'))[2]
    return(marker_list[which(marker_list[,1] == tmp), 2])
  }

  # Remove frequency of "NA" variant
  has_NA <- which(colnames(allele_frequency) == "Freq_NA")
  if (length(has_NA) > 0) {
    allele_frequency <- allele_frequency[,-has_NA]
  }

  # Keep only "Freq" columns
  allele_frequency <- allele_frequency[, c((floor(ncol(allele_frequency)/2)+2):ncol(allele_frequency))]

  # Recode columns with the new variant coding system
  colnames(allele_frequency) <- sapply(colnames(allele_frequency),
                                       recode_allele_frequencies,
                                       marker_list = variant_dataframe)
  allele_frequency <- rbind(colnames(allele_frequency), allele_frequency)

  # Add allele frequencies of missing alleles
  add_missing_allele_frequency <- as.numeric(variant_dataframe$recode[-which(variant_dataframe$recode %in% allele_frequency[1,])])
  add_missing_allele_frequency <- add_missing_allele_frequency[-which(add_missing_allele_frequency == 0)]
  add_missing_allele_frequency_matrix <- matrix(data = c(add_missing_allele_frequency, rep(0, times = length(add_missing_allele_frequency) * (nrow(allele_frequency) - 1))),
                                                ncol = length(add_missing_allele_frequency), nrow = nrow(allele_frequency), byrow = TRUE)


  allele_frequency <- apply(allele_frequency, 2, as.numeric)
  allele_frequency <- cbind(allele_frequency, add_missing_allele_frequency_matrix)
  allele_frequency <- allele_frequency[, order(allele_frequency[1,])]
  allele_frequency <- allele_frequency[-1,]

  return(allele_frequency)
}

#' Prepare data for assignment
#'
#' @param offspring_genotype offspring genotype matrix
#' @param sire_genotype sire genotype matrix
#' @param dam_genotype dam genotype matrix
#' @param ploidy_level ploidy level
#' @param number_cores number of cores for parallel functions
#' @param verbose verbose
#'
#' @return recoded datasets
#'
#' @keywords internal
#' @noRd

prepare_for_assignment <- function(offspring_genotype, sire_genotype, dam_genotype,
                                   ploidy_level,
                                   number_cores,
                                   verbose = FALSE) {
  # DESCRIPTION
  # Recode all genotype matrices for assignment routine
  #
  # INPUTS
  # offspring_genotype : genotype matrix of the offspring in the APIS format
  # sire_genotype : genotype matrix of the sires in the APIS format
  # dam_genotype : genotype matrix of the dams in the APIS format
  # ploidy_level : level of ploidy of the offspring
  # number_cores : number of cores for parallel programming
  # verbose : display information of the screen
  #
  # OUPUTS
  # list of 4 elements : recoded matrix of genotypes (offspring, sire, dam) and recoded allele frequencies

  # Estimation of allele frequencies
  allele_frequency = get_allele_frequencies(genotype_matrix = offspring_genotype,
                                            ploidy_level = ploidy_level)

  # Create the variant recoding dataframe
  variant_recode_dataframe = get_variant_dataframe(offspring_genotype = offspring_genotype,
                                                   sire_genotype = sire_genotype,
                                                   dam_genotype = dam_genotype)

  # Recode genotyping matrices
  if (ploidy_level == 2) {
    if (verbose) {
      cat("recoding offspring genotype matrix", sep = "\n")
    }
    offspring_genotypes_recoded = recode_genotypes_for_assignment_2n(offspring_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)
    if (verbose) {
      cat("recoding sire genotype matrix", sep = "\n")
    }
    sire_genotypes_recoded      = recode_genotypes_for_assignment_2n(sire_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)

    if (verbose) {
      cat("recoding dam genotype matrix", sep = "\n")
    }
    dam_genotypes_recoded       = recode_genotypes_for_assignment_2n(dam_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)

  } else if (ploidy_level == 3) {
    if (verbose) {
      cat("recoding offspring genotype matrix", sep = "\n")
    }
    offspring_genotypes_recoded = recode_genotypes_for_assignment_3n(offspring_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)
    if (verbose) {
      cat("recoding sire genotype matrix", sep = "\n")
    }
    sire_genotypes_recoded      = recode_genotypes_for_assignment_2n(sire_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)

    if (verbose) {
      cat("recoding dam genotype matrix", sep = "\n")
    }
    dam_genotypes_recoded       = recode_genotypes_for_assignment_2n(dam_genotype,
                                                                     variant_dataframe = variant_recode_dataframe,
                                                                     number_cores = number_cores)
  } else {

  }


  # Recoded allele frequencies
  if (verbose) {
    cat("recoding allele frequencies", sep = "\n")
    cat('---------------------------------------------------', sep = '\n')
  }
  allele_frequency_recoded = recode_allele_frequencies_for_assignment(allele_frequency = allele_frequency,
                                                                      variant_dataframe = variant_recode_dataframe)

  # Return
  return(list(offspring_genotypes_recoded, sire_genotypes_recoded, dam_genotypes_recoded, allele_frequency_recoded))
}
