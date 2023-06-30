# THIS FILE CONTAINS ALL THE FUNCTIONS FOR ESTIMATING PROBABILITIES
#
# =====================================================================================
#' Estimate probabilities
#'
#' @param offspring_genotype matrix of the offspring genotypes
#' @param sire_genotype matrix of the sire genotypes
#' @param dam_genotype matrix of the dam genotypes
#' @param allele_frequencies allele frequencies (from get_allele_frequencies() function)
#' @param ploidy_level ploidy level of the offspring
#'
#' @return list of 2 matrices : matrix of mendelian transmission probabilities and matrix of mismatches
#'
#' @useDynLib APIS, .registration = TRUE
#'
#' @keywords internal
#' @noRd

get_individual_probability_2n <- function(offspring_genotype, sire_genotype, dam_genotype,
                                          allele_frequencies,
                                          ploidy_level) {
  # DESCRIPTION
  # Get the mendelian transmission probabilities and mismatches for all the parent pairs for one offspring
  #
  # INPUTS
  # offspring_genotype : genotype matrix of one offspring after recoding
  # sire_genotype : genotype matrix of the sires after recoding
  # dam_genotype : genotype matrix of the dams after recoding
  # allele_frequencies : allele frequency matrix after recoding
  # ploidy_level : ploidy level of the offspring
  #
  # OUTPUTS
  # list of 2 elements : mendelian transmission probabilities and mismatches for all the parental pairs

  # Prepare for Call
  offspring_genotypes_for_C = as.vector(t(offspring_genotype))
  sire_genotypes_for_C      = as.vector(t(sire_genotype))
  dam_genotypes_for_C       = as.vector(t(dam_genotype))
  allele_frequency_for_C    = as.vector(t(allele_frequencies))

  number_sire = nrow(sire_genotype)
  number_dam = nrow(dam_genotype)
  number_marker = ncol(offspring_genotype) / ploidy_level
  number_variant = ncol(allele_frequencies)

  all_probabilities = vector(mode = "numeric", length = number_sire * number_dam)
  all_mismatches    = vector(mode = "integer", length = number_sire * number_dam)

  output_C = .C('get_individual_mendelian_probability_2n',
                as.integer(offspring_genotypes_for_C),
                as.integer(sire_genotypes_for_C),
                as.integer(dam_genotypes_for_C),
                as.integer(number_sire),
                as.integer(number_dam),
                as.integer(number_marker),
                as.integer(number_variant),
                as.double(allele_frequency_for_C),
                as.double(all_probabilities),
                as.integer(all_mismatches))

  return(list(output_C[[9]], output_C[[10]]))
}

#' Estimate probabilities
#'
#' @param offspring_genotype offspring genotype matrix
#' @param sire_genotype sire genotype matrix
#' @param dam_genotype dan genotype matrix
#' @param allele_frequencies allele frequencies
#' @param t_recom recombinaison rate
#'
#' @return list of 2 matrices : matrix of mendelian transmission probabilities and matrix of mismatches
#'
#' @useDynLib APIS, .registration = TRUE
#'
#' @keywords internal
#' @noRd

get_individual_probability_3n <- function(offspring_genotype, sire_genotype, dam_genotype,
                                          allele_frequencies,
                                          ploidy_level,t_recom) {
  # DESCRIPTION
  # get probability matrix and mismatch matrix using C routine for one triploid individual
  #
  # INPUTS
  # offspring_genotype : genotype matrix of one offspring after recoding
  # sire_genotype : genotype matrix of the sires after recoding
  # dam_genotype : genotype matrix of the dams after recoding
  # allele_frequencies : allele frequencies after recoding
  # ploidy_level : ploidy level of the offspring

  # OUPUTS
  # list of 2 elements : matrix of probabilities and matrix of mismatches

  # Prepare for Call
  offspring_genotypes_for_C = as.vector(t(offspring_genotype))
  sire_genotypes_for_C      = as.vector(t(sire_genotype))
  dam_genotypes_for_C       = as.vector(t(dam_genotype))
  allele_frequency_for_C    = as.vector(t(allele_frequencies))

  number_sire = nrow(sire_genotype)
  number_dam = nrow(dam_genotype)
  number_marker = ncol(offspring_genotype) / ploidy_level
  number_variant = ncol(allele_frequencies)

  all_probabilities = vector(mode = "numeric", length = number_sire * number_dam)
  all_mismatches    = vector(mode = "integer", length = number_sire * number_dam)

  output_C = .C('get_individual_mendelian_probability_3n',
                as.integer(offspring_genotypes_for_C),
                as.integer(sire_genotypes_for_C),
                as.integer(dam_genotypes_for_C),
                as.integer(number_sire),
                as.integer(number_dam),
                as.integer(number_marker),
                as.integer(number_variant),
                as.double(allele_frequency_for_C),
                as.double(all_probabilities),
                as.integer(all_mismatches),
                as.double(t_recom)) # ajout JR -- 23/02/2023

  return(list(output_C[[9]], output_C[[10]]))
}

