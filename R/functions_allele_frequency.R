# THIS FILE CONTAINS ALL THE FUNCTIONS FOR ESTIMATING ALLELE FREQUENCIES
#
# =====================================================================================

#' Estimate allele frequencies
#'
#' @param genotype_matrix matrix of genotypes
#' @param ploidy_level level of ploidy
#'
#' @return matrix of allele frequencies
#'
#' @examples
#' data("APIS_offspring")
#' allele_frequencies = get_allele_frequencies(genotype_matrix = APIS_offspring,
#'                                             ploidy_level = 2)
#'
#' @keywords internal
#' @noRd
#'

get_allele_frequencies <- function(genotype_matrix, ploidy_level = 2) {
  # DESCRIPTION
  # Estimate allele frequencies based on genotype matrix
  #
  # INPUTS
  # genotype_matrix : genotype matrix in the APIS format
  # ploidy_level : ploidy level of the individuals of the genotype_matrix
  #
  # OUTPUTS
  # result : dataframe with the count of each allele and their frequency

  genotype_matrix_two_columns <- matrix(NA, nrow = nrow(genotype_matrix), ncol = ploidy_level*ncol(genotype_matrix))

  import_column <- seq(1,ncol(genotype_matrix_two_columns),ploidy_level)

  # Divide each genotype (coded A/A) into 2 columns
  shift <- ploidy_level - 1
  for (i in c(1:ncol(genotype_matrix))) {
    tmp <- strsplit(genotype_matrix[,i], split = '/', fixed = TRUE)
    M <- t(mapply(FUN = function(X) {X}, tmp))
    genotype_matrix_two_columns[,(import_column[i]:(import_column[i]+shift))] <- M
  }

  # List of the different alleles
  variant <- sort(unique(unlist(as.list(apply(genotype_matrix_two_columns, 2, unique)))))

  # Create the results matrix
  result_matrix <- matrix(0, nrow = ncol(genotype_matrix), ncol = length(variant))
  rownames(result_matrix) <- colnames(genotype_matrix)
  colnames(result_matrix) <- variant

  for (n in 1:nrow(result_matrix)) {
    tmp <- table(genotype_matrix_two_columns[,(import_column[n]:(import_column[n]+shift))])
    result_matrix[n,match(names(tmp), colnames(result_matrix))] <- tmp
  }

  # Estimate the allele frequencies
  frequency_matrix <- result_matrix / (rowSums(result_matrix[,which(colnames(result_matrix) != 'NA')]))
  colnames(frequency_matrix) <- paste0('Freq_', colnames(result_matrix))

  # Merge the results
  result <- cbind(result_matrix, tot = rowSums(result_matrix), frequency_matrix)

  # Return the result
  return(result)
}

#' Create variant dataframe
#'
#' @param offspring_genotype matrix of the offspring genotypes
#' @param sire_genotype matrix of the sire genotypes
#' @param dam_genotype matrix of the dam genotypes
#'
#' @return dataframe with the correspondance between input allele and recoded allele
#'
#' @keywords internal
#' @noRd

get_variant_dataframe <- function(offspring_genotype, sire_genotype, dam_genotype) {
  # DESCRIPTION
  # Get all the alleles for all the markers in the population
  #
  # INPUTS
  # offspring_genotype : genotype matrix of the offspring in the APIS format
  # sire_genotype : genotype matrix of the sires in the APIS format
  # dam_genotype : genotype matrix of the dams in the APIS format
  #
  # OUTPUTS
  # variant_dataframe : dataframe with all possible alleles for all the markers and their new code

  variant <- unique(unlist(strsplit(as.vector(rbind(offspring_genotype, sire_genotype, dam_genotype)), '/')))
  NA_variant <- which(variant == "NA")
  if (length(NA_variant) > 0) {
    variant <- variant[-which(variant == "NA")]
  }

  variant_dataframe <- data.frame(variant = as.character(variant),
                                  recode = c(1:length(variant)))
  variant_dataframe$variant <- as.character(variant_dataframe$variant)
  variant_dataframe <- rbind(variant_dataframe, c(as.character("NA"), 0))

  return(variant_dataframe)
}
