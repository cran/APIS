# THIS FILE CONTAINS ALL THE FUNCTIONS FOR ASSIGNMENT POWER
#
# =====================================================================================

#' Assignment power
#'
#' @param sire_genotype matrix of the sire genotypes
#' @param dam_genotype matrix of the dam genotypes
#' @param ploidy_level ploidy level of the parents
#' @param verbose verbose
#'
#' @return the theoretical assignment power calculated with the formula proposed in Vandeputte (2012)
#'
#' @examples
#' data("APIS_sire")
#' data("APIS_dam")
#'
#' P = assignment_power(sire_genotype = APIS_sire, dam_genotype = APIS_dam)
#'
#' @export

assignment_power <- function(sire_genotype, dam_genotype,
                             ploidy_level = 2,
                             verbose = FALSE) {
  # DESCRIPTION
  # This function calculates the theoretical assignment power as proposed in Vandeputte, M (2012)
  #
  # INPUTS
  # sire_genotype : genotype matrix of the sires in the APIS format
  # dam_genotype : genotype matrix of the dams in the APIS format
  # ploidy_level : ploidy level of the parental population
  # verbose : display information on the screen
  #
  # OUTPUTS
  # P : the theoretical assignment power

  parent_population <- rbind(sire_genotype, dam_genotype)

  # Importe the allFreq function and calculate the allele frequencies
  allele_frequency <- get_allele_frequencies(as.matrix(parent_population),
                                             ploidy_level = ploidy_level)

  column <- which(colnames(allele_frequency) == 'tot')
  calculate_frequency <- as.data.frame(allele_frequency[,((column+1):ncol(allele_frequency))])

  test_NA <- which(colnames(allele_frequency) == "Freq_NA")

  if (length(test_NA) != 0) {
    calculate_frequency <- calculate_frequency[,-test_NA]
  }

  total_column <- ncol(calculate_frequency)

  # Calculate Q1 and Q3 for each marker
  calculate_frequency$Q1i <- 1 - 2*rowSums(calculate_frequency[,1:total_column]^2) +
    rowSums(calculate_frequency[,1:total_column]^3) + 2*rowSums(calculate_frequency[,1:total_column]^4) -
    2*rowSums(calculate_frequency[,1:total_column]^2)^2 - 3*rowSums(calculate_frequency[,1:total_column]^5) +
    3*rowSums(calculate_frequency[,1:total_column]^3)*rowSums(calculate_frequency[,1:total_column]^2)

  calculate_frequency$Q3i <- 1 + 4*rowSums(calculate_frequency[,1:total_column]^4) -
    4*rowSums(calculate_frequency[,1:total_column]^5) - 3*rowSums(calculate_frequency[,1:total_column]^6) -
    8*rowSums(calculate_frequency[,1:total_column]^2)^2 + 2*rowSums(calculate_frequency[,1:total_column]^3)^2 +
    8*rowSums(calculate_frequency[,1:total_column]^3)*rowSums(calculate_frequency[,1:total_column]^2)

  # Calculate the global Q1 and Q3
  Q1 <- 1 - prod(1-calculate_frequency$Q1i)
  Q3 <- 1 - prod(1-calculate_frequency$Q3i)

  # Calculate the assignment power
  P <- Q1^(nrow(dam_genotype)+nrow(sire_genotype)-2)*Q3^((nrow(dam_genotype)-1)*(nrow(sire_genotype)-1))

  if (verbose) {
    P_char = substr(P, 1, 5)
    cat(paste0("theoretical assignment power : ", P_char))
    cat("\n")
    cat('---------------------------------------------------', sep = '\n')
  }

  # Return the result
  return(P)
}
