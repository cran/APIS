# THIS FILE CONTAINS ALL THE FUNCTIONS FOR IMPORTING INPUTS
#
# =====================================================================================

# =================================================================================
#
# IMPORT FROM OTHER SOFTWARE FUNCTIONS
#

#' Import from Plink .ped
#'
#' @param ped_file name of the ped file (from Plink)
#' @param no_fid if "no_fid" parameter was used in plink (default : FALSE)
#' @param no_parents if "no_parents" parameter was used in plink (default : FALSE)
#' @param no_sex if "no_sex" parameter was used in plink (default : FALSE)
#' @param no_pheno if "no_pheno" parameter was used in plink (default : FALSE)
#' @param marker_names list of marker names (default : NULL)
#'
#' @return matrix of genotypes for APIS
#'
#' @importFrom data.table fread
#'
#' @export

import_from_ped <- function(ped_file,
                            no_fid = FALSE,
                            no_parents = FALSE,
                            no_sex = FALSE,
                            no_pheno = FALSE,
                            marker_names = NULL) {
  # DESCRIPTION
  # import .ped file in the APIS format
  #
  # INPUTS
  # ped_file : path to the .ped file
  # no_fid : "no_fid" argument from plink
  # no_parents : "no_parents" argument from plink
  # no_sex : "no_sex" argument from plink
  # no_pheno : "no_pheno" argument from plink
  # marker_names : names of the markers
  #
  # OUPUTS
  # ped_file : genotype matrix in the APIS format

  ped_file = as.data.frame(data.table::fread(ped_file, sep = "\t"))
  plink_param = c(!no_parents, !no_parents, !no_sex, !no_pheno)

  if (!no_fid) {
    ped_file = ped_file[, -1]
  }
  rownames(ped_file) = ped_file[, 1] ; ped_file = ped_file[, -1]
  col_remove = sum(plink_param)

  if (col_remove > 0) {
    ped_file = ped_file[, -c(1:col_remove)]
  }

  if (is.null(marker_names)) {
    colnames(ped_file) = paste0('marker_', c(1:ncol(ped_file)))
  } else {
    colnames(ped_file) = marker_names
  }

  ped_file = as.matrix(ped_file)
  storage.mode(ped_file) = "character"

  ped_file = gsub(x = ped_file, pattern = " ", replacement = "/")
  ped_file[which(ped_file == "0/0")] = "NA/NA"

  return(ped_file)
}

#' Import from .vcf
#'
#' @param vcf_file name of the vcf file
#'
#' @return matrix of genotypes for APIS
#'
#' @importFrom data.table fread
#' @export

import_from_vcf <- function(vcf_file) {
  # DESCRIPTION
  # import .vcf file in the APIS format
  #
  # INPUTS
  # ped_file : path to the .vcf file
  #
  # OUPUTS
  # apis_file : genotype matrix in the APIS format

  vcf_file = as.data.frame(data.table::fread(vcf_file))
  vcf_file = vcf_file[, -c(1:2, 4:9)]
  rownames(vcf_file) = vcf_file[, 1] ; vcf_file = vcf_file[, -1]
  apis_file <- t(vcf_file)

  if (!is(apis_file, 'matrix')) {
    apis_file <- as.matrix(apis_file)
  }

  storage.mode(apis_file) <- "character"

  apis_file = gsub(x = apis_file, pattern = "\\|", replacement = "/")
  apis_file = apply(apis_file, 2, function(X) {substr(X, start = 1, stop = 3)})

  apis_file[which(apis_file == "./.")] <- "NA/NA"

  return(apis_file)
}
