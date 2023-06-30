# THIS FILE CONTAINS ALL THE FUNCTIONS FOR DIPLOID ASSIGNMENT USING PARALLEL PROGRAMMING
#
# =====================================================================================

#' Get probabilities
#'
#' @param offspring_genotype matrix of the offspring genotypes recoded (after recode_genotypes_for_assignment() function)
#' @param sire_genotype matrix of the sire genotypes recoded (after recode_genotypes_for_assignment() function)
#' @param dam_genotype matrix of the dam genotypes recoded (after recode_genotypes_for_assignment() function)
#' @param allele_frequencies allele frequencies recoded (after recode_allele_frequencies_for_assignment() function)
#' @param method method of assignment ("mendel" or "exclusion")
#' @param ploidy_level ploidy level
#' @param number_cores number of cores
#'
#' @return list of probabilities and mismatch
#'
#' @useDynLib APIS, .registration = TRUE
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#'
#' @keywords internal
#' @noRd

get_probabilities_parallel_2n <- function(offspring_genotype, sire_genotype, dam_genotype,
                                          allele_frequencies,
                                          method = "mendel",
                                          ploidy_level,
                                          number_cores) {
  # DESCRIPTION
  # Get the mendelian transmission probabilities and mismatches for all the parent pairs
  # for all the diploid offspring individuals using parallel programming
  #
  # INPUTS
  # offspring_genotype : genotype matrix of the offspring after recoding
  # sire_genotype : genotype matrix of the sires after recoding
  # dam_genotype : genotype matrix of the dams after recoding
  # allele_frequencies : allele frequency matrix after recoding
  # method : method of assignment
  # ploidy_level : ploidy level of the offspring
  # number_cores : number of cores for parallel processing
  #
  # OUTPUTS
  # output : the log file

  i = NULL
  iterations = nrow(offspring_genotype)
  new_genotype = NULL

  cl = parallel::makeCluster(number_cores)
  doParallel::registerDoParallel(cl)

  output <- foreach(i=1:iterations, .combine=rbind, .export="get_individual_probability_2n") %dopar% {
    tmp <- get_individual_probability_2n(offspring_genotype = offspring_genotype[i, , drop = FALSE],
                                         sire_genotype = sire_genotype,
                                         dam_genotype = dam_genotype,
                                         allele_frequencies = allele_frequencies,
                                         ploidy_level = ploidy_level)
    all_values = data.frame(sire = rep(rownames(sire_genotype), each = nrow(dam_genotype)),
                            dam = rep(rownames(dam_genotype), times = nrow(sire_genotype)),
                            probability = tmp[[1]],
                            mismatch = tmp[[2]])
    if (method == "mendel" | method == "exclusion") {
      if (method == "mendel"){
        all_values = all_values[order(all_values$probability, decreasing = TRUE), ]
      } else {
        all_values = all_values[order(all_values$mismatch, decreasing = FALSE), ]
      }
      log_file = c(rownames(offspring_genotype)[i],
                   all_values$sire[1], all_values$dam[1], all_values$probability[1], all_values$mismatch[1],
                   all_values$sire[2], all_values$dam[2], all_values$probability[2], all_values$mismatch[2],
                   all_values$sire[3], all_values$dam[3], all_values$probability[3], all_values$mismatch[3],
                   all_values$probability[1] - all_values$probability[2],
                   all_values$probability[2] - all_values$probability[3])
    } else {
      all_values1 = all_values[order(all_values$probability, decreasing = TRUE), ]
      all_values2 = all_values[order(all_values$mismatch, decreasing = FALSE), ]

      log_file = c(rownames(offspring_genotype)[i],
                   all_values1$sire[1], all_values1$dam[1], all_values1$probability[1], all_values1$mismatch[1],
                   all_values1$sire[2], all_values1$dam[2], all_values1$probability[2], all_values1$mismatch[2],
                   all_values1$sire[3], all_values1$dam[3], all_values1$probability[3], all_values1$mismatch[3],
                   all_values1$probability[1] - all_values1$probability[2],
                   all_values1$probability[2] - all_values1$probability[3],
                   all_values2$sire[1], all_values2$dam[1], all_values2$probability[1], all_values2$mismatch[1],
                   all_values2$sire[2], all_values2$dam[2], all_values2$probability[2], all_values2$mismatch[2],
                   all_values2$sire[3], all_values2$dam[3], all_values2$probability[3], all_values2$mismatch[3],
                   all_values2$probability[1] - all_values2$probability[2],
                   all_values2$probability[2] - all_values2$probability[3])

    }

    return(log_file)
  }
  parallel::stopCluster(cl)

  return(output)
}

# =====================================================================================

#' APIS for diploids
#'
#' @param offspring_genotype matrix of the offspring genotypes
#' @param sire_genotype matrix of the sire genotypes
#' @param dam_genotype matrix of the offspring genotypes
#' @param method method : "mendel" i.e. likelihood or "exclusion" (default : "mendel"). Can also be "" to select the method a posteriori.
#' @param exclusion_threshold threshold for "exclusion" method (default : NULL). Override the error parameter if not NULL
#' @param error error accepted (default : 0.05)
#' @param simulation_if_small simulate individuals (TRUE or FALSE)
#' @param number_offspring_simulated number of offspring simulated (default : 500)
#' @param number_cores number of cores
#' @param verbose verbose
#'
#' @return list of 2 elements : a pedigree file and the log file
#'
#' @importFrom cowplot plot_grid
#' @useDynLib APIS, .registration = TRUE
#'
#' @examples
#' data("APIS_offspring")
#' data("APIS_sire")
#' data("APIS_dam")
#'
#' assignment <- APIS_2n(offspring_genotype = APIS_offspring[1:35,1:50],
#'                       sire_genotype = APIS_sire[ ,1:50],
#'                       dam_genotype = APIS_dam[ ,1:50],
#'                       simulation_if_small = FALSE)
#'
#' @export

APIS_2n <- function(offspring_genotype, sire_genotype, dam_genotype,
                    method = "mendel",
                    exclusion_threshold = NULL, error = 0.05,
                    simulation_if_small = FALSE,
                    number_offspring_simulated = max(0, 500 - nrow(offspring_genotype)),
                    number_cores = 2,
                    verbose = FALSE) {
  # DESCRIPTION
  # APIS process for diploid offspring
  #
  # INPUTS
  # offspring_genotype : genotype matrix of one offspring in APIS format
  # sire_genotype : genotype matrix of the sires in APIS format
  # dam_genotype : genotype matrix of the dams in APIS format
  # method : method of assignment
  # exclusion_threshold : number of mismatches allowed in "exclusion" method
  # error : percent of error allowed in "mendel" method
  # simulation_if_small : simulate individuals
  # number_offspring_simulated : number of offspring to simulate
  # number_cores : number of cores for parallel programming
  # verbose : display information on the screen
  #
  # OUTPUTS
  # list of 2 elements : the pedigree and the log file

  if (verbose) {
    print_title()
  }

  # Parental ploidy level
  parental_ploidy_level = length(strsplit(sire_genotype[1, 1], split = "/")[[1]])

  # Check inputs
  check_matrices_format(offspring_genotype, verbose = verbose)
  check_matrices_format(sire_genotype, verbose = verbose)
  check_matrices_format(dam_genotype, verbose = verbose)

  genotypes_checked   = check_genotypes(offspring_genotype, sire_genotype, dam_genotype,
                                        offspring_ploidy_level = 2,
                                        parental_ploidy_level = parental_ploidy_level,
                                        verbose = verbose)
  offspring_genotype  = genotypes_checked[[1]]
  sire_genotype       = genotypes_checked[[2]]
  dam_genotype        = genotypes_checked[[3]]

  if (!is(offspring_genotype, "matrix")) {
    offspring_genotype = t(as.matrix(offspring_genotype))
  }

  check_input_parameters(exclusion_threshold = exclusion_threshold,
                         number_marker = ncol(offspring_genotype),
                         error = error,
                         verbose = verbose)

  # Theoretical assignment power
  P = assignment_power(sire_genotype = sire_genotype,
                       dam_genotype = dam_genotype,
                       ploidy_level = parental_ploidy_level,
                       verbose = verbose)

  # Add simulation
  if (simulation_if_small & number_offspring_simulated > 0) {
    simulation = simulate_offspring(sire_genotype = sire_genotype,
                                    dam_genotype = dam_genotype,
                                    number_offspring = number_offspring_simulated)
    offspring_genotype = rbind(offspring_genotype, simulation$genotypes)
  }

  # Recode genotype matrices
  prepared_data = prepare_for_assignment(offspring_genotype = offspring_genotype,
                                         sire_genotype = sire_genotype,
                                         dam_genotype = dam_genotype,
                                         ploidy_level = 2,
                                         number_cores = number_cores,
                                         verbose = verbose)
  offspring_genotype_recoded  = prepared_data[[1]]
  sire_genotype_recoded       = prepared_data[[2]]
  dam_genotype_recoded        = prepared_data[[3]]
  allele_frequency_recoded    = prepared_data[[4]]

  # Calculate probabilities and mismatches
  if (verbose) {
    cat("estimation of mendelian transmission probabilities and mismatches", sep = "\n")
  }

  log_file = get_probabilities_parallel_2n(offspring_genotype = offspring_genotype_recoded,
                                           sire_genotype = sire_genotype_recoded,
                                           dam_genotype = dam_genotype_recoded,
                                           allele_frequencies = allele_frequency_recoded,
                                           method = method,
                                           ploidy_level = 2,
                                           number_cores = number_cores)
  if (is(log_file, "vector")) {
    log_file = as.data.frame(t(log_file))
  } else {
    log_file = as.data.frame(log_file)
  }

  if (method %in% c("mendel","exclusion")){
    log_file = format_logfile(log_file)
  } else {
    log_file_likelihood = log_file[,1:15]
    log_file_likelihood = format_logfile(log_file_likelihood)

    log_file_exclusion = log_file[,c(1,16:29)]
    log_file_exclusion = format_logfile(log_file_exclusion)
  }

  if (verbose) {
    cat("assignment ...", sep = "\n")
  }

  # if method = mendel
  if(method == "mendel") {
    # Create pedigree
    threshold_mendel = estimate_mendel_threshold(log_file = log_file,
                                                 error = error)
    pedigree = create_pedigree_delta(log_file = log_file,
                                     threshold_delta = threshold_mendel)

    # Plot graphs
    if (simulation_if_small & number_offspring_simulated > 0) {
      simulated_individuals = simulation$pedigree[, 1]
    } else {
      simulated_individuals = NULL
    }

    p1 = plot_delta(log_file = log_file, threshold = threshold_mendel,
                    simulated_individuals = simulated_individuals)
    p2 = plot_probabilities(log_file = log_file, threshold = NULL,
                            simulated_individuals = simulated_individuals)
    p3 = plot_mismatches(log_file = log_file, threshold = NULL,
                         simulated_individuals = simulated_individuals)

    # if method = exclusion
  } else if (method == "exclusion") {
    # Create pedigree
    if (is.null(exclusion_threshold)){
      exclusion_threshold = estimate_exclusion_threshold(log_file = log_file,
                                                         error = error)
    } else if (verbose){
      cat("error override as exclusion_threshold provided", sep = "\n")
    }

    pedigree = create_pedigree_exclusion(log_file = log_file,
                                         threshold_exclusion = exclusion_threshold)

    # Plot graphs
    if (simulation_if_small & number_offspring_simulated > 0) {
      simulated_individuals = simulation$pedigree[, 1]
    } else {
      simulated_individuals = NULL
    }

    p1 = plot_delta(log_file = log_file, threshold = NULL,
                    simulated_individuals = simulated_individuals)
    p2 = plot_probabilities(log_file = log_file, threshold = NULL,
                            simulated_individuals = simulated_individuals)
    p3 = plot_mismatches(log_file = log_file, threshold = exclusion_threshold,
                         simulated_individuals = simulated_individuals)
  } else { # method for selection a posteriori for the shiny app
    threshold_mendel = estimate_mendel_threshold(log_file = log_file_likelihood,
                                                 error = error)
    pedigree_likelihood = create_pedigree_delta(log_file = log_file_likelihood,
                                                threshold_delta = threshold_mendel)

    if (is.null(exclusion_threshold)){
      exclusion_threshold = estimate_exclusion_threshold(log_file = log_file_exclusion,
                                                         error = error)
    } else if (verbose){
      cat("error override as exclusion_threshold provided", sep = "\n")
    }

    pedigree_exclusion = create_pedigree_exclusion(log_file = log_file_exclusion,
                                                   threshold_exclusion = exclusion_threshold)
  }

  if (verbose) {
    cat("assignment completed !", sep = "\n")
    if (method %in% c("mendel","exclusion")){
      print_summary(theoretical_assignment_power = P,
                    assignment_rate = get_assignment_rate(pedigree))
    }
  }

  if (method %in% c("mendel","exclusion")) {
    print(cowplot::plot_grid(p1, p2, p3, nrow = 3))
    if (simulation_if_small & number_offspring_simulated > 0) {
      output = remove_simulation(log_file = log_file,
                                 pedigree = pedigree,
                                 simulated_individuals = simulated_individuals)
      return(output)
    } else {
      return(list(pedigree = pedigree, log_file = log_file))
    }
  } else {
    if (simulation_if_small & number_offspring_simulated > 0) {
      output1 = remove_simulation(log_file = log_file_likelihood,
                                  pedigree = pedigree_likelihood,
                                  simulated_individuals = simulated_individuals)
      output2 = remove_simulation(log_file = log_file_exclusion,
                                  pedigree = pedigree_exclusion,
                                  simulated_individuals = simulated_individuals)
      return(c(output1,output2))
    } else {
      return(list(pedigree_likelihood=pedigree_likelihood,pedigree_exclusion=pedigree_exclusion,
                  log_file_likelihood=log_file_likelihood,log_file_exclusion=log_file_exclusion))
    }
  }
}
