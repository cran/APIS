#' Simulate offspring
#'
#' @param sire_genotype sire genotype
#' @param dam_genotype dam genotype
#' @param number_offspring number of offspring to simulate
#' @param ploidy_level ploidy level of offspring
#' @param sire_contribution sire contribution
#' @param dam_contribution dam contribution
#' @param recombination_rate recombination rate (only important for tri/tetra ploids offspring)
#' @param genotyping_error genotyping error
#'
#' @return list with matrix with simulated offspring and pedigree
#'
#' @examples
#' data("APIS_sire")
#' data("APIS_dam")
#'
#' # For diploide offspring
#' simulate_offspring(sire_genotype=APIS_sire, dam_genotype=APIS_dam,
#'                    number_offspring=10,
#'                    ploidy_level = 2,
#'                    sire_contribution = 1, dam_contribution = 1,
#'                    recombination_rate = 0.5,
#'                    genotyping_error = 0.01)
#'
#' # For triploide offspring
#' simulate_offspring(sire_genotype=APIS_sire, dam_genotype=APIS_dam,
#'                    number_offspring=10,
#'                    ploidy_level = 3,
#'                    sire_contribution = 1, dam_contribution = 2,
#'                    recombination_rate = 0.5,
#'                    genotyping_error = 0.01)
#'
#' @export

simulate_offspring = function(sire_genotype, dam_genotype,
                               number_offspring,
                               ploidy_level = 2,
                               sire_contribution = 1, dam_contribution = 1,
                               recombination_rate = 0.5,
                               genotyping_error = 0.01) {

  if (sire_contribution + dam_contribution != ploidy_level) {
    stop('sire_contribution + dam_contribution must be equal to the ploidy level')
  }
  if (sire_contribution>2 | dam_contribution>2){
    stop('sire_contribution and dam_contribution must be equal inferior to 2. Method not implemented above.')
  }
  n_sire = nrow(sire_genotype)
  n_dam = nrow(dam_genotype)
  n_marker = ncol(sire_genotype)

  mat_off = matrix(nrow=number_offspring,ncol = n_marker)
  colnames(mat_off) = colnames(sire_genotype)
  rownames(mat_off) = paste0("Simulated_off_",1:number_offspring)

  simulated_pedigree = matrix(NA, nrow = number_offspring, ncol = 3)
  colnames(simulated_pedigree) = c('individual', 'sire', 'dam')
  simulated_pedigree[,1] = rownames(mat_off)

  all_allele = list()
  for (a in 1:n_marker){
    all_allele[[a]] = unique(c(sire_genotype[,a],dam_genotype[,a]))
    all_allele[[a]] = unique(unlist(strsplit(x = c(sire_genotype[,a],dam_genotype[,a]),split = "/",fixed = TRUE)))
    test_na = which(all_allele[[a]]=='NA')
    if (length(test_na)>0){
      all_allele[[a]] = all_allele[[a]][-test_na]
    }
  }

  id_selected_sire = sample(x = 1:n_sire,size = number_offspring,replace = T)
  id_selected_dam = sample(x = 1:n_dam,size = number_offspring,replace = T)

  simulated_pedigree[,2] = rownames(sire_genotype)[id_selected_sire]
  simulated_pedigree[,3] = rownames(dam_genotype)[id_selected_dam]

  for (k in 1:number_offspring){
    print(k)
    id_sire = id_selected_sire[k]
    id_dam = id_selected_dam[k]

    sire_selected = strsplit(x = sire_genotype[id_sire,],split = "/",fixed = TRUE)
    dam_selected = strsplit(x = dam_genotype[id_dam,],split = "/",fixed = TRUE)

    which_allele_sire = sample(x = 1:2,size = n_marker,replace = TRUE)
    which_allele_dam = sample(x = 1:2,size = n_marker,replace = TRUE)

    if (sire_contribution==2){
      which_allele_sire2 = c()
      for (l in 1:n_marker){
        if (which_allele_sire[l]==1){
          id2=sample(1:2,size=1,replace=F,prob=c(recombination_rate,1-recombination_rate))
        } else {
          id2=sample(1:2,size=1,replace=F,prob=c(1-recombination_rate,recombination_rate))
        }
        which_allele_sire2 = c(which_allele_sire2,id2)
      }
    }
    if (dam_contribution==2){
      which_allele_dam2 = c()
      for (l in 1:n_marker){
        if (which_allele_dam[l]==1){
          id2=sample(1:2,size=1,replace=F,prob=c(recombination_rate,1-recombination_rate))
        } else {
          id2=sample(1:2,size=1,replace=F,prob=c(1-recombination_rate,recombination_rate))
        }
        which_allele_dam2 = c(which_allele_dam2,id2)
      }
    }
    offspring_genotype = c()
    if (dam_contribution==1 & sire_contribution==1){
      for (n in 1:n_marker){
        sire_allele = mutate_allele(allele = sire_selected[[n]][which_allele_sire[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele = mutate_allele(allele = dam_selected[[n]][which_allele_dam[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        off_allele = paste0(sort(c(sire_allele, dam_allele)), collapse = "/")
        if (length(grep(pattern = "NA",x = off_allele))>0){
          off_allele="NA/NA"
        }
        offspring_genotype[n] = off_allele
      }
    } else if (dam_contribution==2 & sire_contribution==1){
      for (n in 1:n_marker){
        sire_allele = mutate_allele(allele = sire_selected[[n]][which_allele_sire[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele = mutate_allele(allele = dam_selected[[n]][which_allele_dam[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele2 = mutate_allele(allele = dam_selected[[n]][which_allele_dam2[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        off_allele = paste0(sort(c(sire_allele, dam_allele,dam_allele2)), collapse = "/")
        if (length(grep(pattern = "NA",x = off_allele))>0){
          off_allele="NA/NA/NA"
        }
        offspring_genotype[n] = off_allele
      }
    } else if(dam_contribution==1 & sire_contribution==2){
      for (n in 1:n_marker){
        sire_allele = mutate_allele(allele = sire_selected[[n]][which_allele_sire[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        sire_allele2 = mutate_allele(allele = sire_selected[[n]][which_allele_sire2[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele = mutate_allele(allele = dam_selected[[n]][which_allele_dam[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        off_allele = paste0(sort(c(sire_allele,sire_allele2, dam_allele)), collapse = "/")
        if (length(grep(pattern = "NA",x = off_allele))>0){
          off_allele="NA/NA/NA"
        }
        offspring_genotype[n] = off_allele
      }
    } else if(dam_contribution==2 & sire_contribution==2){
      for (n in 1:n_marker){
        sire_allele = mutate_allele(allele = sire_selected[[n]][which_allele_sire[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        sire_allele2 = mutate_allele(allele = sire_selected[[n]][which_allele_sire2[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele = mutate_allele(allele = dam_selected[[n]][which_allele_dam[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        dam_allele2 = mutate_allele(allele = dam_selected[[n]][which_allele_dam2[n]],possib_allele = all_allele[[n]],error_rate = genotyping_error)
        off_allele = paste0(sort(c(sire_allele,sire_allele2, dam_allele,dam_allele2)), collapse = "/")
        if (length(grep(pattern = "NA",x = off_allele))>0){
          off_allele="NA/NA/NA/NA"
        }
        offspring_genotype[n] = off_allele
      }
    }
    mat_off[k,]=offspring_genotype
  }
  return(list(genotypes = mat_off, pedigree = simulated_pedigree))
}

#' Simulate mutation
#'
#' @param allele parental allele
#' @param possib_allele possible allele for this marker
#' @param error_rate error rate
#'
#' @return mutated allele or not mutated allele
#'
#' @importFrom stats runif
#'
#' @keywords internal
#' @noRd

mutate_allele=function(allele,possib_allele,error_rate=0.1){
  if (length(possib_allele)>0){
    if (runif(1,0,1)<error_rate){
      new_allele = sample(x = possib_allele,size = 1)
      return(new_allele)
    } else {
      return(allele)
    }
  } else {
    return(allele)
  }
}

#' Remove simulated offspring from results
#'
#' @param log_file log_file
#' @param pedigree pedigree
#' @param simulated_individuals simulated offspring
#'
#' @return void
#'
#' @keywords internal
#' @noRd

remove_simulation <- function(log_file, pedigree, simulated_individuals) {
  # DESCRIPTION
  # remove simulated individuals from the APIS outputs
  #
  # INPUTS
  # log_file : the log file
  # pedigree : the pedigree
  # simulated_individuals : names of the simulated individuals
  #
  # OUTPUTS
  # list of 2 elements : the new log file and pedigree without the simulated individuals

  new_log = subset(log_file, !(log_file$offspring %in% simulated_individuals))
  new_pedigree = subset(pedigree, !(pedigree$offspring %in% simulated_individuals))

  return(list(log_file = new_log, pedigree = new_pedigree))
}
