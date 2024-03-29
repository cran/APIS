% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_APIS_3n.R
\name{APIS_3n}
\alias{APIS_3n}
\title{APIS for triploids}
\usage{
APIS_3n(
  offspring_genotype,
  sire_genotype,
  dam_genotype,
  method = "mendel",
  exclusion_threshold = NULL,
  error = 0.05,
  simulation_if_small = FALSE,
  number_offspring_simulated = max(0, 500 - nrow(offspring_genotype)),
  number_cores = 2,
  verbose = FALSE,
  t_recom = 0.5
)
}
\arguments{
\item{offspring_genotype}{matrix of the offspring genotypes}

\item{sire_genotype}{matrix of the sire genotypes}

\item{dam_genotype}{matrix of the dam genotypes}

\item{method}{method : "mendel" i.e. likelihood or "exclusion" (default : "mendel"). Can also be "" to select the method a posteriori.}

\item{exclusion_threshold}{threshold for "exclusion" method (default : NULL). Override the error parameter if not NULL}

\item{error}{error accepted (default : 0.05)}

\item{simulation_if_small}{simulate individuals (TRUE or FALSE) (default : TRUE)}

\item{number_offspring_simulated}{number of offspring simulated (default : 500)}

\item{number_cores}{number of cores}

\item{verbose}{verbose}

\item{t_recom}{recombination rate}
}
\value{
list of 2 elements : a pedigree file and the log file
}
\description{
APIS for triploids
}
\examples{
data("APIS_offspring3n")
data("APIS_sire")
data("APIS_dam")
assignment <- APIS_3n(offspring_genotype = APIS_offspring3n[1:35,1:50],
                      sire_genotype = APIS_sire[ ,1:50],
                      dam_genotype = APIS_dam[ ,1:50],
                      simulation_if_small = FALSE)

}
