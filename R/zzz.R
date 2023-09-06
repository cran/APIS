.onAttach <- function(libname, pkgname) {
  package_citation <- "Griot et al. (2020). APIS: An auto-adaptive parentage inference software that tolerates missing parents"
  doi = "https://doi.org/10.1111/1755-0998.13103"
  # message("Thank you for using APIS!")
  # message("To acknowledge our work, please cite the package:")
  # message(package_citation)
  # message(doi)
  packageStartupMessage("Thank you for using APIS!")
  packageStartupMessage("To acknowledge our work, please cite the package:")
  packageStartupMessage(package_citation)
  packageStartupMessage(doi)
}
