.onAttach <- function(libname, pkgname) {
  package_citation <- "Griot et al. (2020). APIS: An auto-adaptive parentage inference software that tolerates missing parents"
  doi = "https://doi.org/10.1111/1755-0998.13103"
  packageStartupMessage("Thank you for using APIS!")
  packageStartupMessage("To acknowledge our work, please cite the package:")
  packageStartupMessage(package_citation)
  packageStartupMessage(doi)
  packageStartupMessage("And for the triploid version:")
  packageStartupMessage("Roche et al. (2024). APIS: an updated parentage assignment software managing triploids induced from diploid parents")
  packageStartupMessage("https://doi.org/10.1093/g3journal/jkae143")
}
