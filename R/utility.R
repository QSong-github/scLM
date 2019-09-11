
#' @title Poisson vector.
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer.
#' @param same row of the data in dataList trydata
#' @return A numeric vector of log density.
#' @export main function
InstalledPackage <- function(package)
{
  available <- suppressMessages(suppressWarnings(sapply(package, require, quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
  missing <- package[!available]
  if (length(missing) > 0) return(FALSE)
  return(TRUE)
}


#' @title Poisson vector.
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer.
#' @param same row of the data in dataList trydata
#' @return A numeric vector of log density.
#' @export main function
CRANChoosen <- function()
{
  return(getOption("repos")["CRAN"] != "@CRAN@")
}


#' check the avaliablity of neccessary packages
#' @title Poisson vector.
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer.
#' @param same row of the data in dataList trydata
#' @return A numeric vector of log density.
#' @export main function
UsePackage <- function(package, defaultCRANmirror = "http://cran.at.r-project.org") {
  if(!InstalledPackage(package))
  {
    if(!CRANChoosen())
    {
      chooseCRANmirror(graphics=FALSE,ind=87)
      if(!CRANChoosen())
      {
        options(repos = c(CRAN = defaultCRANmirror))
      }
    }

    suppressMessages(suppressWarnings(install.packages(package)))
    if(!InstalledPackage(package)) return(FALSE)
  }
  return(TRUE)
}


for(library in c('MASS','lattice','zic','rpart'))
{
  if(!UsePackage(library))
  {
    stop("Error!", library)
  }
}

if(!InstalledPackage('mpath')){
  install.packages('https://cran.r-project.org/src/contrib/mpath_0.3-17.tar.gz',repos=NULL)}

if(!InstalledPackage('Matrix')){
  install.packages('https://cran.r-project.org/src/contrib/Matrix_1.2-17.tar.gz',repos=NULL)}


for( p in c('pscl','glmnet','parallel','purrr','Matrix')){
  if(!require(p,character.only = TRUE)) install.packages(p, repos='http://cran.us.r-project.org')
  library(p,character.only = TRUE)
}

