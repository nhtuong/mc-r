
#'@name check.installed.mc
#'@aliases check.installed.mc
#'@export check.installed.mc
#'@docType methods
#'@title Checking package installation
#'@description Check whether a packages is installed
#'@param pkg name of package
#'@return A logical value indicating whether the package is installed
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
check.installed.mc<-function(pkg){
  return(is.element(pkg, installed.packages()[,1]))
}



#'@name library.mc
#'@aliases library.mc
#'@export library.mc
#'@docType methods
#'@title Loading and Listing of Packages
#'@description On-the-fly load or install a package
#'@param pkg name of package
#'@return A list of attached packages
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
library.mc<-function(pkg,repos="cran"){
  if(!check.installed.mc(pkg)){
    
    if (repos=="bioc"){
      source("http://www.bioconductor.org/biocLite.R")
      biocLite(pkg)
    } else if (repos=="cran"){
      install.packages(pkg)
    } else {
      install.packages(pkg, repos = repos, type="source")
    }
    
    
  }
  
  library(pkg,character.only=TRUE)
  
  
}
