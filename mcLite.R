#'@docType package
#'Package: mc
#'Type: Package
#'Title: Commonly used function for Nutriomique Team (INSERM U872)
#'Version: 0.1
#'Date: 2013-06-07
#'Author: Aurelie Cotillard, Edi Prifti, Hoai Tuong Nguyen
#'Maintainer: Hoai Tuong Nguyen <hoai-tuong.nguyen@inserm.fr>
#'Description: Statistical and datamining tools for metagenomic data analysis.
#'License: PPL


#'@docType methods
#'@title Checking package installation
#'@description Check whether a packages is installed
#'@param pkg name of package
#'@return A logical value indicating whether the package is installed
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
#'@aliases check.installed.mc
#'@rdname check.installed.mc
#'@export check.installed.mc
check.installed.mc<-function(pkg){
  return(is.element(pkg, installed.packages()[,1]))
}

#'@docType methods
#'@title Loading and Listing of Packages
#'@description On-the-fly load or install a package
#'@param pkg name of package
#'@return A list of attached packages
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
#'@aliases library.mc
#'@rdname library.mc
#'@export library.mc
library.mc<-function(pkg){
if(!check.installed.mc(pkg))
  install.packages(pkg)
  print(pkg)
  library(pkg,character.only=TRUE)
}

#Load dependencies
library.mc("xtable")


#'@docType methods
#'@title Data Input
#'@description Read a very large data file
#'@param file the name of the file which the data are to be read from.
#'@param header a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.
#'@param sep  the field separator character. Values on each line of the file are separated by this character. If sep = "" (the default for read.table) the separator is 'white space', that is one or more spaces, tabs, newlines or carriage returns. A (character) name of the column in a \code{\link[base]{data.frame}} which contains the addresses
#'@param nrow number of rows
#'@return A data frame (\code{\link[base]{data.frame}}) containing a representation of the data in the file
#'@author Hoai Tuong Nguyen
#'@examples
#'library(mc)
#'data<-read.table.mc("http://statistics.vn/data/doesgenes.txt",header=T,sep=";",nrow=1000)
#'@seealso \code{\link[utils]{read.table}}
#'@aliases read.table.mc
#'@rdname read.table.mc
#'@export read.table.mc
read.table.mc<-function(file,header=FALSE,sep="",nrow=-1){
  tab5rows <- read.table(file, nrows = 5,sep=sep)
  classes <- sapply(tab5rows, class)
  tabAll <- read.table(file,  header=header, colClasses=classes,sep=sep,nrows=nrow,comment.char = "")
  return(tabAll)
}





#'@docType methods
#'@title Object Summaries
#'@description Summarize an numeric table, save the output to a table, export the output to Latex code
#'@param object an object for which a summary is desired.
#'@param latex a logical value indicating whether output to latex is called
#'@return Table of of the value returned by summary (and output Latex code if parameter 'latex' is TRUE)
#'@author Hoai Tuong Nguyen
#'@examples
#'data(mtcars)
#'sum<-summary.numeric.mc(mtcars,latex=T)
#'@seealso \code{\link[base]{summary}}
#'@aliases summary.numeric.mc
#'@rdname summary.numeric.mc
#'@export summary.numeric.mc
summary.numeric.mc<-function(object,latex=FALSE){
  classes<-sapply(1:ncol(object), function(x) class(object[,x]))
  summary.numeric<-sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x]))) 
  if (length(which(is.na(object)))>0){
    summary.numeric<-sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x]))) 
    tmp <- object.frame()
    for(i in seq(along=summary.numeric)) for(j in 1:length(summary.numeric[[i]]))
      tmp[i,j] <- ifelse(is.na(summary.numeric[[i]][j]),"0",summary.numeric[[i]][j])
    summary.numeric<-tmp
    colnames(summary.numeric)<-c(names(summary(1)),"NA")
    summary.numeric[which(is.na(summary.numeric[,7])),7]<-"0"    
  } else {
    summary.numeric<-t(data.frame(sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x])))))
    colnames(summary.numeric)<-names(summary(1))
  }
  rownames(summary.numeric)<-colnames(object)[which(classes=="numeric")]
  if(latex){
    print(xtable(summary.numeric))
  }
  return(summary.numeric)
}



#'@docType methods
#'@title Normality Test
#'@description Perform a normality test for variables 
#'@param m variable (list, matrix, data frame...)
#'@param alpha p-value threshold
#'@return logical value indicating whether variable is normally distributed
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' normality.mc(mtcars)
#'@aliases normality.mc
#'@rdname normality.mc
#'@export normality.mc
normality.mc<-function(m,alpha=0.05){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m), function(x) shapiro.test(m[,x])$p.value<=alpha))
  else return(shapiro.test(m)$p.value<=alpha)
}


#'@docType methods
#'@title Object Classes
#'@description Get class of variable
#'@param variable (list, matrix, data frame...)
#'@return class of variables or of columns of matrix/data frame
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' class.mc(mtcars)
#'@aliases class.mc
#'@rdname class.mc
#'@export class.mc
class.mc<-function(m){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m),function(x) class(m[,x])))
  else return(class(m))
}