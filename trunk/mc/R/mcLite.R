#Package: mc
#Title: Commonly used functions for Nutriomics Team (INSERM U872)
#Version: 1.0
#Date: 2013-11-26
#Author: Aurelie Cotillard, Edi Prifti, Hoai Tuong Nguyen (A-Z order)
#Maintainer: Hoai Tuong Nguyen <hoai-tuong.nguyen@inserm.fr>
#Description: Statistical and datamining tools for omics data analysis.
#License: PPL

#Examples
if(FALSE){
  
  #load "mc" package
  library(mc)
  
  #load "xtable" package, automatically install the package if it does not exist, then load it
  library.mc("xtable")
  
  #read a large file
  #data<-read.table.mc("http://statistics.vn/data/doesgenes.txt",header=T,sep=";",nrow=1000)
  
  #get statistics on the columns of a matrix/data.frame and export the results as table to Latex codes
  attach(mtcars)
  sum<-summary.numeric.mc(mtcars,latex=T)
  
  #test the normality of a (list of) numeric variable(s)
  attach(mtcars)
  normality.mc(mtcars)
  
  #get class type for a (list of) variable(s)
  attach(mtcars)
  class.mc(mtcars)
  
  #plot the correlation, add lowess line to plot and write the output
  output.dir="../results"
  attach(swiss)
  reg.plot.mc(Fertility,Agriculture,
              type="lowess",
              pch=ifelse(swiss$Examination>10, 0, 1),         
              subjects=as.vector(rownames(swiss)),
              title="CORRELATION - LOWESS",
              xlab="Fertility",ylab="Agriculture",
              legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Examination>10","Examination<=10"),
                                   col=c("black","black")),
              imgfile=sprintf("%s/lw_swiss-Fertility-Agriculture.pdf",output.dir),
              pointsfile=sprintf("%s/lw_swiss-Fertility-Agriculture.csv",output.dir))  
  
  #draw a boxplot for two classes, add results of t-test to plot, write the output
  output.dir="../results"
  attach(lung)
  pdf(sprintf("%s/lung_factors_by_sex_boxplot.pdf",output.dir))
  outfile<-sprintf("%s/lung_factors_by_sex_t-test.csv",output.dir)
  par(mfrow = c(4, 4))
  lapply(c(1:4,6:10),function(x) boxplot.class.mc(data=lung,x,
                                                  class=lung$sex,
                                                  xlab="Sex (0=Female, 1=Male)",
                                                  outfile=outfile))
  dev.off()
  
  
}



#'@name check.installed.mc
#'@aliases check.installed.mc
#'@export check.installed.mc
#'@docType methods
#'@title Checking package installation
#'@description Check whether a packages is installed
#'@param pkg name of package
#'@return A logical value indicating whether the package is installed
#'@author Hoai Tuong Nguyen
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


#loading dependencies
library.mc("corrplot")
library.mc("bnlearn")
library.mc("FunNet")
library.mc("zoo")
library.mc("Hmisc")
library.mc("quantreg")
library.mc("gdata","cran")
library.mc("WriteXLS","cran")
library.mc("xtable")
library.mc("samr","bioc")
library.mc("ggplot2","bioc")
library.mc("reshape")
library.mc("scales")
library.mc("hgu95av2","bioc")
library.mc("hgu95av2.db","bioc")
library.mc("illuminaHumanv3.db","bioc")
library.mc("GO.db","bioc")
library.mc("sqldf")

#'@name read.table.mc
#'@aliases read.table.mc
#'@export read.table.mc
#'@docType methods
#'@title Data Input
#'@description Read a very large data file
#'@param file the name of the file which the data are to be read from.
#'@param header a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.
#'@param sep  the field separator character. Values on each line of the file are separated by this character. If sep = "" (the default for read.table) the separator is 'white space', that is one or more spaces, tabs, newlines or carriage returns. A (character) name of the column in a \code{\link[base]{data.frame}} which contains the addresses
#'@param nrow number of rows
#'@return A data frame (\code{\link[base]{data.frame}}) containing a representation of the data in the file
#'@author Hoai Tuong Nguyen
#'@seealso \code{\link[utils]{read.table}}
read.table.mc<-function(file,header=FALSE,sep="",nrow=-1){
  #read 5 first rows to get class names of column
  tab5rows <- read.table(file, nrows = 5,sep=sep)
  #get class names
  classes <- sapply(tab5rows, class)
  #get data frame with specific parameter
  tabAll <- read.table(file,  header=header, colClasses=classes,sep=sep,nrows=nrow,comment.char = "")
  return(tabAll)
}



#'@name summary.numeric.mc
#'@aliases summary.numeric.mc
#'@export summary.numeric.mc
#'@docType methods
#'@title Object Summaries
#'@description Summarize an numeric table, save the output to a table, export the output to Latex code
#'@param object an object for which a summary is desired.
#'@param latex a logical value indicating whether output to latex is called
#'@return Table of of the value returned by summary (and output Latex code if parameter 'latex' is TRUE)
#'@author Hoai Tuong Nguyen
#'@seealso \code{\link[base]{summary}}
summary.numeric.mc<-function(object,latex=FALSE){
  #get classes of columns
  classes<-sapply(1:ncol(object), function(x) class(object[,x]))
  #get summaries for numeric columns
  summary.numeric<-sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x]))) 
  #transform summaries into table format
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
  #export table into Latex code
  if(latex){
    print(xtable(summary.numeric))
  }
  return(summary.numeric)
}



#'@name normality.mc
#'@aliases normality.mc
#'@export normality.mc
#'@docType methods
#'@title Normality Test
#'@description Perform a normality test for variables 
#'@param m variable (list, matrix, data frame...)
#'@param alpha p-value threshold
#'@return logical value indicating whether variable is normally distributed
#'@author Hoai Tuong Nguyen
normality.mc<-function(m,alpha=0.05){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m), function(x) shapiro.test(m[,x])$p.value<=alpha))
  else return(shapiro.test(m)$p.value<=alpha)
}



#'@name class.mc
#'@aliases class.mc
#'@export class.mc
#'@docType methods
#'@title Object Classes
#'@description Get class of variable
#'@param variable (list, matrix, data frame...)
#'@return class of variables or of columns of matrix/data frame
#'@author Hoai Tuong Nguyen
class.mc<-function(m){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m),function(x) class(m[,x])))
  else return(class(m))
}


#'@name reg.plot.mc
#'@aliases reg.plot.mc
#'@export reg.plot.mc
#'@docType methods
#'@title X-Y Plotting
#'@description Plot a pair of variables and add regression line (linear or lowess) to plot
#'@param x a numeric vector 
#'@param y a numeric vector 
#'@param type type of regression line 
#'@param pch type points
#'@param subjects list of labels for points
#'@param title main title of plot
#'@param xlab a title for the x axis
#'@param ylab a title for the y axis
#'@param col list of colors for points
#'@param legend.topleft legend at the top-left of plot
#'@param legend.topright legend at the top-right of plot
#'@param imgfile image output filename
#'@param pointsfile points output filename
#'@author Hoai Tuong Nguyen
reg.plot.mc<-function(x,y,separator=NULL,type="none",quantile="outter",arrows=0,pch,subjects=NULL,title="CORRELATION",xlab="X",ylab="Y",col,legend.topleft,legend.topright,legend.bottomleft,legend.bottomright,imgfile=NULL,pointsfile=NULL){
  
  source("http://www.r-statistics.com/wp-content/uploads/2010/04/Quantile.loess_.r.txt")
  
  #correlation
  r=rcorr(x,y,type="pearson")[[1]][1,2]
  rho=rcorr(x,y,type="spearman")[[1]][1,2]
  
  #output plot
  if (!missing(imgfile))
    pdf(file=imgfile)
  
  #colors of points
  if (missing(col))
    col=rep("black",length(x),)
  
  inter=NULL
  lm=NULL
  lw=NULL
  lmlow=NULL
  lmup=NULL
  
  #Regression line
  if (type=="lowess"){
    if (quantile!="inter"){
      lw<-lowess(x,y)
      
      up<-which(x>lw$x & y>lw$y)
    } else {
      lmlow<-Quantile.loess(x, y,the.quant = 0.25)
      lw<-lowess(x,y)
      lmup<-Quantile.loess(x, y,the.quant = 0.75)
      
      
      up<-which(x>lw$x & y>lw$y)
      
      inter<-which(y>=fitted(lmlow) & y<=fitted(lmup))
      
    }
  } else if (type=="lm") {
    if (quantile!="inter"){
      if (is.null(separator)){
        lm<-lm(y~x)      
        up<-which(y>fitted(lm))        
      }else{
        fitted<-separator[2]*x+separator[1]
        up<-which(y>fitted)
      }

    } else {
      
      lmlow<-rq(y ~ x, tau = 0.25)
      lm<-lm(y~x)
      lmup<-rq(y ~ x, tau = 0.75)
      
      
      
      up<-which(y>fitted(lm))
      
      
      
      
      inter<-which(y>=fitted(lmlow) & y<=fitted(lmup))
      #inter<-which(residuals(lm(y ~x))>quantile(residuals(lm(y ~x)),c(0.25,0.5,0.75))[1]&residuals(lm(y ~x))<quantile(residuals(lm(y ~x)),c(0.25,0.5,0.75))[3])
    }
    
  }
  
  #output points
  if (type=="none")
    pointsfile="none"
  if (pointsfile!="none"){
    p<-rep("0",length(x))
    p[up]<-"1"
    p[inter]<-"-"
    
    p.class<-as.vector(p)
    p.class[which(p.class=="0")]<-0
    p.class[which(p.class=="1")]<-1
    p.class[which(p.class=="-")]<-NA
    p.class<-as.numeric(p.class)
    
    col.class<-as.vector(col)
    col.class[which(col.class=="blue")]<-0
    col.class[which(col.class=="red")]<-1
    col.class[which(col.class=="yellow")]<-NA
    col.class<-as.numeric(col.class)
    
    
  rgc<- cor(p.class,col.class,use="complete.obs")

  
  #main plot
  plot(x,y,pch=pch,main=sprintf("%s\n%s vs %s \nr.idx=%0.4f; rho.idx=%0.4f; r.gc=%0.4f",title,ylab,xlab,r,rho,rgc),xlab=xlab,ylab=ylab,col=col,cex=0.4)
  
  if (type=="lm" & !is.null(lmlow))
    abline(lmlow,col="blue",lty = 2)
  
  if (type=="lm" & !is.null(lmup))
    abline(lmup,col="green",lty = 2)
  
  if (type=="lowess" & !is.null(lmlow))
    points(lmlow$y.loess ~ lmlow$x, type = "l", col = "blue")
  
  if (type=="lowess" & !is.null(lmup))
    points(lmup$y.loess ~ lmup$x, type = "l", col = "green")
  
  if (type=="lm"){
    if (is.null(separator))
      abline(lm)
    else abline(separator)
  }
    
  if (type=="lowess")
    lines(lw,col=3)
  
  
  coltxt=rep("black",length(x),)
  coltxt[inter]<-"red"
  #labels
  text(x, y, subjects, cex=0.25,pos=1,offset=0.2,col=coltxt)
  
  
  if (arrows>0){
    nrow<-length(x)/(arrows+1)
    s<-1:nrow
    iter<-arrows
    while (iter>0){
      #arrows(x[s], y[s], x[s+nrow], y[s+nrow], col= ifelse(y[s]>y[s+nrow],"blue","red"),length=0.03)
      arrows(x[s], y[s], x[s+nrow], y[s+nrow], col= iter+20,length=0.03)
      s<-s+nrow
      iter<-iter-1
    }
      
  }
  
  
  #Legend  
  if (!missing(legend.topleft))
    legend("topleft", title=legend.topleft$title,pch=legend.topleft$pch, legend = legend.topleft$label, col = legend.topleft$col, cex=0.4)
  if (!missing(legend.topright))
    legend("topright", title=legend.topright$title,pch=legend.topright$pch, legend = legend.topright$label, col = legend.topright$col, cex=0.4)
  if (!missing(legend.bottomleft))
    legend("bottomleft", title=legend.topleft$title,pch=legend.bottomleft$pch, legend = legend.bottomleft$label, col = legend.bottomleft$col, cex=0.4)
  if (!missing(legend.bottomright))
    legend("bottomright", title=legend.bottomright$title,pch=legend.bottomright$pch, legend = legend.bottomright$label, col = legend.bottomright$col, cex=0.4)
  

    
    if (!missing(subjects))
      write.table(rbind(c("Name",xlab,ylab,"Levels"),cbind(subjects,x,y,p)),file=pointsfile,col.names=F,row.names=F,sep=";",quote=F)
    else 
      write.table(rbind(c(xlab,ylab,"Levels"),cbind(x,y,p)),file=pointsfile,col.names=F,row.names=T,sep=";",quote=F)
  }
  
  #end of output plot
  if (!missing(imgfile))
    dev.off()
  
}



#'@name risk.level.mc
#'@aliases risk.level.mc
#'@export risk.level.mc
#'@docType methods
#'@title X-Y Plotting
#'@description Plot a pair of variables and add regression line (linear or lowess) to plot
#'@param x a numeric vector 
#'@param y a numeric vector 
#'@param type type of regression line 
#'@param pch type points
#'@param subjects list of labels for points
#'@param title main title of plot
#'@param xlab a title for the x axis
#'@param ylab a title for the y axis
#'@param col list of colors for points
#'@param legend.topleft legend at the top-left of plot
#'@param legend.topright legend at the top-right of plot
#'@param imgfile image output filename
#'@param pointsfile points output filename
#'@author Hoai Tuong Nguyen
risk.level.mc<-function(x,y,type="lm",pch,subjects=NULL,title="CORRELATION",xlab="X",ylab="Y",col,legend.topleft,legend.topright,imgfile=NULL,pointsfile=NULL){
  
  #correlation
  r=rcorr(x,y,type="pearson")[[1]][1,2]
  rho=rcorr(x,y,type="spearman")[[1]][1,2]
  
  #output plot
  if (!missing(imgfile))
    pdf(file=imgfile)
  
  #colors of points
  if (missing(col))
    col=rep("black",length(x),)
  
  #main plot
  plot(x,y,pch=pch,main=sprintf("%s\n%s vs %s \nr=%0.4f; rho=%0.4f",title,ylab,xlab,r,rho),xlab=xlab,ylab=ylab,col=col)
  
  #labels
  text(x, y, subjects, cex=0.3,pos=1,offset=0.2)
  
  inter=NULL
  
  #Regression line
  if (type=="lowess"){
    lw<-lowess(x,y)
    lines(lw,col=3)
    up<-which(x>lw$x & y>lw$y)
  } else if (type=="lm") {
    if (quantile!="inter"){
      lm<-lm(y~x)
      abline(lm)
      up<-which(y>fitted(lm))
    } else {
      
      lmlow<-rq(y ~ x, tau = 0.25)
      lm<-lm(y~x)
      lmup<-rq(y ~ x, tau = 0.75)
      
      abline(lmlow,col="red",lty = 2)
      abline(lm)
      abline(lmup,col="green",lty = 2)
      up<-which(y>fitted(lm))
      inter<-which(y>fitted(lmlow) & y<fitted(lmup))
    }
    
  }
  
  #Legend  
  if (!missing(legend.topleft))
    legend("topleft", title=legend.topleft$title,pch=legend.topleft$pch, legend = legend.topleft$label, col = legend.topleft$col, cex=0.4)
  if (!missing(legend.topright))
    legend("topright", title=legend.topright$title,pch=legend.topright$pch, legend = legend.topright$label, col = legend.topright$col, cex=0.4)
  
  
  #output points
  if (!missing(pointsfile)){
    p<-rep("0",length(x))
    p[up]<-"1"
    p[inter]<-"-"
    if (!missing(subjects))
      write.table(rbind(c("Name",xlab,ylab,"Levels"),cbind(subjects,x,y,p)),file=pointsfile,col.names=F,row.names=F,sep=";",quote=F)
    else 
      write.table(rbind(c(xlab,ylab,"Levels"),cbind(x,y,p)),file=pointsfile,col.names=F,row.names=T,sep=";",quote=F)
  }
  
  #end of output plot
  if (!missing(imgfile))
    dev.off()
  
}


#'@name boxplot.class.mc
#'@aliases boxplot.class.mc
#'@export boxplot.class.mc
#'@docType methods
#'@title Box Plots 
#'@description Draw a boxplot for a column of data frame
#'@param data a data frame
#'@param x index of column
#'@param class vector of classes
#'@param outfile output filename
#'@author Hoai Tuong Nguyen
boxplot.class.mc<-function(data,x,type="auto",class,xlab,ylab,outfile=NULL){
  
  
  if(missing(ylab))
    ylab=names(data)[x]
  if(missing(xlab))
    xlab=names(class)
  nskip1 <-!sum(!is.na(data[,x]))<=3
  #nskip2 <-!sum(data[1,x]==data[,x]) == nrow(data)
  print(data[,x])
  print(nskip1)
  print(nlevels(as.factor(data[,x])))

  if (!is.na(data[,x]) & nskip1 & (nlevels(as.factor(data[,x]))>1)){
    
    boxplot(data[,x]~class,ylab=ylab,xlab=xlab) 
    if(type=="t")
      t.res<-t.test(data[,x]~class,ylab=names(data)[x])
    if(type=="mwu")
      t.res<-wilcox.test(data[,x]~class,ylab=names(data)[x])
    if(type=="auto")
      if (shapiro.test(data[,x])$p.value<=0.05)
        t.res<-t.test(data[,x]~class,ylab=names(data)[x])
    else t.res<-wilcox.test(data[,x]~class,ylab=names(data)[x])
    
    cor.res<-cor(data[,x],class,use="complete.obs")
    
    title(main=sprintf("t=%0.2f; r=%.3f; p=%0.2e\n%s",t.res$statistic,cor.res,t.res$p.value,ifelse(t.res$p.value<=0.001,"***",ifelse(t.res$p.value<=0.01 & t.res$p.value>0.001,"**",ifelse(t.res$p.value<=0.05 & t.res$p.value>0.01,"*","")))),cex=0.5)  
    
    out<-cbind(colnames(data)[x],t.res$statistic,cor.res,t.res$p.value,
               ifelse(t.res$p.value<=0.001,"***",ifelse(t.res$p.value<=0.01 & t.res$p.value>0.001,"**",ifelse(t.res$p.value<=0.05 & t.res$p.value>0.01,"*",""))),
               shapiro.test(data[,x])$p.value<=0.05)
    
  } else 
    out<-cbind(colnames(data)[x],"NA","NA","NA","NA")
  
  if(!missing(outfile)){
    write.table(out,outfile,col.names=F,row.names=F,append=T,quote=F,sep=";")
  }   
}


#'@name yule.Q.mc
#'@aliases yule.Q.mc
#'@export yule.Q.mc
#'@docType methods
#'@title Correlation between nominal variables 
#'@description A measure of correlation between nominal variables
#'@param x a nominal vector
#'@param y a nominal vector
#'@author Hoai Tuong Nguyen
yule.Q.mc=function(x,y){(table(x,y)[1,1]*table(x,y)[2,2]-table(x,y)[1,2]*table(
  x,y)[2,1])/(table(x,y)[1,1]*table(x,y)[2,2]+table(x,y)[1,2]*table(x,y)[2,1])
}

library.mc("samr","bioc")

#'@name samt.mc
#'@aliases samt.mc
#'@export samt.mc
#'@docType methods
#'@title Adapted version of samr 
#'@description An adapted version of samr for gene differtiation analysis
#'@param df a microarray dataframe 
#'@param label vector of label
#'@param nperms number of permutations
#'@param logged2 logical. Expression level has been transformed by logorith base 2
#'@author Hoai Tuong Nguyen
samt.mc<-function(df,label,nperms=100,logged2=T,type="Two class unpaired",seq=F){
  if(logged2)
    df<-log2(df)
  data<-list(x=df,y=label, geneid=as.character(1:nrow(df)),genenames=rownames(df), logged2=TRUE)
  if (!seq){
    if (type=="Two class unpaired")
      samr.obj<-samr(data,  resp.type=type, nperms=nperms,testStatistic="wilcoxon")
    else samr.obj<-samr(data,  resp.type=type, nperms=nperms)
  }else {
    if (type=="Two class unpaired")
      samr.obj<-SAMseq(x=df,y=label,  resp.type=type, nperms=nperms,testStatistic="wilcoxon",fdr.output = 0.20)
    else samr.obj<-SAMseq(x=df,y=label,  resp.type=type, nperms=nperms,fdr.output = 0.20)
  }
  
  delta.table <- samr.compute.delta.table(samr.obj)
  delta<-0.45
  siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
  
  pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
  return(list(probeID=rbind(siggenes.table$genes.up,siggenes.table$genes.lo)[,2],
              UpLow=c(rep("UP",siggenes.table$ngenes.up),rep("LOW",siggenes.table$ngenes.lo)),
              Stat=rbind(siggenes.table$genes.up,siggenes.table$genes.lo)[,4],
              RawpValue=pv[rbind(siggenes.table$genes.up,siggenes.table$genes.lo)[,2]],
              FoldChange=rbind(siggenes.table$genes.up,siggenes.table$genes.lo)[,7],
              FDR=rbind(siggenes.table$genes.up,siggenes.table$genes.lo)[,8],
              nsiggenes=siggenes.table$ngenes.up+siggenes.table$ngenes.lo))
}

#'@name samt.opt.mc
#'@aliases samt.opt.mc
#'@export samt.opt.mc
#'@docType methods
#'@title Optimized version of samr 
#'@description An optimized version of samr for gene differtiation analysis
#'@param df a microarray dataframe 
#'@param label vector of label
#'@param nmax.perms number of maximum permutations
#'@param logged2 logical. Expression level has been transformed by logorith base 2
#'@author Hoai Tuong Nguyen
samt.opt.mc<-function(df,label,nmax.perms=100, logged2=T,type="Two class unpaired",seq=F){
  if(logged2)
    df<-log2(df)
  
  bins<-seq(10,nmax.perms,10)
  
  out.log<-do.call(rbind,sapply(bins,function(t){
    samres<-samt.mc(df, label,t,logged2,type,seq=seq)
    list(c(nperms=t,nsiggenes=samres$nsiggenes,samres=samres))
  }
  ))
  #nperms.opt<-as.numeric(out.log[which(unlist(out.log[,2])==max(unlist(out.log[,2]))),1])
  
  #samres.opt<-out.log[which(unlist(out.log[,2])==max(unlist(out.log[,2]))),3:ncol(out.log)]
  #print(paste("Optimal number of permuations:",nperms.opt))
  #print(out.log)
  plot(out.log[,1:2])
  
  return(out.log)
}




#'@name samr.opt.cluster.mc
#'@aliases samr.opt.cluster.mc
#'@export samr.opt.cluster.mc
#'@docType methods
#'@title Optimized version of samr with clusters
#'@description An optimized version of samr with cluster for gene differtiation analysis
#'@param df a microarray dataframe 
#'@param label vector of label
#'@param nmax.perms number of maximum permutations
#'@param logged2 logical. Expression level has been transformed by logorith base 2
#'@author Hoai Tuong Nguyen
samr.opt.cluster.mc<-function(df,class,label=NULL,file,index,cluster,type="Two class unpaired",logged2=T,nmax.perms=200,seq=F){
  print(cluster)
  print(is.null(nrow(df)))
  if (!is.null(nrow(df)))
    if ( nrow(df)>=10){
      
      if (exists("class")&is.null(label)){
        print("OK")
        label<-as.vector(class[colnames(df),4])
        label[label==0]<-2
      }
      
      sig<-samt.opt.mc(df, 
                       label,
                       nmax.perms=nmax.perms,
                       logged2=logged2,
                       type=type,
                       seq=seq)
      
      i.opt<-which.max(unlist(sig[,2]))
      
      if (max(unlist(sig[,2]))>0){
        sig.opt<-sig[i.opt,]
        res.sig<-data.frame(sig.opt$samres.probeID,
                            rep(index,length(sig.opt$samres.probeID)),
                            sig.opt$samres.UpLow,
                            sig.opt$samres.Stat,
                            sig.opt$samres.RawpValue,
                            sig.opt$samres.FoldChange,
                            sig.opt$samres.FDR)
        colnames(res.sig)<-c("probeID", "Index", "UpLow", "Stat",    "RawpValue",   "FoldChange", "FDR(%)")
        write.table(res.sig,file=sprintf("%s_%s.csv",file,cluster),col.names=T,row.names=F,quote=F,sep="\t")    
      }
      return(sig)
    }
  
}







#'@name get.nona.mc
#'@aliases get.nona.mc
#'@export get.nona.mc
#'@docType methods
#'@title Filtering missing value genes 
#'@description Filtering missing value genes
#'@param df a microarray dataframe 
#'@param rate percentage of acceptable non-missing values per gene
#'@author Hoai Tuong Nguyen
get.nona.mc<-function(df,rate){
  index.nona<-apply(df, 1, function(x) sum(is.na(x)))/ncol(df)<=(1 - rate/100)
  return(df[index.nona,])
}


#'@name addsigText.mc
#'@aliases addsigText.mc
#'@export addsigText.mc
#'@docType methods
#'@title Adding significant bar
#'@description Adding a significant bar to the barplot
#'@param x0,y0   coordinates of points from which to draw
#'@param x1,y1   coordinates of points to which to draw
#'@param lab label to be added
#'@author Hoai Tuong Nguyen
addsigText.mc<-function(x0,y0,x1,y1,lab="",offset){
  arrows(x0=x0, y0=max(y0,y1)+offset, x1 = x1, y1 = max(y0,y1)+offset,code=0)
  text(x=mean(c(x0,x1)),y=max(y0,y1)+offset,label=lab,po=3,xpd=TRUE) 
  text(x=x0, y=max(y0,y1)+offset,label="|",xpd=TRUE) 
  text(x=x1, y=max(y0,y1)+offset,label="|",xpd=TRUE)
}

#'@name barplotStar.mc
#'@aliases barplotStar.mc
#'@export barplotStar.mc
#'@docType methods
#'@title Bar Plots with significant bar 
#'@description Creates a bar plot with vertical or horizontal bars and significant stars.
#'@param df  matrix of values of columns
#'@param title title of bar plot
#'@param col color of two columns
#'@param xlab  a label for the x axis
#'@param ylab  a label for the y axis
#'@param labels series of columns' names
#'@author Hoai Tuong Nguyen
barplotStar.mc<-function(df,title=" ",col=c("darkblue","red"),xlab="",ylab="",labels){
  barX <-barplot(df, 
                 ylim=c(0,max(df)+max(df)/5),
                 main=title,
                 xlab=xlab, 
                 ylab=ylab,
                 col=col, 
                 names.arg=labels,
                 beside=TRUE)
  #text(x=barX,y=df,label=df,po=3,xpd=TRUE) 
  sapply(1:ncol(barX), function(x) addsigText.mc(x0=barX[1,x],y0=df[1,x],x1=barX[2,x],y1=df[2,x],lab="***",offset=max(df)/10))
  
}


#'@name heatmap.mc
#'@aliases heatmap.mc
#'@export heatmap.mc
#'@docType methods
#'@title Draw a Heat Map 
#'@description Creates a Heat Map from a matrix of continous values
#'@param df  matrix of values of columns
#'@param col color of highest value
#'@author Hoai Tuong Nguyen
heatmap.mc<-function(df,col){
  df.m <- melt(df)
  colnames(df.m)<-c("Name","variable","value")
  df.m <- ddply(df.m, .(variable), transform,   rescale = rescale(value))
  (p <- ggplot(df.m, aes(variable, Name))
   + geom_tile(aes(fill = rescale), colour = "white")
   + scale_fill_gradient(low = "white",  high = col))
  return(p)
}


#'@name FunNet.mc
#'@aliases FunNet.mc
#'@export FunNet.mc
#'@docType methods
#'@title FunNet adapted version 
#'@description FunNet adapted version with customized output filenames
#'@author Hoai Tuong Nguyen
FunNet.mc<-function (wd = "", org = "hsa", two.lists = TRUE, up.frame = NULL, 
                     down.frame = NULL, genes.frame = NULL, restrict = FALSE, 
                     ref.list = NULL, logged = FALSE, discriminant = FALSE, go.bp = TRUE, 
                     go.cc = TRUE, go.mf = TRUE, kegg = TRUE, annot.method = "specificity", 
                     annot.details = TRUE, direct = FALSE, enriched = TRUE, fdr = NA, 
                     build.annot.net = TRUE, coexp.matrix = NULL, coexp.method = "spearman", 
                     estimate.th = FALSE, hard.th = NA, soft.th = NA, topological = FALSE, 
                     keep.sign = FALSE, level = NA, annot.clust.method = "umilds", 
                     annot.prox.measure = "unilat.pond.norm.mean", test.recovery = FALSE, 
                     test.robust = FALSE, replace.annot = NA, random.annot = FALSE, 
                     build.gene.net = FALSE, gene.clust.method = "hclust", gene.net.details = FALSE, 
                     gene.clusters = NA, alpha = 0.05, RV = 0.9, sigma = NA, keep.rdata = FALSE, 
                     zip = TRUE,fileprefix) 
{
  if (org == "HS") {
    org <- "hsa"
  }
  if (org == "MM") {
    org <- "mmu"
  }
  if (org == "RN") {
    org <- "rno"
  }
  if (org == "SC") {
    org <- "sce"
  }
  parameter.list <- list(analysis.date = date(), package.version = .funnet.version, 
                         org = org, annot.date = annot.date, annot.method = annot.method, 
                         annot.clust.method = annot.clust.method, annot.prox.measure = annot.prox.measure, 
                         direct = direct, enriched = enriched, fdr = fdr, two.lists = two.lists, 
                         restrict = restrict, go.bp = go.bp, go.cc = go.cc, go.mf = go.mf, 
                         kegg = kegg, discriminant = discriminant, logged = logged, 
                         annot.details = annot.details, estimate.th = estimate.th, 
                         hard.th = hard.th, soft.th = soft.th, coexp.method = coexp.method, 
                         topological = topological, keep.sign = keep.sign, build.annot.net = build.annot.net, 
                         level = level, test.recovery = test.recovery, test.robust = test.robust, 
                         replace.annot = replace.annot, random.annot = random.annot, 
                         build.gene.net = build.gene.net, gene.clust.method = gene.clust.method, 
                         gene.clusters = gene.clusters, gene.net.details = gene.net.details, 
                         alpha = alpha, RV = RV, sigma = sigma, keep.rdata = keep.rdata, 
                         zip = zip)
  .check.parameters(parameter.list, coexp.matrix, up.frame, 
                    down.frame, ref.list, genes.frame)
  if (discriminant) {
    two.lists <- TRUE
    restrict <- TRUE
  }
  if (!is.null(up.frame) & !is.null(down.frame) & is.null(genes.frame)) {
    two.lists <- TRUE
  }
  if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
    two.lists <- FALSE
    discriminant <- FALSE
  }
  cat(paste("\n\tFunNet started at: ", date(), sep = ""))
  cat(paste("\n\t\tUsing annotations updated on: ", annot.date, 
            sep = ""))
  if (wd != "") {
    setwd(wd)
  }
  results.dir <- paste(fileprefix,"_", format(Sys.time(), "%Y_%b_%d_%H-%M-%S"), 
                       sep = "")
  dir.create(paste(getwd(), "/", results.dir, sep = ""))
  dir.create(paste(getwd(), "/", results.dir, "/html", sep = ""))
  dir.create(paste(getwd(), "/", results.dir, "/images", sep = ""))
  try(write.table(as.matrix(print(parameter.list)), col.names = F, 
                  file = paste(getwd(), "/", results.dir, "/", "parameters_list.txt", 
                               sep = ""), sep = "\t"))
  wd <- getwd()
  locus.name <- annot.base[[org]]$locus.name[, 1:2]
  locus.symbol <- annot.base[[org]]$locus.name[, c(1, 3)]
  rownames(locus.name) <- locus.name[, 1]
  rownames(locus.symbol) <- locus.symbol[, 1]
  if (two.lists) {
    up.down <- .filter.genes(up.frame = up.frame, down.frame = down.frame, 
                             two.lists = TRUE, locus.name = locus.name, logged = logged)
    up.frame <- up.down$up.frame
    down.frame <- up.down$down.frame
    rm(up.down)
    if (discriminant) {
      ref.list <- c(as.character(up.frame[, 1]), as.character(down.frame[,1]))
    }
    else if (restrict & !discriminant) {
      ref.list <- .filter.genes(restrict = TRUE, ref.list = ref.list, 
                                locus.name = locus.name)
    }
    else if (!restrict & !discriminant) {
      ref.list <- NULL
    }
  }
  else {
    genes.frame <- .filter.genes(genes.frame = genes.frame, 
                                 two.lists = FALSE, locus.name = locus.name)
    if (restrict) {
      ref.list <- .filter.genes(restrict = TRUE, ref.list = ref.list, 
                                locus.name = locus.name)
    }
    else {
      ref.list <- NULL
    }
  }
  cat(paste("\n\tSaving start-up environment... ", format(Sys.time(), 
                                                          "%X"), sep = ""))
  save(up.frame, down.frame, ref.list, genes.frame, parameter.list, 
       coexp.matrix, locus.name, file = paste(getwd(), "/", 
                                              results.dir, "/", "start-up_environment.RData", sep = ""), 
       compress = T)
  save(up.frame, down.frame, ref.list, genes.frame, parameter.list, 
       coexp.matrix, locus.name, file = paste(getwd(), "/", 
                                              "start-up_environment.RData", sep = ""), compress = T)
  if (estimate.th) {
    datas <- NULL
    if (!is.null(up.frame) & !is.null(down.frame)) {
      datas <- rbind(up.frame, down.frame)
    }
    if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
      datas <- genes.frame
    }
    rownames(datas) <- datas[, 1]
    datas <- datas[, 2:ncol(datas)]
    cat("\n\tHard thresholding...\n")
    try(hard.th <- .PickHardThreshold(datExpr1 = t(datas), 
                                      coexp.method = coexp.method))
    try(write.table(hard.th$tablou, file = paste(getwd(), 
                                                 "/", results.dir, "/", coexp.method, "_hard_threshold.txt", 
                                                 sep = ""), append = FALSE, col.names = TRUE, , row.names = F, 
                    sep = "\t"))
    try(write(paste("\n\nHard threshold estimate: ", hard.th$estimate, 
                    "\n", sep = ""), file = paste(getwd(), "/", results.dir, 
                                                  "/", coexp.method, "_hard_threshold.txt", sep = ""), 
              append = TRUE))
    try(write(paste("Do not trust this automated estimation without checking it!\n", 
                    "Do not hesitate to select another threshold depending on the associated connectivity values.\n", 
                    "Then please restart FunNet interaction analysis with your selected threshold.\n", 
                    sep = ""), file = paste(getwd(), "/", results.dir, 
                                            "/", coexp.method, "_hard_threshold.txt", sep = ""), 
              append = TRUE))
    cat("\n\tSoft thresholding...\n")
    try(soft.th <- .PickSoftThreshold(datExpr1 = t(datas), 
                                      coexp.method = coexp.method))
    try(write.table(soft.th$tablou, file = paste(getwd(), 
                                                 "/", results.dir, "/", coexp.method, "_soft_threshold.txt", 
                                                 sep = ""), append = FALSE, col.names = TRUE, , row.names = F, 
                    sep = "\t"))
    try(write(paste("\n\nSoft threshold estimate: ", soft.th$estimate, 
                    "\n", sep = ""), file = paste(getwd(), "/", results.dir, 
                                                  "/", coexp.method, "_soft_threshold.txt", sep = ""), 
              append = TRUE))
    try(write(paste("Do not trust this automated estimation without checking it!\n", 
                    "Do not hesitate to select another threshold depending on the associated connectivity values.\n", 
                    "Then please restart FunNet interaction analysis with your selected threshold.\n", 
                    sep = ""), file = paste(getwd(), "/", results.dir, 
                                            "/", coexp.method, "_soft_threshold.txt", sep = ""), 
              append = TRUE))
    print("Estimation of the co-expression threshold finished!")
    print("Please restart FunNet with your chosen threshold.")
    if (!keep.rdata) {
      try(unlink(paste(getwd(), "/", results.dir, "/", 
                       list.files(path = paste(getwd(), "/", results.dir, 
                                               "/", sep = ""), pattern = "[:print:]*.RData"), 
                       sep = ""), recursive = TRUE))
    }
    if (zip) {
      try(unlink(paste(getwd(), "/", results.dir, "/html", 
                       sep = ""), recursive = TRUE))
      try(unlink(paste(getwd(), "/", results.dir, "/images", 
                       sep = ""), recursive = TRUE))
      try(system(command = paste("zip -r9q ", results.dir, 
                                 ".zip ", "./", results.dir, "/*", sep = "")))
      try(unlink(paste(getwd(), "/", results.dir, sep = ""), 
                 recursive = TRUE))
    }
    options(show.error.messages = FALSE)
    stop()
  }
  if (is.null(coexp.matrix) & (build.annot.net | build.gene.net)) {
    cat(paste("\n\tComputing co-expression matrix... ", format(Sys.time(), 
                                                               "%X"), sep = ""))
    datas <- NULL
    if (!is.null(up.frame) & !is.null(down.frame)) {
      datas <- rbind(up.frame, down.frame)
    }
    if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
      datas <- genes.frame
    }
    rownames(datas) <- datas[, 1]
    datas <- datas[, 2:ncol(datas)]
    
    if (coexp.method %in% c("spearman", "pearson", "kendall")) {
      coexp.matrix <- rcorr(t(datas), type = coexp.method)$r
      
    }
    else if (coexp.method == "euclid") {
      coexp.matrix <- 1 - (as.matrix(dist(datas))/max(dist(datas), 
                                                      na.rm = TRUE))
    }
    try(save(coexp.matrix, file = paste(getwd(), "/", results.dir, 
                                        "/", "brut_coexp_matrix.RData", sep = ""), compress = T))
  }
  if (build.annot.net | build.gene.net) {
    sign.matrix <- coexp.matrix/abs(coexp.matrix)
    coexp.matrix <- abs(coexp.matrix)
    if (annot.clust.method %in% c("umilds", "spectral")) {
      if (!is.na(hard.th)) {
        coexp.matrix[coexp.matrix >= hard.th] <- 1
        coexp.matrix[coexp.matrix < hard.th] <- 0
      }
      else if (!is.na(soft.th)) {
        coexp.matrix <- coexp.matrix^soft.th
      }
      if (topological) {
        coexp.matrix <- 1 - .TOMdist(adjmat1 = coexp.matrix)
      }
      else if (keep.sign) {
        coexp.matrix <- coexp.matrix * sign.matrix
      }
    }
    else if (annot.clust.method == "ucknn" & !is.na(hard.th)) {
      coexp.matrix[coexp.matrix < hard.th] <- 0
    }
    try(save(coexp.matrix, sign.matrix, parameter.list, file = paste(getwd(), 
                                                                     "/", results.dir, "/", "gene_adj_matrix.RData", sep = "")))
  }
  if (kegg) {
    terms.name <- KEGG.terms.name
    rownames(terms.name) <- terms.name[, 1]
    file.annot <- annot.base[[org]]$KEGG.file.annot
    taxoname <- "KEGG"
    if (two.lists) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = FALSE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = NA, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = FALSE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = NA, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  go.name <- c("GO Biological Process", "GO Cellular Component", 
               "GO Molecular Function")
  terms.name <- GO.terms.name
  rownames(terms.name) <- terms.name[, 1]
  if (go.bp) {
    file.annot <- annot.base[[org]]$GO.DIR.BP.file.annot
    taxoname <- go.name[1]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (go.cc) {
    file.annot <- annot.base[[org]]$GO.DIR.CC.file.annot
    taxoname <- go.name[2]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (go.mf) {
    file.annot <- annot.base[[org]]$GO.DIR.MF.file.annot
    taxoname <- go.name[3]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (build.gene.net & !is.null(coexp.matrix)) {
    try(clusters <- .build.coexp.net(coexp.matrix = coexp.matrix, 
                                     locus.name = locus.name, locus.symbol = locus.symbol, 
                                     gene.clust.method = gene.clust.method, gene.clusters = gene.clusters))
    try(save(clusters, parameter.list, file = paste(getwd(), 
                                                    "/", results.dir, "/", "co-expression_clusters.RData", 
                                                    sep = ""), compress = T))
    net.matrix <- coexp.matrix
    rownames(net.matrix) <- locus.symbol[rownames(net.matrix), 
                                         2]
    colnames(net.matrix) <- locus.symbol[colnames(net.matrix), 
                                         2]
    try(.cyto.sym(net.matrix = net.matrix, file.net = paste(getwd(), 
                                                            "/", results.dir, "/", "co-expression_net.txt", sep = ""), 
                  diagonal = FALSE, thresh = NULL))
    rm(net.matrix)
    try(centrality <- .genes.centrality(adj.matrix = coexp.matrix, 
                                        clusters = clusters, taxoname = taxoname, locus.symbol = locus.symbol, 
                                        results.dir = results.dir, coexp = TRUE))
    if (two.lists) {
      up.down <- rbind(matrix(1, nrow(up.frame), 1), matrix(0, 
                                                            nrow(down.frame), 1))
      rownames(up.down) <- c(as.character(up.frame[, 1]), 
                             as.character(down.frame[, 1]))
      try(write.table(cbind(rownames(clusters$gene.connect), 
                            as.vector(locus.symbol[rownames(clusters$gene.connect), 
                                                   2]), as.vector(locus.name[rownames(clusters$gene.connect), 
                                                                             2]), up.down, clusters$gene.connect, centrality[rownames(clusters$gene.connect), 
                                                                                                                             ]), file = paste(getwd(), "/", results.dir, 
                                                                                                                                              "/", "co-expression_net_info.txt", sep = ""), 
                      sep = "\t", col.names = c("geneid", "symbol", 
                                                "name", "up(1)_down(0)", colnames(clusters$gene.connect), 
                                                colnames(centrality)), row.names = F))
      try(rm(clusters, up.down, centrality))
    }
    else {
      try(write.table(cbind(rownames(clusters$gene.connect), 
                            as.vector(locus.symbol[rownames(clusters$gene.connect), 
                                                   2]), as.vector(locus.name[rownames(clusters$gene.connect), 
                                                                             2]), clusters$gene.connect, centrality[rownames(clusters$gene.connect), 
                                                                                                                    ]), file = paste(getwd(), "/", results.dir, 
                                                                                                                                     "/", "co-expression_net_info.txt", sep = ""), 
                      sep = "\t", col.names = c("geneid", "symbol", 
                                                "name", colnames(clusters$gene.connect), colnames(centrality)), 
                      row.names = F))
      try(rm(clusters, centrality))
    }
    cat(paste("\n\tCo-expression net building finished... ", 
              date(), sep = ""))
    rm()
  }
  if (!keep.rdata) {
    try(unlink(paste(getwd(), "/", results.dir, "/", list.files(path = paste(getwd(), 
                                                                             "/", results.dir, "/", sep = ""), pattern = "[:print:]*.RData"), 
                     sep = ""), recursive = TRUE))
  }
  if (zip) {
    try(system(command = paste("zip -r9q ", results.dir, 
                               ".zip ", "./", results.dir, "/*", sep = "")))
    try(unlink(paste(getwd(), "/", results.dir, sep = ""), 
               recursive = TRUE))
  }
  cat(paste("\n\tEnd  of treatment at: ", date(), "\n", sep = ""))
  rm()
}


#'@name signs.mc
#'@aliases signs.mc
#'@export signs.mc
#'@docType methods
#'@title Get signs of changes 
#'@description Get signs of changes for a variable at two timepoints
#'@param t0  values of variable at t0
#'@param tn  values of variable at tn
#'@author Hoai Tuong Nguyen
signs.mc<-function(t0,tn){
  delta<-tn-t0
  delta[delta<0]<-"-"
  delta[delta>0]<-"+"
  delta[delta==0]<-"="
  return(delta)  
}


#'@name cause.mc
#'@aliases cause.mc
#'@export cause.mc
#'@docType methods
#'@title Causality analysis
#'@description Causality analysis
#'@param x0  values of variable at x0
#'@param xn  values of variable at xn
#'@param y0  values of variable at y0
#'@param yn  values of variable at yn
#'@author Hoai Tuong Nguyen
cause.mc<-function(x0,xn,y0,yn){
  signs.x<-signs.mc(x0,xn)
  signs.y<-signs.mc(y0,yn)  
  signs.sum<-summary(as.factor(paste(signs.x,signs.y,sep="")))
  max.probs<-signs.sum/length(signs.x)
  max.probs.names<-names(max.probs)
  return(list(signs.sum=signs.sum,max.probs=max.probs,max.probs.names=max.probs.names))  
}




#'@name updateAnnot.mc
#'@aliases updateAnnot.mc
#'@export updateAnnot.mc
#'@docType methods
#'@title Update Annotation Database
#'@description Causality analysis
#'@param dbfile  actual DB
#'@param record  new record
#'@author Hoai Tuong Nguyen
updateAnnot.mc<-function(dbfile,record){
  Annot.DB<-get(load(file=dbfile))
  colnames(record)<-c("PROBEID_ILMN","ENTREZID","SYMBOL","GENENAME")
  new.record.idx<-which(!record[,1]%in%as.vector(Annot.DB[,1]))
  if (length(new.record.idx)>0)
    Annot.DB<-rbind(Annot.DB,record[new.record.idx,])
  return(Annot.DB)
}


#'@name bn.bst.mc
#'@aliases bn.bst.mc
#'@export bn.bst.mc
#'@docType methods
#'@title Learning Bayesian networks with bootstraping
#'@description Learning Bayesian networks with bootstraping
#'@param data dataframe
#'@param n  number of sample
#'@author Hoai Tuong Nguyen
bn.bst.mc<-function(data,n){
  
  net.bst<-do.call(rbind,sapply(1:n,function(x){
    print(paste("Learning iteration:",x))
    net<-iamb(data[sample(nrow(data),replace=T),]) 
    list(net$arcs)
  }))
  
  #nodes<-nodes(net.bst)
  nodes<-union(unique(net.bst[,1]),unique(net.bst[,2]))
  nodes<-nodes[order(nodes)]
  
  cooc.mat<-sapply(1:length(nodes), function(i)
    sapply(1:length(nodes), function(j){
      noc<-nrow(net.bst[net.bst[,1]==nodes[i]&net.bst[,2]==nodes[j],])
      list(ifelse(is.null(noc),0,noc))
    }
    )
  )
  
  cooc.mat<-t(cooc.mat)
  rownames(cooc.mat)<-colnames(cooc.mat)<-nodes
  
  cooc.mat<-as.numeric.mat.mc(cooc.mat)/n
  list<-mat2list.mc(cooc.mat)
  
  #toplist<-list[list[,3]>0.5,]
  toplist<-list
  
  from<-as.character(toplist[,1])
  to<-as.character(toplist[,2])
  weight<-as.numeric(toplist[,3])
  toplist.df<-data.frame(cbind(from,to,as.numeric(weight)),stringsAsFactors=FALSE)
  colnames(toplist.df)<-c("f","t","weight")
  
  return(toplist.df)
}


#'@name as.numeric.mat.mc
#'@aliases as.numeric.mat.mc
#'@export as.numeric.mat.mc
#'@docType methods
#'@title Convert matrix to numeric matrix
#'@description Convert matrix to numeric matrix
#'@param data dataframe
#'@param dec  decimal of number
#'@author Hoai Tuong Nguyen
as.numeric.mat.mc<-function(data,dec=","){
  rn<-rownames(data)
  cn<-colnames(data)
  data<-data.frame(sapply(1:ncol(data),function(x) type.convert(as.character(data[,x]),dec=dec)))
  rownames(data)<-rn
  colnames(data)<-cn
  return(data)
}


#'@name mat2list.mc
#'@aliases mat2list.mc
#'@export mat2list.mc
#'@docType methods
#'@title Convert matrix to list
#'@description Convert matrix to list
#'@param data dataframe
#'@return a list with weight
#'@author Hoai Tuong Nguyen
mat2list.mc<-function(data,alpha=0){
  
  index.mc<-function(k,m){
    i=(k-1)%%m+1
    j=(k-1)%/%m+1
    return(list(i,j))
  }
  
  getWeight.mc<-function(x,data){
    i=index.mc(x,nrow(data))[[1]]
    j=index.mc(x,nrow(data))[[2]]
    w=as.numeric(data[i,j])
    return(list(i=rownames(data)[i],j=colnames(data)[j],w=as.numeric(data[i,j])))
  }
  
  
  
  #return(t(sapply(1:(nrow(data)*ncol(data)),function(x) getWeight.mc(x,data))))
  res<-t(sapply(which(data>alpha),function(x) getWeight.mc(x,data)))
  
  from<-as.character(res[,1])
  to<-as.character(res[,2])
  weight<-as.numeric(res[,3])
  res.df<-cbind(from,to,weight)
  res.df<-data.frame(res.df,stringsAsFactors=FALSE)
  colnames(res.df)=c("f","t", "weight")
  
  
  
  return(res.df)
}

#'@name list2mat.mc
#'@aliases list2mat.mc
#'@export list2mat.mc
#'@docType methods
#'@title Convert list to matrix
#'@description Convert list to matrix
#'@param list dataframe
#'@return a matrix
#'@author Hoai Tuong Nguyen
list2mat.mc<-function(list,nodes=NULL){
  if (is.null(nodes)){
    nodes<-union(unique(as.character(list[,1])),unique(as.character(list[,2])))
    nodes<-nodes[order(nodes)]
  }
  cooc.mat<-sapply(1:length(nodes), function(i)
    sapply(1:length(nodes), function(j){
      noc<-mean(as.numeric(list[list[,1]==nodes[i]&list[,2]==nodes[j],3]))
      list(ifelse(is.na(noc),0,noc))
    }
    )
  )
  cooc.mat<-t(cooc.mat)
  rownames(cooc.mat)<-colnames(cooc.mat)<-nodes
  
  
  return(cooc.mat)
}




normalizeList.mc<-function(list){
  df<-data.frame(list,stringsAsFactors=FALSE)
  colnames(df)=c("f","t", "weight" ,"method" ,"param" , "dbfile")
  list<-sqldf('SELECT "f","t", AVG(weight) as "weight","method" ,"param" , "dbfile" FROM df GROUP BY "f","t","method" ,"param" , "dbfile" having AVG(weight) >0.5')
  print("after normalized")
  print(list[,1:3])
  return(list)
}

normalizeList.simple.mc<-function(list,alpha=0.5){
  df<-data.frame(list,stringsAsFactors=FALSE)
  colnames(df)=c("f","t", "weight")
  list<-sqldf(sprintf('SELECT "f","t", AVG(weight) as "weight" FROM df GROUP BY "f","t" having AVG(weight) >%f',alpha))
  list<-list[order(list[,3],decreasing=TRUE),]
  return(list)
}


#'@name impute.mc
#'@aliases impute.mc
#'@export impute.mc
#'@docType methods
#'@title Imputing missing data
#'@description Imputing missing data
#'@param data dataframe
#'@param m  number of imputations
#'@author Hoai Tuong Nguyen

impute.mc<-function(data,m){
  library.mc("Amelia",repos="http://r.iq.harvard.edu")
  idvars<-which(sapply(1:ncol(data),function(x) (sum(is.na(data[,x]))==0 | nlevels(as.factor(data[,x]))==1)))
  print(idvars)
  data.imputed<-amelia(data,idvars=idvars,m=m)
  return(data.imputed$imputations[m])
}



#'@aliases getIgraph.mc
#'@export getIgraph.mc
#'@docType methods
#'@title Getting igraph object
#'@description Getting igraph object
#'@param nodes list of nodes
#'@param edges  list of edges
#'@param ncolors  colors of nodes
#'@param ecolors  colors of edges
#'@param eweights  weights of edges
#'@author Hoai Tuong Nguyen
getIgraph.mc<-function(nodes,edges,ncolors,ecolors,eweights,directed=TRUE,layout=layout.kamada.kawai){
  library("igraph")
  g<-graph.data.frame(edges, directed=directed, vertices=cbind(nodes,nodes))  
  V(g)$label<-as.character(nodes)
  #if(exists("ncolors"))
  V(g)$color <- rainbow(length(nodes))[edge.betweenness.community(g)$membership*4]
  #if(exists("ecolors"))
  #  E(g)$color <- ecolors
  #if(exists("eweights"))
  #  E(g)$weight <-eweights
  
  g$layout <- layout 
  
  #layout.fruchterman.reingold
  #layout.kamada.kawai
  #layout.fruchterman.reingold.grid
  #layout.lgl
  #layout.graphopt
  
  
  return(g)  
}


#'@name igraph2gephi.mc
#'@aliases igraph2gephi.mc
#'@export igraph2gephi.mc
#'@docType methods
#'@title Convertin igraph object to Gephi object
#'@description Convertin igraph object to Gephi object
#'@param g igraph object
#'@param filepath  output file
#'@return Gephi object
#'@author Hoai Tuong Nguyen
igraph2gephi.mc <- function(g, filepath="converted_graph.gexf")
{
  require(igraph)
  require(rgexf)
  
  # gexf nodes require two column data frame (id, label)
  # check if the input vertices has label already present
  # if not, just have the ids themselves as the label
  #print("OK1")
  #if(is.null(V(g)$label))
  #  V(g)$label <- as.character(V(g))
  
  # similarily if edges does not have weight, add default 1 weight
  if(is.null(E(g)$weight))
    E(g)$weight <- rep.int(1, ecount(g))
  
  nodes <- data.frame(cbind(V(g)$label, V(g)$label),stringsAsFactors=FALSE)
  edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))
  #print("OK")
  # combine all node attributes into a matrix (and take care of & for xml)
  #vAttrNames <- setdiff(list.vertex.attributes(g), "label") 
  nodesAtt <- data.frame(sapply(list.vertex.attributes(g), function(attr) sub("&", "&#038;",get.vertex.attribute(g, attr))))
  #print("OK3")
  # combine all edge attributes into a matrix (and take care of & for xml)
  eAttrNames <- setdiff(list.edge.attributes(g), "weight") 
  #print("OK4")
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&#038;",get.edge.attribute(g, attr))))
  #print("OK5")
  gAttrNames <- setdiff(list.graph.attributes(g), "layout") 
  # combine all graph attributes into a meta-data
  graphAtt <- data.frame(sapply(gAttrNames, function(attr) sub("&", "&#038;",get.graph.attribute(g, attr))))
  #print(graphAtt)
  # generate the gexf object
  output <- write.gexf(nodes, edges, 
                       edgesWeight=E(g)$weight,
                       edgesAtt = edgesAtt,
                       nodesAtt = nodesAtt,
                       meta=c(list(creator="Gopalakrishna Palem", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
  
  #print(output, filepath, replace=T)
}



#'@aliases updateKDB.mc
#'@export updateKDB.mc
#'@docType methods
#'@title Updating knowledge database
#'@description Updating knowledge database
#'@param arcs list of edges
#'@param method  method used to get edges
#'@param param  parameters of method
#'@param dbfiles  path storing data files
#'@author Hoai Tuong Nguyen
updateKDB.mc<-function(newKDB.or,nodes=NULL){
  KDB.mat.or<-getKDB.mat.mc()
  KDB.or<-mat2list.mc(KDB.mat.or)
  KDB.mat.or<-as.numeric.mat.mc(list2mat.mc(KDB.or,nodes),".")
  
  newKDB.mat.or<-as.numeric.mat.mc(list2mat.mc(newKDB.or,nodes),".")
  
  meanKDB.mat.or<-(KDB.mat.or+newKDB.mat.or)/2
  save(meanKDB.mat.or,file=sprintf("%s/INTEGRATION/KDB.mat.RData",input.dir)) 
  return(meanKDB.mat.or)
}

#'@aliases getKDB.mc
#'@export getKDB.mc
#'@docType methods
#'@title Getting knowledge database
#'@description Getting knowledge database
#'@author Hoai Tuong Nguyen
getKDB.mc<-function(){
  KDB<-get(load(file=sprintf("%s/INTEGRATION/KDB.RData",input.dir)))
}

#'@aliases saveKDB.mc
#'@export saveKDB.mc
#'@docType methods
#'@title Saving knowledge database
#'@description Getting knowledge database
#'@author Hoai Tuong Nguyen
saveKDB.mc<-function(){
  save(KDB,file=sprintf("%s/INTEGRATION/KDB.RData",input.dir))
}

#'@aliases getKDB.mat.mc
#'@export getKDB.mat.mc
#'@docType methods
#'@title Getting matrix of knowledge database
#'@description Getting matrix of knowledge database
#'@author Hoai Tuong Nguyen
getKDB.mat.mc<-function(){
  return(get(load(file=sprintf("%s/INTEGRATION/KDB.mat.RData",input.dir))))
}





#'@aliases neo4j.query.mc
#'@export neo4j.query.mc
#'@docType methods
#'@title Executing a query of Cipher language
#'@description Executing a query of Cipher language
#'@param query
#'@author Hoai Tuong Nguyen
neo4j.query.mc <- function(querystring) {
  library('bitops')
  library('RCurl')
  #sudo apt-get install libcurl4-openssl-dev
  library('RJSONIO')  
  h = basicTextGatherer()
  curlPerform(url = "localhost:7474/db/data/ext/CypherPlugin/graphdb/execute_query", 
              postfields = paste("query", curlEscape(querystring), 
                                 sep = "="), writefunction = h$update, verbose = FALSE)
  result <- fromJSON(h$value())
  data <- data.frame(t(sapply(result$data, unlist)))
  names(data) <- result$columns
  return(data)
}





#'@aliases errplot.mc
#'@export errplot.mc
#'@docType methods
#'@title Plotting scatterplot with error bar
#'@description Plotting scatterplot with error bar
#'@author Hoai Tuong Nguyen
errplot.mc<-function(x,y,err,xlab="x",ylab="y",title=""){
  library("Hmisc")
  df <- data.frame(x,y)
  plot(df,xlab=xlab,ylab=ylab,main=title)
  errbar( x, y, y + err, y - err )
  lines(df)
}

#'@aliases corplot.mc
#'@export corplot.mc
#'@docType methods
#'@title Plotting scatterplot with error bar
#'@description Plotting scatterplot with error bar
#'@author Hoai Tuong Nguyen
corplot.mc<-function(df){
  dat.scaled<-scale(df,center=TRUE,scale=TRUE);
  corrplot(cor(dat.scaled), order = "hclust")
}

#'@aliases lsplot.mc
#'@export lsplot.mc
#'@docType methods
#'@title Plotting linear separator
#'@description Plotting linear separator
#'@author Hoai Tuong Nguyen
lsplot.mc<-function(x,y,labels,xlab="X",ylab="Y",title=""){
  cols<-rep("blue",length(x))
  cols[y>=fitted(lm(y~x))]<-"red"
  plot(x,y,xlab=xlab,ylab=ylab,col=cols, pch=20,
       main=title)
  abline(lm(y~x))
  text(x,y,labels=labels,pos=3,cex=.5)  
}
