#Package: mc
#Title: Omics Data Analysis
#Version: 1.0
#Date: 2014-07-25
#Author: Hoai Tuong Nguyen, David Dernoncourt
#Maintainer: Hoai Tuong Nguyen <hoai-tuong.nguyen@inserm.fr>, David Dernoncourt <me@daviddernoncourt.com>
#Description: Statistical and datamining tools for omics data analysis.
#License: GPL

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
library.mc<-function(pkg,repos="bioc"){
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
library.mc("igraph")
library.mc("RCytoscape")

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



#'@aliases weighted.var.mc
#'@export weighted.var.mc
#'@docType methods
#'@title Calculate weighted variance
#'@description Calculate weighted variance, original function posted by Dr. Gavin Simpson: https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
#'@param x vector of values
#'@param w weights of each values
#'@param na.rm boolean, whether or not remove NAs
#'@author David Dernoncourt
weighted.var.mc = function(x,w=NA,na.rm=FALSE) {
  if(is.na(w[1])) {w=rep(1,length(x));}
  if(na.rm) {
    w = w[i <- !is.na(x)];
    x = x[i];
  }
  sum.w = sum(w);
  return(sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2));
}


#'@aliases stability.CW.mc
#'@export stability.CW.mc
#'@docType methods
#'@title Calculate weighted consistency stability measure
#'@description Calculate weighted consistency stability measure and various variants. Based on Somol 2010 "Evaluating Stability and Comparing Output of Feature Selectors that Optimize Feature Subset Cardinality" (doi:10.1109/tpami.2010.34)
#'@param S selection matrix, 1 row  = 1 selection, 1 column = 1 variable. Value 0 if not selected, 1 if selected.
#'@param type of stability measure. Cs: consistency. CW: weighted consistency. CWrel: relative weigthed consistency.
#'@author David Dernoncourt
stability.CW.mc = function(S,type='CWrel'){
  # S: matrix containing one row per trained FS, with as many columns as features. 1 = feature is selected in this run, 0 = feature isn't selected
  # type: what kind of score we want
  validType=match(type,c('Cf','Cs','CW','CWrel','CWPR')); # http://stackoverflow.com/questions/1169248/r-function-for-testing-if-a-vector-contains-a-given-element
  if(is.na(validType) || length(validType)!=1) {
    type='CWrel';
    warning('Invalid type provided. Defaulting to CWrel.');
  }
  
  freqs=apply(S,2,sum); # get number of occurences of every feature
  n=nrow(S); # get maximum occurence (= number of rows)
  cardY=ncol(S); # get number of features (= number of cols)
  
  # compute C(f) for all features
  Cf=(freqs-1)/(n-1);
  Cf[freqs==0]=0; # fixes the negative consistensy of features never selected
  if(type=='Cf') {return(Cf);}
  
  if(type=='Cs'){
    X=rep(0,cardY); X[freqs>0]=1; # subset of Y representing all features that appear anywhere in S
    cardX=sum(X); # cardinality of X
    Cs=1/cardX*sum(Cf); # that's C(S) already! :)
    return(Cs);
  }
  
  N=sum(S); # total number of occurences of any feature in S
  if(type=='CW'){
    CW=sum(freqs/N*(freqs-1)/(n-1));
    return(CW);
  }
  
  D=N%%cardY; # D = N mod |Y|
  H=N%%n;
  
  if(type=='CWrel'){
    CWrel=(cardY*(N-D+sum(freqs*(freqs-1)))-N^2+D^2) / (cardY*(H^2+n*(N-H)-D)-N^2+D^2);
    return(CWrel);
  }
  
  CWPR=(cardY*(N-D+sum(freqs*(freqs-1)))-N^2+D^2) / (cardY*(N*(n-1)+N-D)-N^2+D^2);
  return(CWPR);
}


#'@aliases artif.data.getDeltaForBE.mc
#'@export artif.data.getDeltaForBE.mc
#'@docType methods
#'@title Calculate coefficient to obtain target Bayes error
#'@description Calculate a delta coefficient so as to obtain target Bayes error when using delta*mu to build artificial data (assuming variance=1)
#'@param wanted_error Bayes error we want (proportion, between 0 and 1)
#'@param mus A vector of means for separation (one class would have delta*mus, the other -delta*mus)
#'@author David Dernoncourt
artif.data.getDeltaForBE.mc = function(wanted_error,mus){
  targetZ=qnorm(1-wanted_error);
  delta=targetZ/sqrt(t(mus)%*%mus);
  return(as.double(delta));
}

#'@aliases artif.data.getBE.mc
#'@export artif.data.getBE.mc
#'@docType methods
#'@title Calculate Bayes error from a vector of means
#'@description Calculate Bayes error obtained when using vector of means delta*mu to build artificial data (assuming variance=1)
#'@param mus A vector of means for separation (one class would have mus, the other -mus)
#'@author David Dernoncourt
artif.data.getBE.mc = function(mus){
  return(as.double(1-pnorm(1*sqrt(t(mus)%*%mus))));
}

#'@aliases artif.data.gen.mu.mc
#'@export artif.data.gen.mu.mc
#'@docType methods
#'@title Generate a vector of means
#'@description Generate a vector of means appropriate to build artificial data
#'@param m How many non-zero variables
#'@param style Distribution for drawing the mus (beware, some of those are weird)
#'@param BayesError Bayes error we want (proportion, between 0 and 1). If you leave NA, no adjustment will be made and you can end up with quite extreme mu values
#'@param styleParam Optional parameter for style. Exponent for style exponential and linear. Useless for unique and uniform.
#'@param D Number of variables. Make >m if you want some variables with ?==0
#'@author David Dernoncourt
artif.data.gen.mu.mc = function (m,style="linear",BayesError=NA,styleParam=2,D=m){
  validStyle=match(style,c("unique","uniform","exponential","linear"));
  if(is.na(validStyle) || length(validStyle)!=1) stop("Invalid style");
  
  if(style=="unique") {mu=rep(1,m);}
  else if(style=="uniform") {mu=1:m;}
  else if(style=="exponential") {
    mu=rexp(m,1);mu=mu^styleParam; # random mus from exponential distribution # curve(dexp(x,1), 0, 5,xlab="?",ylab="Proba density")
    if(styleParam<1){warning("styleParam should be >=1 when using exponential distrib")}
  }
  else if(style=="linear") {
    mu=sample(0:1000000,size=m,replace=TRUE,prob=seq(1,0,length=1000001));mu=(mu/1000000)^styleParam;
    if(styleParam<1){warning("styleParam should be >=1 when using linear distrib")}
  }
  
  mu=sort(mu,decreasing=TRUE); # always return sorted mu, this is required by lots of things so disabling this will break stuff
  
  if(!is.na(BayesError)) {mu=mu*artif.data.getDeltaForBE.mc(BayesError,mu);}
  
  if(D>m) {mu=c(mu,rep(0,D-m));}
  
  return(mu);
}

#'@aliases artif.data.gen.model1.mc
#'@export artif.data.gen.model1.mc
#'@docType methods
#'@title Generate simple artificial data
#'@description Generate simple artificial data: model 1: in one class ?, in the other class -?
#'@param N Number of observations (rows)
#'@param D Number of variables (columns)
#'@param mu ? for C1, -? for C2. NB: can (should!) be a vector, typically the output of artif.data.gen.mu.mc().
#'@param m Number of relevant variables (in case you feel like trimming mu here)
#'@param sigma Standard deviation. You'll probably want to leave this to 1 if you care about the Bayes Error of artif.data.gen.mu.mc()
#'@author David Dernoncourt
artif.data.gen.model1.mc = function(N=8,D=5,mu=1,m=length(mu),sigma=1){
  # prevents a mu longer that doesn't match number of variables (unless mu is a single value)
  if(length(mu)>m){length(mu)=m; cat("Warning: mu was stripped because it was longer than m");}
  else if(length(mu)<m && length(mu)>1){stop("mu length > 1 and doesn't match m yet is smaller");}
  data = matrix(rnorm((N*D),mean=0,sd=sigma), nr=N, nc=D);
  
  C1n=round(N/2); # number of observations in class C1
  #C1vars = matrix(rnorm(C1n*m,mean=mu,sd=sigma), nr=C1n, nc=m);
  C1vars = t(matrix(rnorm(C1n*m,mean=mu,sd=sigma), nr=m, nc=C1n));
  
  C2n=N-C1n; # number of observations in class C2
  #C2vars = matrix(rnorm(C2n*m,mean=-mu,sd=sigma), nr=C2n, nc=m);
  C2vars = t(matrix(rnorm(C2n*m,mean=-mu,sd=sigma), nr=m, nc=C2n));
  
  class = c(rep(1,C1n),rep(0,C2n));
  
  data[1:C1n,1:m]=C1vars;
  data[(C1n+1):N,1:m]=C2vars;
  
  return(list(X=data,y=class));
}

#'@aliases artif.data.gen.model4.mc
#'@export artif.data.gen.model4.mc
#'@docType methods
#'@title Generate simple artificial data with correlation
#'@description Generate simple artificial data with correlation: model 4: like model 1 but with blocks of correlated variables
#'@param N Number of observations (rows)
#'@param D Number of variables (columns)
#'@param mu ? for C1, -? for C2. NB: can (should!) be a vector, typically the output of artif.data.gen.mu.mc().
#'@param m Number of relevant variables (in case you feel like trimming mu here)
#'@param sigma Standard deviation.
#'@param corBlockSize Size of correlated variables blocks
#'@param corStrength Correlation between correlated variables (within blocks)
#'@author David Dernoncourt
artif.data.gen.model4.mc = function(N=8,D=5,mu=1,m=length(mu),sigma=1,corBlockSize=2,corStrength=0.5){
  nBlocks=floor(D/corBlockSize);
  corMatrix=diag(c(rep(1-corStrength,floor(D/corBlockSize)*corBlockSize),rep(1.0,D-floor(D/corBlockSize)*corBlockSize)));
  
  data = artif.data.gen.model1.mc(N=N,D=D,m=m,sigma=sigma,mu=mu);
  
  for(i in 0:(nBlocks-1)) {
    corMatrix[(i*corBlockSize+1):((i+1)*corBlockSize),(i*corBlockSize+1):((i+1)*corBlockSize)]=
      corMatrix[(i*corBlockSize+1):((i+1)*corBlockSize),(i*corBlockSize+1):((i+1)*corBlockSize)]+
      matrix(corStrength,nrow=corBlockSize,ncol=corBlockSize);
  }
  cholMat=chol(corMatrix);
  
  data$X=data$X%*%cholMat;
  
  return(data);
}


#'@aliases beeboxplot.mc
#'@export beeboxplot.mc
#'@docType methods
#'@title Plotting barplot with data points
#'@description Plotting barplot with data points
#'@author Hoai Tuong Nguyen
library.mc("beeswarm")
beeboxplot.mc<-function(x,class,xlab="",ylab="",main="",col=4,pch=16){
  boxplot(x ~ class, 
          outline = TRUE,     ## avoid double-plotting outliers, if any
          main = main,
          xlab=xlab,
          ylab=ylab)
  beeswarm(x ~ class, 
           col = col, pch = pch, add = TRUE)
}






#'@aliases getMultilayerNet.mc
#'@export getMultilayerNet.mc
#'@docType methods
#'@title Creating multilayer networks
#'@description Creating multilayer networks
#'@author Hoai Tuong Nguyen
getMultilayerNet.mc<-function(gl,offset.x,offset.y){
  
  g<-gl[which.max(sapply(1:length(gl), function(x) length(unlist(edgeL(gl[[x]])))))][[1]]
  
  #  create a CytoscapeWindow, after first making sure that no prior window of the same name
  #  name exists already.  (CytoscapeConnections are cheap; create them whenever you need them.)
  cy = CytoscapeConnection()
  window.title = 'vig1'
  if (window.title %in% as.character (getWindowList (cy)))
    deleteWindow (cy, window.title)
  cw = new.CytoscapeWindow (window.title, g)
  displayGraph (cw)
  #layoutNetwork(cw, 'jgraph-spring')
  layoutNetwork(cw, 'force-directed')
  
  pos<-getNodePosition (cw, nodes(g))
  
  
  g<-gl[[1]]
  
  nodes.list<-list()    
  nodes.list[[1]]<-nodes(g)
  
  for (i in 2:length(gl)){
    #g.list[i]<-gl[[i]]
    nodes.list[[i]]<-paste(nodes.list[[1]],"_",i,sep="")
    g<-addNode(nodes.list[[i]],g)
    g<-addEdge(nodes.list[[i-1]], nodes.list[[i]], g,c(rep(0,length(nodes.list[[i]]))))  
    igi <- igraph.from.graphNEL(gl[[i]])
    eli<-get.edgelist(igi)
    g<- addEdge(paste(eli[,1],"_",i,sep=""), paste(eli[,2],"_",i,sep=""), g)  
  }
  
  
  
  
  #   g.i<-igraph.from.graphNEL(g)
  # 
  #   V(g.i)$shape<-c(V(g1.i)$shape,V(g2.i)$shape,V(g3.i)$shape)  
  #   V(g.i)$color<-c(V(g1.i)$color,V(g2.i)$color,V(g3.i)$color)  
  #   E(g.i)$label<-c(E(g1.i)$label,rep("",length(nodes1)),rep("",length(nodes1)),E(g2.i)$label,E(g3.i)$label)
  #   E(g.i)$weight<-c(E(g1.i)$weight,rep(0,length(nodes1)),rep(0,length(nodes1)),E(g2.i)$weight,E(g3.i)$weight)
  #   
  #   
  #   g<-igraph.to.graphNEL(g.i)
  #   
  #   g <- initNodeAttribute (g, "shape", "char", "circle")
  #   g <- initNodeAttribute (g, "color", "char", "black")
  #   
  #   g = initEdgeAttribute (g, "weight", "numeric", 1.0)
  #   g = initEdgeAttribute (g, "size", "numeric", 1.0)
  #   g = initEdgeAttribute (g, "label", "char", "")
  #   
  
  window.title = 'vig2'
  if (window.title %in% as.character (getWindowList (cy)))
    deleteWindow (cy, window.title)
  cw2 = new.CytoscapeWindow (window.title, g)
  displayGraph (cw2)
  layoutNetwork (cw2, 'jgraph-spring')
  
  
  for (i in 1:length(gl)){
    setNodePosition(cw2, nodes.list[[i]], unlist(pos)[grep(".x",names(unlist(pos)))]+offset.x*(i-1), 
                    unlist(pos)[grep(".y",names(unlist(pos)))]+offset.y*ifelse(i%%2==0,0,2))    
  }  
  
  
  
}



#'@name corrplot.mc
#'@aliases corrplot.mc
#'@export corrplot.mc
#'@docType methods
#'@title A visualization of a correlation matrix 
#'@description A graphical display of a correlation matrix, confidence interval.
#'@author Hoai Tuong Nguyen (adapted from "corrplot" package)
corrplot.mc<-function (corr, method = c("circle", "square", "ellipse", "number", 
                                        "shade", "color", "pie"), type = c("full", "lower", "upper"), 
                       add = FALSE, col = NULL, bg = "white", title = "", is.corr = TRUE, 
                       diag = TRUE, outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL, 
                       addCoef.col = NULL, addCoefasPercent = FALSE, order = c("original", 
                                                                               "AOE", "FPC", "hclust", "alphabet"), hclust.method = c("complete", 
                                                                                                                                      "ward", "single", "average", "mcquitty", "median", "centroid"), 
                       addrect = NULL, rect.col = "black", rect.lwd = 2, tl.pos = NULL, 
                       tl.cex = 1, tl.col = "red", tl.offset = 0.4, tl.srt = 90, 
                       cl.pos = NULL, cl.lim = NULL, cl.length = NULL, cl.cex = 0.8, 
                       cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5, addshade = c("negative", 
                                                                                           "positive", "all"), shade.lwd = 1, shade.col = "white", 
                       p.mat = NULL, sig.level = 0.05, sig = c("p-value"),insig = c("pch", "p-value", 
                                                                                    "blank", "n"), p.cex = 1, pch = 4, pch.col = "black", pch.cex = 3, 
                       plotCI = c("n", "square", "circle", "rect"), lowCI.mat = NULL, 
                       uppCI.mat = NULL, ...) 
{
  method <- match.arg(method)
  type <- match.arg(type)
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)
  plotCI <- match.arg(plotCI)
  insig <- match.arg(insig)
  if (!is.matrix(corr) & !is.data.frame(corr)) 
    stop("Need a matrix or data frame!")
  if (is.null(addgrid.col)) {
    addgrid.col <- ifelse(method == "color" | method == "shade", 
                          "white", "grey")
  }
  if (any(corr < cl.lim[1]) | any(corr > cl.lim[2])) 
    stop("color limits should cover matrix")
  if (is.null(cl.lim)) {
    if (is.corr) 
      cl.lim <- c(-1, 1)
    if (!is.corr) 
      cl.lim <- c(min(corr), max(corr))
  }
  intercept <- 0
  zoom <- 1
  if (!is.corr) {
    if (max(corr) * min(corr) < 0) {
      intercept <- 0
      zoom <- 1/max(abs(cl.lim))
    }
    if (min(corr) >= 0) {
      intercept <- -cl.lim[1]
      zoom <- 1/(diff(cl.lim))
    }
    if (max(corr) <= 0) {
      intercept <- -cl.lim[2]
      zoom <- 1/(diff(cl.lim))
    }
    corr <- (intercept + corr) * zoom
  }
  cl.lim2 <- (intercept + cl.lim) * zoom
  int <- intercept * zoom
  if (min(corr) < -1 - .Machine$double.eps || max(corr) > 1 + 
        .Machine$double.eps) {
    stop("The matrix is not in [-1, 1]!")
  }
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
  }
  n <- nrow(corr)
  m <- ncol(corr)
  min.nm <- min(n, m)
  ord <- 1:min.nm
  if (!order == "original") {
    ord <- corrMatOrder(corr, order = order, hclust.method = hclust.method)
    corr <- corr[ord, ord]
  }
  if (is.null(rownames(corr))) 
    rownames(corr) <- 1:n
  if (is.null(colnames(corr))) 
    colnames(corr) <- 1:m
  getPos.Dat <- function(mat) {
    x <- matrix(1:n * m, n, m)
    tmp <- mat
    if (type == "upper") 
      tmp[row(x) > col(x)] <- Inf
    if (type == "lower") 
      tmp[row(x) < col(x)] <- Inf
    if (type == "full") 
      tmp <- tmp
    if (!diag) 
      diag(tmp) <- Inf
    Dat <- tmp[is.finite(tmp)]
    ind <- which(is.finite(tmp), arr.ind = TRUE)
    Pos <- ind
    Pos[, 1] <- ind[, 2]
    Pos[, 2] <- -ind[, 1] + 1 + n
    return(list(Pos, Dat))
  }
  Pos <- getPos.Dat(corr)[[1]]
  n2 <- max(Pos[, 2])
  n1 <- min(Pos[, 2])
  nn <- n2 - n1
  newrownames <- as.character(rownames(corr)[(n + 1 - n2):(n + 
                                                             1 - n1)])
  m2 <- max(Pos[, 1])
  m1 <- min(Pos[, 1])
  mm <- m2 - m1
  newcolnames <- as.character(colnames(corr)[m1:m2])
  DAT <- getPos.Dat(corr)[[2]]
  len.DAT <- length(DAT)
  assign.color <- function(DAT) {
    newcorr <- (DAT + 1)/2
    newcorr[newcorr == 1] <- 1 - 1e-10
    col.fill <- col[floor(newcorr * length(col)) + 1]
  }
  col.fill <- assign.color(DAT)
  isFALSE = function(x) identical(x, FALSE)
  isTRUE = function(x) identical(x, TRUE)
  if (isFALSE(tl.pos)) {
    tl.pos <- "n"
  }
  if (is.null(tl.pos) | isTRUE(tl.pos)) {
    if (type == "full") 
      tl.pos <- "lt"
    if (type == "lower") 
      tl.pos <- "ld"
    if (type == "upper") 
      tl.pos <- "td"
  }
  if (isFALSE(cl.pos)) {
    cl.pos <- "n"
  }
  if (is.null(cl.pos) | isTRUE(cl.pos)) {
    if (type == "full") 
      cl.pos <- "r"
    if (type == "lower") 
      cl.pos <- "b"
    if (type == "upper") 
      cl.pos <- "r"
  }
  if (outline) 
    col.border <- "black"
  if (!outline) 
    col.border <- col.fill
  if (!add) {
    par(mar = mar, bg = "white")
    plot.new()
    xlabwidth <- ylabwidth <- 0
    for (i in 1:50) {
      xlim <- c(m1 - 0.5 - xlabwidth, m2 + 0.5 + mm * cl.ratio * 
                  (cl.pos == "r"))
      ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"), 
                n2 + 0.5 + ylabwidth)
      plot.window(xlim + c(-0.2, 0.2), ylim + c(-0.2, 0.2), 
                  asp = 1, xaxs = "i", yaxs = "i")
      x.tmp <- max(strwidth(newrownames, cex = tl.cex))
      y.tmp <- max(strwidth(newcolnames, cex = tl.cex))
      if (min(x.tmp - xlabwidth, y.tmp - ylabwidth) < 1e-04) 
        break
      xlabwidth <- x.tmp
      ylabwidth <- y.tmp
    }
    if (tl.pos == "n" | tl.pos == "d") 
      xlabwidth <- ylabwidth <- 0
    if (tl.pos == "td") 
      ylabwidth <- 0
    if (tl.pos == "ld") 
      xlabwidth <- 0
    laboffset <- strwidth("W", cex = tl.cex) * tl.offset
    xlim <- c(m1 - 0.5 - xlabwidth - laboffset, m2 + 0.5 + 
                mm * cl.ratio * (cl.pos == "r")) + c(-0.35, 0.15)
    ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"), 
              n2 + 0.5 + ylabwidth * abs(sin(tl.srt * pi/180)) + 
                laboffset) + c(-0.15, 0.35)
    if (.Platform$OS.type == "windows") {
      windows.options(width = 7, height = 7 * diff(ylim)/diff(xlim))
    }
    plot.window(xlim = xlim, ylim = ylim, asp = 1, xlab = "", 
                ylab = "", xaxs = "i", yaxs = "i")
  }
  laboffset <- strwidth("W", cex = tl.cex) * tl.offset
  symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                         len.DAT), bg = bg, fg = bg)
  if (method == "circle" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, bg = col.fill, 
            circles = 0.9 * abs(DAT)^0.5/2, fg = col.border)
  }
  if (method == "ellipse" & plotCI == "n") {
    ell.dat <- function(rho, length = 99) {
      k <- seq(0, 2 * pi, length = length)
      x <- cos(k + acos(rho)/2)/2
      y <- cos(k - acos(rho)/2)/2
      return(cbind(rbind(x, y), c(NA, NA)))
    }
    ELL.dat <- lapply(DAT, ell.dat)
    ELL.dat2 <- 0.85 * matrix(unlist(ELL.dat), ncol = 2, 
                              byrow = TRUE)
    ELL.dat2 <- ELL.dat2 + Pos[rep(1:length(DAT), each = 100), 
                               ]
    polygon(ELL.dat2, border = col.border, col = col.fill)
  }
  if (method == "number" & plotCI == "n") {
    text(Pos[, 1], Pos[, 2], font = 2, col = col.fill, labels = round((DAT - 
                                                                         int) * ifelse(addCoefasPercent, 100, 1)/zoom, ifelse(addCoefasPercent, 
                                                                                                                              0, 2)))
  }
  if (method == "pie" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, 
                                                           len.DAT) * 0.85)
    pie.dat <- function(theta, length = 100) {
      k <- seq(pi/2, pi/2 - theta, length = 0.5 * length * 
                 abs(theta)/pi)
      x <- c(0, cos(k)/2, 0)
      y <- c(0, sin(k)/2, 0)
      return(cbind(rbind(x, y), c(NA, NA)))
    }
    PIE.dat <- lapply(DAT * 2 * pi, pie.dat)
    len.pie <- unlist(lapply(PIE.dat, length))/2
    PIE.dat2 <- 0.85 * matrix(unlist(PIE.dat), ncol = 2, 
                              byrow = TRUE)
    PIE.dat2 <- PIE.dat2 + Pos[rep(1:length(DAT), len.pie), 
                               ]
    polygon(PIE.dat2, border = "black", col = col.fill)
  }
  if (method == "shade" & plotCI == "n") {
    addshade <- match.arg(addshade)
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = addgrid.col)
    shade.dat <- function(w) {
      x <- w[1]
      y <- w[2]
      rho <- w[3]
      x1 <- x - 0.5
      x2 <- x + 0.5
      y1 <- y - 0.5
      y2 <- y + 0.5
      dat <- NA
      if ((addshade == "positive" || addshade == "all") & 
            rho > 0) {
        dat <- cbind(c(x1, x1, x), c(y, y1, y1), c(x, 
                                                   x2, x2), c(y2, y2, y))
      }
      if ((addshade == "negative" || addshade == "all") & 
            rho < 0) {
        dat <- cbind(c(x1, x1, x), c(y, y2, y2), c(x, 
                                                   x2, x2), c(y1, y1, y))
      }
      return(t(dat))
    }
    pos_corr <- rbind(cbind(Pos, DAT))
    pos_corr2 <- split(pos_corr, 1:nrow(pos_corr))
    SHADE.dat <- matrix(na.omit(unlist(lapply(pos_corr2, 
                                              shade.dat))), byrow = TRUE, ncol = 4)
    segments(SHADE.dat[, 1], SHADE.dat[, 2], SHADE.dat[, 
                                                       3], SHADE.dat[, 4], col = shade.col, lwd = shade.lwd)
  }
  if (method == "square" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = abs(DAT)^0.5, 
            bg = col.fill, fg = col.border)
  }
  if (method == "color" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = col.border)
  }
  symbols(Pos, add = TRUE, inches = FALSE, bg = NA, squares = rep(1, 
                                                                  len.DAT), fg = addgrid.col)
  if (plotCI != "n") {
    if (is.null(lowCI.mat) || is.null(uppCI.mat)) 
      stop("Need lowCI.mat and uppCI.mat!")
    if (!order == "original") {
      lowCI.mat <- lowCI.mat[ord, ord]
      uppCI.mat <- uppCI.mat[ord, ord]
    }
    pos.lowNew <- getPos.Dat(lowCI.mat)[[1]]
    lowNew <- getPos.Dat(lowCI.mat)[[2]]
    pos.uppNew <- getPos.Dat(uppCI.mat)[[1]]
    uppNew <- getPos.Dat(uppCI.mat)[[2]]
    if (!(method == "circle" || method == "square")) 
      stop("method shoud be circle or square if draw confidence interval!")
    k1 <- (abs(uppNew) > abs(lowNew))
    bigabs <- uppNew
    bigabs[which(!k1)] <- lowNew[!k1]
    smallabs <- lowNew
    smallabs[which(!k1)] <- uppNew[!k1]
    sig <- sign(uppNew * lowNew)
    if (plotCI == "circle") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(bigabs)^0.5/2, 
              bg = ifelse(sig > 0, col.fill, col[ceiling((bigabs + 
                                                            1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                               col.fill, col[ceiling((bigabs + 1) * length(col)/2)]))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(smallabs)^0.5/2, 
              bg = ifelse(sig > 0, bg, col[ceiling((smallabs + 
                                                      1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                         col.fill, col[ceiling((smallabs + 1) * length(col)/2)]))
    }
    if (plotCI == "square") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(bigabs)^0.5, bg = ifelse(sig > 
                                                                       0, col.fill, col[ceiling((bigabs + 1) * length(col)/2)]), 
              fg = ifelse(sig > 0, col.fill, col[ceiling((bigabs + 
                                                            1) * length(col)/2)]))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(smallabs)^0.5, 
              bg = ifelse(sig > 0, bg, col[ceiling((smallabs + 
                                                      1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                         col.fill, col[ceiling((smallabs + 1) * length(col)/2)]))
    }
    if (plotCI == "rect") {
      rect.width <- 0.25
      rect(pos.uppNew[, 1] - rect.width, pos.uppNew[, 2] + 
             smallabs/2, pos.uppNew[, 1] + rect.width, pos.uppNew[, 
                                                                  2] + bigabs/2, col = col.fill, border = col.fill)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + DAT/2, pos.lowNew[, 1] + rect.width, pos.lowNew[, 
                                                                                                             2] + DAT/2, col = "black", lwd = 1)
      segments(pos.uppNew[, 1] - rect.width, pos.uppNew[, 
                                                        2] + uppNew/2, pos.uppNew[, 1] + rect.width, 
               pos.uppNew[, 2] + uppNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + lowNew/2, pos.lowNew[, 1] + rect.width, 
               pos.lowNew[, 2] + lowNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - 0.5, pos.lowNew[, 2], 
               pos.lowNew[, 1] + 0.5, pos.lowNew[, 2], col = "grey70", 
               lty = 3)
    }
  }
  if (!is.null(p.mat) & !insig == "n") {
    if (!order == "original") 
      p.mat <- p.mat[ord, ord]
    pos.pNew <- getPos.Dat(p.mat)[[1]]
    pNew <- getPos.Dat(p.mat)[[2]]
    ind.p <- which(pNew > (sig.level))
    if (insig == "pch") {
      points(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
             pch = pch, col = pch.col, cex = pch.cex, lwd = 2)
    }
    if (insig == "p-value") {
      text(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
           round(pNew[ind.p], 2), col = pch.col)
    }
    if (insig == "blank") {
      symbols(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
              inches = FALSE, squares = rep(1, length(pos.pNew[, 
                                                               1][ind.p])), fg = addgrid.col, bg = bg, add = TRUE)
    }
    
  }
  
  if (!is.null(p.mat) & !sig == "n") {
    options("scipen"=-100, "digits"=2)
    r.mat <- corr[ord, ord]
    rNew <- getPos.Dat(r.mat)[[2]]    
    if (sig == "p-value") {
      text(pos.pNew[, 1], pos.pNew[, 2], 
           sprintf("r= %s\n p=%s",format(rNew, scientific=TRUE),format(pNew, scientific=TRUE)), col = pch.col,cex = p.cex)
    }
    
    
  }
  
  if (cl.pos != "n") {
    colRange <- assign.color(cl.lim2)
    ind1 <- which(col == colRange[1])
    ind2 <- which(col == colRange[2])
    colbar <- col[ind1:ind2]
    if (is.null(cl.length)) 
      cl.length <- ifelse(length(colbar) > 20, 11, length(colbar) + 
                            1)
    labels <- seq(cl.lim[1], cl.lim[2], length = cl.length)
    at <- seq(0, 1, length = length(labels))
    if (cl.pos == "r") {
      vertical <- TRUE
      xlim <- c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
      ylim <- c(n1 - 0.5, n2 + 0.5)
    }
    if (cl.pos == "b") {
      vertical <- FALSE
      xlim <- c(m1 - 0.5, m2 + 0.5)
      ylim <- c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn * 
                  0.02)
    }
    colorlegend(colbar = colbar, labels = round(labels, 2), 
                offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex, 
                xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }
  if (tl.pos != "n") {
    ylabwidth2 <- strwidth(newrownames, cex = tl.cex)
    xlabwidth2 <- strwidth(newcolnames, cex = tl.cex)
    pos.xlabel <- cbind(m1:m2, n2 + 0.5 + laboffset)
    pos.ylabel <- cbind(m1 - 0.5, n2:n1)
    if (tl.pos == "td") {
      if (type != "upper") 
        stop("type should be \"upper\" if tl.pos is \"dt\".")
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
    }
    if (tl.pos == "ld") {
      if (type != "lower") 
        stop("type should be \"lower\" if tl.pos is \"ld\".")
      pos.xlabel <- cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
    }
    if (tl.pos == "d") {
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
      pos.ylabel <- pos.ylabel[1:min(n, m), ]
      symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE, 
              bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1, 
                                                                       length(pos.ylabel[, 1])))
      text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n, 
                                                                     m)], col = tl.col, cex = tl.cex, ...)
    }
    else {
      text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames, 
           srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5, 
                                                     0), c(0, 0)), col = tl.col, cex = tl.cex, offset = tl.offset, 
           ...)
      text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames, 
           col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset, 
           ...)
    }
  }
  title(title, ...)
  if (!is.null(addCoef.col) & (!method == "number")) {
    text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = round((DAT - 
                                                                  int) * ifelse(addCoefasPercent, 100, 1)/zoom, ifelse(addCoefasPercent, 
                                                                                                                       0, 2)))
  }
  if (type == "full" & plotCI == "n" & !is.null(addgrid.col)) 
    rect(m1 - 0.5, n1 - 0.5, m2 + 0.5, n2 + 0.5, border = addgrid.col)
  if (!is.null(addrect) & order == "hclust" & type == "full") {
    corrRect.hclust(corr, k = addrect, method = hclust.method, 
                    col = rect.col, lwd = rect.lwd)
  }
  invisible(corr)
}



#'@aliases filterCPSS.mc
#'@export filterCPSS.mc
#'@docType methods
#'@title Performs Complementary Pairs Stability Selection (CPSS)
#'@description Performs Complementary Pairs Stability Selection (CPSS), which is a particular case of ensemble FS method
#'@param X A matrix of features, where one row = one sample / one column = one feature
#'@param y A vector of outcomes (coded as 0 or 1)
#'@param subFilter Which filter to use
#'@param subthreshold Threshold for the subFilter to select features
#'@param PSLimit Limit for the probability to be selected. Default 0.5 = we select features which are selected in 50% of internal runs
#'@param nPairs Number of pairs in the ensemble.
#'@param subExtraParam A list containing more parameters for the subFilter (mostly added for comprehensiveness reasons, shouldn't be much useful here)
#'@param verbose Boolean, whether or not to output more details (default TRUE - it's not that much)
#'@author David Dernoncourt

# based on Shah 2011 - Variable selection with error control - another look at stability selection - http://arxiv.org/abs/1105.5578

filterCPSS.mc = function(X,y,subFilter,subThreshold,PSLimit=0.50,nPairs=100,subExtraParam=NA,verbose=TRUE) {
  # nombre d'observations
  N=length(y);
  NHalf1=round(N/2);
  NHalf2=N-NHalf1;
  # matrice pour stocker les nPairs*2 sélections (0: non sélectionné, 1: sélectionné)
  selectionMatrix=matrix(0,nrow=nPairs*2,ncol=ncol(X));
  
  for(i in 1:nPairs) {
    # on crée la paire d'index
    indexes=list(NULL);
    repeat {
      indexes[[1]]=sample(1:N,NHalf1,replace=F);
      indexes[[2]]=(1:N)[-indexes[[1]]];
      if(length(which(y[indexes[[1]]]==0))>0 && length(which(y[indexes[[1]]]==1))>0 && length(which(y[indexes[[2]]]==0))>0 && length(which(y[indexes[[2]]]==1))>0) {
        break;
      }
    }
    for(j in 1:2) {
      filtered=doFilter.mc(X=X[indexes[[j]],],y=y[indexes[[j]]],filter=subFilter,threshold=subThreshold,extraParam=subExtraParam,verbose=verbose);
      selectionMatrix[(i-1)*2+j,filtered$Xselected]=1;
    }
  }
  featureScores=apply(selectionMatrix,2,sum)/(nPairs*2);
  featureRanks=rank(-featureScores);
  Xselected=featureScores>=PSLimit;
  # il faut s'assure que l'output contient au moins 2 variables
  i=0.01;
  while(length(which(Xselected))<2) {
    Xselected=featureScores>=(PSLimit-i);
    i=i+1;
  }
  if(verbose) {
    cat('PSlimit:',PSLimit,'| selection length:',length(which(Xselected)),'out of',length(Xselected));
  }
  return(list(
    featureScores=featureScores,
    featureRanks=featureRanks,
    Xselected=Xselected));
}


#'@aliases filterEnsemble.mc
#'@export filterEnsemble.mc
#'@docType methods
#'@title Performs ensemble feature selection
#'@description Ensemble feature selection, which can use all filters supported by doFilter.mc
#'@param X A matrix of features, where one row = one sample / one column = one feature
#'@param y A vector of outcomes (coded as 0 or 1)
#'@param subFilter Which filter to use. Can be a vector of filters (for hybrid ensemble), but this might then be buggy.
#'@param subthreshold Threshold for the subFilter to select features
#'@param PSLimit Limit for the probability to be selected. Default 0.5 = we select features which are selected in 50% of internal runs. Note that in any case the function also return scores for all features so you can do your own thresholding (like, top 100) outside.
#'@param ensembleRuns Number of resamplings in the ensemble. Default is 20, which is usually a good tradeoff between accuracy/stability and computation time.
#'@param aggreType Type of agggregation: 'stab' (stability selection), 'avgScore' (average score), avgRank, expRank, bestScore, bestRank
#'@param subExtraParam A list containing more parameters for the subFilter (mostly added for comprehensiveness reasons, shouldn't be much useful here)
#'@param verbose Boolean, whether or not to output more details (default TRUE - it's not that much)
#'@author David Dernoncourt
filterEnsemble.mc = function(X,y,subFilter,subthreshold,PSLimit=0.50,ensembleRuns=20,aggreType='avgScore',subExtraParam=NA,verbose=TRUE) {
  # nombre d'observations
  N=length(y);
  # matrice pour stocker les ensembleRuns sélections (0: non sélectionné, 1: sélectionné), scores et rangs
  selectionMatrix=scoreMatrix=rankMatrix=matrix(0,nrow=ensembleRuns,ncol=ncol(X));
  
  for(i in 1:ensembleRuns) {
    thisSubFilter=subFilter[i%%length(subFilter)+1];
    # on choisi les indexes (avec remise)
    indexes=NA;
    # en s'assurant qu'on a au moins un y=0 et un y=1 dans le tirage
    repeat {
      indexes=sample(1:N,N,replace=T);
      if(length(which(y[indexes]==0))>0 && length(which(y[indexes]==1))>0) {
        break;
      }
    }
    # on filtre et on enregistre le résultat
    filtered=doFilter.mc(X=X[indexes,],y=y[indexes],filter=thisSubFilter,threshold=subthreshold,extraParam=subExtraParam,verbose=verbose);
    selectionMatrix[i,filtered$Xselected]=1;
    scoreMatrix[i,]=filtered$featureScores;
    rankMatrix[i,]=filtered$featureRanks;
  }
  maxRank=matrix(max(rankMatrix),nrow=nrow(rankMatrix),ncol=ncol(rankMatrix));
  if(aggreType=='stab') {
    featureScores=apply(selectionMatrix,2,sum)/ensembleRuns;
  } else if(aggreType=='avgScore') {
    featureScores=apply(scoreMatrix,2,sum)/ensembleRuns;
  } else if(aggreType=='bestScore') {
    featureScores=apply(scoreMatrix,2,max);
  } else if(aggreType=='avgRank') {
    rankMatrix=(maxRank-rankMatrix)/maxRank;
    featureScores=apply(rankMatrix,2,sum)/ensembleRuns;
  } else if(aggreType=='bestRank') {
    rankMatrix=(maxRank-rankMatrix)/maxRank;
    featureScores=apply(rankMatrix,2,max);
  } else if(aggreType=='expRank') {
    rankMatrix=(maxRank-rankMatrix)^2/maxRank^2;
    featureScores=apply(rankMatrix,2,sum)/ensembleRuns;
  }
  featureRanks=rank(-featureScores);
  Xselected=featureScores>=PSLimit;
  # il faut s'assurer que l'output contient au moins 2 variables
  i=0.01;
  while(length(which(Xselected))<2) {
    Xselected=featureScores>=(PSLimit-i);
    i=i+0.01;
  }
  if(verbose) {cat("limit: ",PSLimit," | selection length: ",length(which(Xselected)));}
  return(list(
    featureScores=featureScores,
    featureRanks=featureRanks,
    Xselected=Xselected));
}



#'@aliases doFilter.mc
#'@export doFilter.mc
#'@docType methods
#'@title Applies a feature selection method
#'@description Wrapper for various binary feature selection methods
#'@param X A matrix of features, where one row = one sample / one column = one feature
#'@param y A vector of outcomes (coded as 0 or 1)
#'@param filter Which filter to use
#'@param threshold Filter threshold (useful for the few feature selection methods which don't output feature scores)
#'@param extraParam A list containing more filter parameters (mostly used for ensemble feature selection)
#'@param verbose Boolean, whether or not to output more details (default TRUE - it's not that much)
#'@author David Dernoncourt
doFilter.mc = function(X,y,filter,threshold,extraParam=NA,verbose=TRUE,
                       subFilter="ttest", # internal filter for techniques such as CPSS
                       wrapThreshold=0.50, # threshold for the technique such as CPSS
                       ensembleRuns=5, # number of bootstraps if performing ensemble
                       ensembleAggreg='avgScore', # aggregation used for ensemble
                       hybridWeight=NA, # manual weight for hybrid ensemble
                       ...
) {
  if(verbose) {
    cat(filter);
  }
  if(filter=="ttest") {
    library.mc('multtest');
    featureScores = mt.teststat(t(X), y, ...);
    featureRanks = rank(-abs(featureScores));
    Xselected = featureRanks<=threshold;
  } else if(filter=="CMAelasticnet" || filter=="CMAsvmRfe" || filter=="CMArf" || filter=="CMAshrinkcat"
            || filter=="CMAboosting" || filter=="CMAgolub" || filter=="CMAwelch" || filter=="CMAwilcox" || filter=="CMAlasso") {
    library.mc('CMA');
    if(filter=="CMAelasticnet") {
      #library.mc('glmpath');
      selection=GeneSelection(X=X, y=as.factor(y), method='elasticnet', trace=verbose, ...);
    } else if(filter=="CMAsvmRfe") {
      library.mc('e1071');
      selection=GeneSelection(X=X, y=as.factor(y), method='rfe', trace=verbose, ...);
    } else if(filter=="CMArf") {
      library.mc('randomForest');
      selection=GeneSelection(X=X, y=as.factor(y), method='rf', trace=verbose, ...);
    } else if(filter=="CMAshrinkcat") {
      selection=GeneSelection(X=X, y=as.factor(y), method='shrinkcat', trace=verbose);
    } else if(filter=="CMAboosting") {
      selection=GeneSelection(X=X, y=as.factor(y), method='boosting', trace=verbose, ...);
    } else if(filter=="CMAgolub") {
      selection=GeneSelection(X=X, y=as.factor(y), method='golub', trace=verbose);
    } else if(filter=="CMAwelch") {
      selection=GeneSelection(X=X, y=as.factor(y), method='welch.test', trace=verbose);
    } else if(filter=="CMAwilcox") {
      selection=GeneSelection(X=X, y=as.factor(y), method='wilcox.test', trace=verbose);
    } else if(filter=="CMAlasso") {
      selection=GeneSelection(X=X, y=as.factor(y), method='lasso', trace=verbose, ...);
    }
    selection=toplist(selection,k=ncol(X),iter=1,show=FALSE);
    
    featureRanks=featureScores=rep(NA,ncol(X));
    featureScores[selection$index]=selection$importance;
    selection$rank=1:ncol(X);
    featureRanks[selection$index]=selection$rank;
    Xselected = featureRanks<=threshold;
  } else if(filter=="sdaT") {
    library.mc('sda');
    selection = sda.ranking(X, y, diagonal=TRUE, ...);
    featureRanks=featureScores=rep(NA,ncol(X));
    featureScores[selection[1:ncol(X),"idx"]]=selection[1:ncol(X),"score"];
    featureRanks = rank(-featureScores);
    Xselected = featureRanks<=threshold;
  } else if(filter=="MI") { #my manual mutual information
    tmpDiscr=discretize(as.data.frame(X));
    featureScores=rep(NA,ncol(X));
    tmpLoop=1:length(featureScores);
    for(i in tmpLoop){
      featureScores[i]=mutinformation(X=tmpDiscr[,i],Y=y,method="emp")
    }
    featureRanks=rank(-featureScores);
    Xselected = featureRanks<=threshold;
  } else if(filter=="random") { # not a real filter, we just pick variables randomly! (useful to get baseline stability)
    featureScores = runif(ncol(X),0,10);
    featureRanks = rank(-abs(featureScores));
    Xselected = featureRanks<=threshold;
  } else if(filter=="firstD") { # we just pick the first [threshold] variables (useful to get best var in artificial data)
    featureScores = ncol(X):1;
    featureRanks = rank(-abs(featureScores));
    Xselected = featureRanks<=threshold;
  } else if(filter=="ReliefF") {
    library.mc('CORElearn');
    tmpDataFrame=as.data.frame(X);
    tmpDataFrame[,'y']=y;
    tmpObject=attrEval(y~.,tmpDataFrame,estimator='ReliefFequalK', ...);
    featureScores=abs(tmpObject);
    names(featureScores)=NULL;
    featureRanks=rank(-featureScores);
    Xselected = featureRanks<=threshold;
  } else if(filter=="none") {
    featureRanks=featureScores=rep(NA,ncol(X));Xselected = rep(TRUE,ncol(X));
  } else if(filter=="CPSS" || filter=="CPSS2") {
    tmpOut=filterCPSS.mc(X=X,y=y,subFilter=extraParam[['subFilter']],subThreshold=extraParam[['subThreshold']],PSLimit=threshold,nPairs=extraParam[['nPairs']]);
    featureRanks=tmpOut$featureRanks;
    featureScores=tmpOut$featureScores;
    Xselected=tmpOut$Xselected;
    if (filter=='CPSS2') { # dans CPSS2 au lieu de garder les variables tq pSelected>wrapThreshold, on garde les threshold meilleures variables
      Xselected=(featureRanks<=threshold);
    }
  } else if(filter=="Ensemble" || filter=="Ensemble2") {
    tmpOut=filterEnsemble.mc(X=X,y=y,subFilter=extraParam[['subFilter']],subthreshold=extraParam[['subThreshold']],PSLimit=threshold,ensembleRuns=extraParam[['ensembleRuns']],aggreType=extraParam[['ensembleAggreg']]);
    featureRanks=tmpOut$featureRanks;
    featureScores=tmpOut$featureScores;
    Xselected=tmpOut$Xselected;
    if(filter=="Ensemble2") { # dans Ensemble2 au lieu de garder les variables tq pSelected>wrapThreshold, on garde les threshold meilleures variables
      Xselected=(featureRanks<=threshold);
    }
  } else {stop(paste("Filter",filter,"doesn't exist"));}
  

    ## below = todo (except last part = return output! ;))
    ## Tuong (26/08/2014): I removed "else" to debug an error 
    if(filter=="quickWrapper") { # our basic wrapper
      tmpOut = quickWrapper(X=X,y=y,classifier=subFilter,threshold=threshold);
      featureScores = tmpOut$featureScores;
      featureRanks = rank(-abs(featureScores));
      Xselected = rep(FALSE,ncol(X));
      Xselected[tmpOut$selectedIndexes] = TRUE;
    } else if(filter=="HybridEnsemble") {
      tmpOut=filterHybridEnsemble(X=X,y=y,subFilter=subFilter,subthreshold=threshold,PSLimit=wrapThreshold,ensembleRuns=ensembleRuns,aggreType=ensembleAggreg,manualWeight=hybridWeight);
      featureRanks=tmpOut$featureRanks;
      featureScores=tmpOut$featureScores;
      #Xselected=tmpOut$Xselected;
      Xselected=(featureRanks<=threshold);
    } else if(filter=="Boosting") {
      tmpOut=filterBoosting(X=X,y=y,subFilter=NA,stepSize=10,threshold=threshold,aggreType=NA);
      featureRanks=tmpOut$featureRanks;
      featureScores=tmpOut$featureScores;
      Xselected=tmpOut$Xselected;
    } else if(filter=="Boosting2") {
      tmpOut=filterBoosting2(X=X,y=y,subFilter='ttest',subthreshold=10,boostRuns=ensembleRuns,threshold=threshold,aggreType=ensembleAggreg,scoreWeightsOn=TRUE);
      featureRanks=tmpOut$featureRanks;
      featureScores=tmpOut$featureScores;
      Xselected=tmpOut$Xselected;
    }
    #####

  
  return(list(
    featureScores=featureScores,
    featureRanks=featureRanks,
    Xselected=Xselected));
}


#'@aliases doClassify.mc
#'@export doClassify.mc
#'@docType methods
#'@title Train and applies a classification method
#'@description Wrapper for various binary classification methods
#'@param Xtrain A matrix of features for traning, where one row = one sample / one column = one feature
#'@param ytrain A vector of outcomes for traning (coded as 0 or 1)
#'@param Xtest A matrix of features for testing, where one row = one sample / one column = one feature
#'@param classifier Which classifier to use
#'@param probabilities Boolean, whether or not to return probability of class==1 instead of class (default FALSE)
#'@param threshold Filter threshold (useful for the few feature selection methods which don't output feature scores)
#'@param extraParam A list containing more filter parameters (mostly used for ensemble feature selection)
#'@param verbose Boolean, whether or not to output more details (default TRUE - it's not that much)
#'@author David Dernoncourt
doClassify.mc = function(Xtrain,ytrain,Xtest,classifier,probabilities=FALSE,extraParam=list(),verbose=TRUE,
	...
){
  if(classifier=="sdaDDA") {
    # trains
  	sda.fit = sda(as.matrix(Xtrain), ytrain, diagonal=TRUE);
  	# predicts
  	if(probabilities==FALSE) {
  	  ynew = predict(sda.fit, as.matrix(Xtest))$class;
  	} else {
  	  posteriors = predict(sda.fit, as.matrix(Xtest))$posterior;
  	  ynew=posteriors[,2];
  	}
  } else if(classifier=="sdaLDA") {
  	sda.fit = sda(as.matrix(Xtrain), ytrain, diagonal=FALSE);
    if(probabilities==FALSE) {
      ynew = predict(sda.fit, as.matrix(Xtest))$class;
    } else {
      posteriors = predict(sda.fit, as.matrix(Xtest))$posterior;
      ynew=posteriors[,2];
    }
  } else if(classifier=="knn") {
    if(probabilities==FALSE) {
      stop('This kNN cannot output probabilities');
    }
    library.mc('class');
    if(is.null(extraParam[['nNeighbors']])) {extraParam[['nNeighbors']]=3;}
    ynew = knn(Xtrain, Xtest, ytrain, k=extraParam[['nNeighbors']]);
  } else if(classifier=="RF") {
    rf.fit = randomForest(x=as.matrix(Xtrain), y=as.factor(ytrain));
	if(is.null(extraParam[['RFclasswt']])) {extraParam[['RFclasswt']]=NULL;}
    if(probabilities==FALSE) {
  	  ynew = predict(rf.fit, as.matrix(Xtest), classwt=extraParam[['RFclasswt']]);
    } else {
      posteriors = predict(rf.fit, as.matrix(Xtest), type="prob", classwt=extraParam[['RFclasswt']], corr.bias=FALSE);
      ynew=posteriors[,2];
    }
  } else if(classifier=="cforest") {
    library.mc('party');
    rf.fit = cforest(formula=ytrain~., data=as.data.frame(cbind(Xtrain,ytrain)));
    # predict behaves strangely, alway outputs probabilities but as vector if type not set and as messy list if type="prob"...
    posteriors = as.vector(predict(rf.fit, Xtest));
    if(probabilities==FALSE) {
      ynew=rep(0,length(posteriors));
      ynew[which(posteriors)>0.5]=1;
    } else {
      ynew=round(posteriors,3);
    }
    names(ynew)=rownames(Xtest);
  } else if(classifier=="SVM.linear") {
    if(probabilities==FALSE) {
      svm.e1071.mdl = svm(Xtrain, as.factor(ytrain), kernel="linear");
    	ynew = predict(svm.e1071.mdl, Xtest);
    } else {
      svm.e1071.mdl = svm(Xtrain, as.factor(ytrain), kernel="linear", probability=TRUE);
      posteriors = predict(svm.e1071.mdl, Xtest, probability=TRUE);
      ynew=attr(posteriors,"probabilities")[,2];
    }
  } else if(classifier=="SVM.rbf") {
    if(probabilities==FALSE) {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="rbfdot");
      ynew = predict(svm.kernlab.mdl,Xtest);
    } else {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="rbfdot",prob.model=TRUE);
      posteriors = predict(svm.kernlab.mdl, Xtest, type="probabilities");
      ynew=posteriors[,2];
    }
  } else if(classifier=="SVM.anova") {
    if(probabilities==FALSE) {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="anovadot");
      ynew = predict(svm.kernlab.mdl,Xtest);
    } else {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="anovadot",prob.model=TRUE);
      posteriors = predict(svm.kernlab.mdl, Xtest, type="probabilities");
      ynew=posteriors[,2];
    }
  } else if(classifier=="SVM.laplace") {
    if(probabilities==FALSE) {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="laplacedot");
      ynew = predict(svm.kernlab.mdl,Xtest);
    } else {
      svm.kernlab.mdl = ksvm(Xtrain, as.factor(ytrain),kernel="laplacedot",prob.model=TRUE);
      posteriors = predict(svm.kernlab.mdl, Xtest, type="probabilities");
      ynew=posteriors[,2];
    }
  } else if(classifier=="nnet") {
    myNetwork = nnet(Xtrain, ytrain,linout=TRUE,size=round(ncol(Xtrain)/2),trace=FALSE,MaxNWts=20000);cat(".");
    if(probabilities==FALSE) {
  	  ynew = round(predict(myNetwork,Xtest));
    } else {
      ynew = predict(myNetwork,Xtest);
    }
  } else if(classifier=="neuralMLP") {
  	if(probabilities==FALSE) {
      stop('neuralMLP should be used with probabilities');
    }
    library.mc('neural');
	if(is.null(extraParam[['neurons']])) {extraParam[['neurons']]=c(5,3);}
	if(is.null(extraParam[['it']])) {extraParam[['it']]=50;}
	if(is.null(extraParam[['alpha']])) {extraParam[['alpha']]=0.2;}
	if(is.null(extraParam[['thresh']])) {extraParam[['thresh']]=0;}
	if(is.null(extraParam[['actfns']])) {extraParam[['actfns']]=c();}
	tmpNetwork=mlptrain(inp=as.matrix(Xtrain),neurons=extraParam[['neurons']],out=as.matrix(ytrain),alfa=extraParam[['alpha']],it=extraParam[['it']],online=TRUE,permute=TRUE,thresh=extraParam[['thresh']],actfns=extraParam[['actfns']],visual=FALSE);
	ynew=mlp(inp=as.matrix(Xtest),weight=tmpNetwork$weight,dist=tmpNetwork$dist,neurons=tmpNetwork$neurons,actfns=extraParam[['actfns']]);
  } else if(classifier=="glmnet") {
    library.mc('glmnet');
    tmpFit=cv.glmnet(x=Xtrain,y=ytrain,family='binomial',alpha=0.5);
    ynew=predict(tmpFit,type='response',newx=Xtest);
    if(probabilities==FALSE) {
      ynew = round(ynew);
    }
  } else if(classifier=="GNB") {
    library.mc('klaR');
    # http://stackoverflow.com/questions/9157626/naive-bayes-in-r
    nb.res = NaiveBayes(x=Xtrain,grouping=as.factor(ytrain));
    ynew = as.integer(predict(nb.res, Xtest)$class)-1;
  } else {
    stop(paste("invalid classifier:",classifier));
  }
  return(ynew);
}
