**Table of Contents:**





# "mc" R package: Omics Data Analysis #

## Bioconductor-like synopsis ##

Version: 1.3 (Released on October 01, 2014)

Description: Statistical and datamining tools for omics data analysis.

Author(s): **Hoai Tuong Nguyen, David Dernoncourt, Meriem Abdennour**

Maintainer(s): **Hoai Tuong Nguyen** <[hoai-tuong.nguyen@inserm.fr](mailto:hoai-tuong.nguyen@inserm.fr)>, **David Dernoncourt** <[me@daviddernoncourt.com](mailto:me@daviddernoncourt.com)>

## Idea ##
You are researchers working with OMICS (genomics, metagenomics, metabolomics, etc) data. You are not keen on coding, but want to get quickly the relevant results from your data. Here, we make and share the ready-to-run codes to generate the ready-to-publish results.

## Join us ##
To join us or inquire for further informations, please inbox our maintainer(s).

## Citation ##
Please kindly cite us as **"mc" R package** within publications that make use of methods inspired by this work.


# Getting started #

## Installation ##
To install this package, install R (>3.0) start R and enter:

```
source(pipe("wget -O - https://mc-r.googlecode.com/svn/trunk/mc/R/mcLite.R"))
```

**Note:** We use the function "pipe" to deal with "https" protocol.
For Windows user:
```
source("https://mc-r.googlecode.com/svn/trunk/mc/R/mcLite.R")
```

## Examples ##

### Create correlation heatmap ###
The package "[corrplot](http://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)" offers various options to plot a correlation heatmap. However, there is no options to add p-value and correlation coefficent to the heatmap. The function "corrplot.mc" allows you to do so.

http://data2papers.files.wordpress.com/2014/07/corr_heatmap.png?w=840&h=630

**Legend:** Correlation heatmap. Cross symbols are added at non-significant correlations.

```
library.mc("corrplot")
corrplot.mc(cor(mtcars), type = "lower", method = "pie", diag=F, sig.level = 0.05, sig = "p-value",p.cex=0.8) 
```

**Note:** You may want to look futher information from "corrplot" project on CRAN for more options and examples: http://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html



### Create a multilayer graph ###

You would like to create a multilayer graph, for instance, to compare the changes of graph through the time.

https://data2papers.files.wordpress.com/2014/07/multilayer_network.png?w=740&h=429

**Legend:** Multilayer graph with three layers. Dot-lines represent the reference between layers.


The function "getMultilayerNet.mc" is designed for this purpose. Here is an example:

```
library.mc("igraph")
library.mc("RCytoscape") 

#Create 3 separate graphs
##Create graph #1
g1 = new ("graphNEL", edgemode = "directed")
g1 = initNodeAttribute (g1, "nodeType", "char", "undefined")
g1 = initNodeAttribute (g1, "label", "char", "undefined")
g1 = initEdgeAttribute (g1, "edgeType", "char", "undefined")
g1 = initEdgeAttribute (g1, "weight", "numeric", 1.0)

## add nodes graph #1
g1 = addNode ("A", g1)
g1 = addNode ("B", g1)
g1 = addNode ("C", g1)
g1 = addNode ("D", g1)

## copy nodes for graph #2 & graph #3
g2<-g3<-g1

# add edges for graph #1
g1 = addEdge ("A", "B", g1)
g1 = addEdge ("A", "C", g1)
g1 = addEdge ("B", "D", g1)
g1 = addEdge ("C", "D", g1)

## add edges for graph #2
g2 = addEdge ("B", "C", g2)
g2 = addEdge ("A", "C", g2)

## add edges for graph #3
g3 = addEdge ("A", "C", g3)


##Prepare input for getMultilayerNet.mc
gl<-list(g1,g2,g3)
#gl<-list(g1,g2,g3,g2,g1)
###offset from the first graph
offset.x=500
offset.y=-500
###create multilayer graph
getMultilayerNet.mc(gl,offset.x,offset.y)
```

**Note:** Please enable Cytoscape 2.8 RPC plugin before. You perhaps want to get to know how to do this. Please check out this link: http://rcytoscape.systemsbiology.net/versions/current/index.html


### Load & install a package ###
To use a package, just execute the function "library.mc". It loads automatically the library if existed, otherwise it installs automatically from the global servers such as CRAN (repos="cran"), Bioconductor (repos="bioc"):
```
#syntax 1: library.mc("name-of-new-package") 
#syntax 2: library.mc("name-of-new-package",repos="cran") # repos=c("cran","bioc")
library.mc("timecourse")
```

**Tips:** Most of packages on CRAN are also available on Bioconductor. Hence, you have more change to install successfully with the option <repos="bioc">. In the function "library.mc", this option is set as default.