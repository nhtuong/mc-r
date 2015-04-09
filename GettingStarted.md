### Installation ###
To install this package, start R and enter:

```
source(pipe("wget -O - https://dl.dropboxusercontent.com/s/rlj5vhgso52pb6j/mcLite.R"))
```

**Note:** We use the function "pipe" to deal with "https" protocol.

### Examples ###
#### Load & install a package ####
To use a package, just execute the function "library.mc". It loads automatically the library if existed, otherwise it installs automatically from the global servers such as CRAN (repos="cran"), Bioconductor (repos="bioc"):
```
#Syntax 1: library.mc("name-of-new-package") 
#Syntax 2: library.mc("name-of-new-package",repos="cran") # repos=c("cran","bioc")
library.mc("timecourse")
```

**Tips:** Most of packages on CRAN are also available on Bioconductor. Hence, you have more change to install successfully with the option <repos="bioc">. In the function "library.mc", this option is set as default.

#### Create a multilayer graph ####

You would like to create a multilayer graph, for instance, to compare the changes of graph through the time.

https://data2papers.files.wordpress.com/2014/07/multilayer_network.png?w=740&h=429

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
g1 = addEdge ("A", "B", g)
g1 = addEdge ("A", "C", g)
g1 = addEdge ("B", "D", g)
g1 = addEdge ("C", "D", g)

## add edges for graph #2
g2 = addEdge ("B", "C", g2)
g2 = addEdge ("A", "C", g2)

## add edges for graph #3
g3 = addEdge ("A", "C", g3)


##Prepare input for getMultilayerNet.mc
gl<-list(g1,g2,g3)
###offset from the first graph
offset.x=500
offset.y=-500
###create multilayer graph
getMultilayerNet.mc(gl,500,-500)
```

**Note:** Please enable Cytoscape 2.8 RPC plugin before. You perhaps want to get to know how to do this. Please check out this link: http://rcytoscape.systemsbiology.net/versions/current/index.html