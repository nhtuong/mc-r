\docType{methods}
\name{samr.opt.cluster.mc}
\alias{samr.opt.cluster.mc}
\title{Optimized version of samr with clusters}
\usage{
  samr.opt.cluster.mc(df, class, label = NULL, file, index,
    cluster, type = "Two class unpaired", logged2 = T,
    nmax.perms = 200, seq = F)
}
\arguments{
  \item{df}{a microarray dataframe}

  \item{label}{vector of label}

  \item{nmax.perms}{number of maximum permutations}

  \item{logged2}{logical. Expression level has been
  transformed by logorith base 2}
}
\description{
  An optimized version of samr with cluster for gene
  differtiation analysis
}
\author{
  Hoai Tuong Nguyen
}

