\docType{methods}
\name{samt.opt.mc}
\alias{samt.opt.mc}
\title{Optimized version of samr}
\usage{
  samt.opt.mc(df, label, nmax.perms = 100, logged2 = T,
    type = "Two class unpaired", seq = F)
}
\arguments{
  \item{df}{a microarray dataframe}

  \item{label}{vector of label}

  \item{nmax.perms}{number of maximum permutations}

  \item{logged2}{logical. Expression level has been
  transformed by logorith base 2}
}
\description{
  An optimized version of samr for gene differtiation
  analysis
}
\author{
  Hoai Tuong Nguyen
}

