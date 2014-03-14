\docType{methods}
\name{stability.CW.mc}
\alias{stability.CW.mc}
\title{Calculate weighted consistency stability measure}
\usage{
stability.CW.mc(S, type = "CWrel")
}
\arguments{
  \item{S}{selection matrix, 1 row = 1 selection, 1 column
  = 1 variable. Value 0 if not selected, 1 if selected.}

  \item{type}{of stability measure. Cs: consistency. CW:
  weighted consistency. CWrel: relative weigthed
  consistency.}
}
\description{
Calculate weighted consistency stability measure and
various variants. Based on Somol 2010 "Evaluating Stability
and Comparing Output of Feature Selectors that Optimize
Feature Subset Cardinality" (doi:10.1109/tpami.2010.34)
}
\author{
David Dernoncourt
}

