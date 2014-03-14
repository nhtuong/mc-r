\docType{methods}
\name{artif.data.getDeltaForBE.mc}
\alias{artif.data.getDeltaForBE.mc}
\title{Calculate coefficient to obtain target Bayes error}
\usage{
artif.data.getDeltaForBE.mc(wanted_error, mus)
}
\arguments{
  \item{wanted_error}{Bayes error we want (proportion,
  between 0 and 1)}

  \item{mus}{A vector of means for separation (one class
  would have delta*mus, the other -delta*mus)}
}
\description{
Calculate a delta coefficient so as to obtain target Bayes
error when using delta*mu to build artificial data
(assuming variance=1)
}
\author{
David Dernoncourt
}

