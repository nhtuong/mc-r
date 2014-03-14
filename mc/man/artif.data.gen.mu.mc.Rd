\docType{methods}
\name{artif.data.gen.mu.mc}
\alias{artif.data.gen.mu.mc}
\title{Generate a vector of means}
\usage{
artif.data.gen.mu.mc(m, style = "linear", BayesError = NA, styleParam = 2,
  D = m)
}
\arguments{
  \item{m}{How many non-zero variables}

  \item{style}{Distribution for drawing the mus (beware,
  some of those are weird)}

  \item{BayesError}{Bayes error we want (proportion,
  between 0 and 1). If you leave NA, no adjustment will be
  made and you can end up with quite extreme mu values}

  \item{styleParam}{Optional parameter for style. Exponent
  for style exponential and linear. Useless for unique and
  uniform.}

  \item{D}{Number of variables. Make >m if you want some
  variables with µ==0}
}
\description{
Generate a vector of means appropriate to build artificial
data
}
\author{
David Dernoncourt
}

