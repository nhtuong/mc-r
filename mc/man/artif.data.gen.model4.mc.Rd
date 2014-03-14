\docType{methods}
\name{artif.data.gen.model4.mc}
\alias{artif.data.gen.model4.mc}
\title{Generate simple artificial data with correlation}
\usage{
artif.data.gen.model4.mc(N = 8, D = 5, mu = 1, m = length(mu),
  sigma = 1, corBlockSize = 2, corStrength = 0.5)
}
\arguments{
  \item{N}{Number of observations (rows)}

  \item{D}{Number of variables (columns)}

  \item{mu}{µ for C1, -µ for C2. NB: can (should!) be a
  vector, typically the output of artif.data.gen.mu.mc().}

  \item{m}{Number of relevant variables (in case you feel
  like trimming mu here)}

  \item{sigma}{Standard deviation.}

  \item{corBlockSize}{Size of correlated variables blocks}

  \item{corStrength}{Correlation between correlated
  variables (within blocks)}
}
\description{
Generate simple artificial data with correlation: model 4:
like model 1 but with blocks of correlated variables
}
\author{
David Dernoncourt
}

