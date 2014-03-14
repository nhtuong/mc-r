\docType{methods}
\name{artif.data.gen.model1.mc}
\alias{artif.data.gen.model1.mc}
\title{Generate simple artificial data}
\usage{
artif.data.gen.model1.mc(N = 8, D = 5, mu = 1, m = length(mu),
  sigma = 1)
}
\arguments{
  \item{N}{Number of observations (rows)}

  \item{D}{Number of variables (columns)}

  \item{mu}{µ for C1, -µ for C2. NB: can (should!) be a
  vector, typically the output of artif.data.gen.mu.mc().}

  \item{m}{Number of relevant variables (in case you feel
  like trimming mu here)}

  \item{sigma}{Standard deviation. You'll probably want to
  leave this to 1 if you care about the Bayes Error of
  artif.data.gen.mu.mc()}
}
\description{
Generate simple artificial data: model 1: in one class µ,
in the other class -µ
}
\author{
David Dernoncourt
}

