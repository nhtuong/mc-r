\docType{methods}
\name{reg.plot.mc}
\alias{reg.plot.mc}
\title{X-Y Plotting}
\usage{
reg.plot.mc(x, y, separator = NULL, type = "none", quantile = "outter",
  arrows = 0, pch, subjects = NULL, title = "CORRELATION", xlab = "X",
  ylab = "Y", col, legend.topleft, legend.topright, legend.bottomleft,
  legend.bottomright, imgfile = NULL, pointsfile = NULL)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{a numeric vector}

  \item{type}{type of regression line}

  \item{pch}{type points}

  \item{subjects}{list of labels for points}

  \item{title}{main title of plot}

  \item{xlab}{a title for the x axis}

  \item{ylab}{a title for the y axis}

  \item{col}{list of colors for points}

  \item{legend.topleft}{legend at the top-left of plot}

  \item{legend.topright}{legend at the top-right of plot}

  \item{imgfile}{image output filename}

  \item{pointsfile}{points output filename}
}
\description{
Plot a pair of variables and add regression line (linear or
lowess) to plot
}
\author{
Hoai Tuong Nguyen
}

