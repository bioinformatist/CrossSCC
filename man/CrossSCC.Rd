% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CrossSCC.R
\docType{package}
\name{CrossSCC}
\alias{CrossSCC}
\title{CrossSCC: An Single-Cell Clustering analysis framework Crossing batches of data}
\usage{
CrossSCC(
  m,
  markers,
  ncores = 16,
  verbose = R.utils::Verbose(threshold = -1, timestamp = TRUE),
  show.progress.bar = TRUE,
  log.file = NULL
)
}
\arguments{
\item{m}{a matrix of single cell expression values.}

\item{markers}{a dataframe contains marker information. Try examples for details.}

\item{ncores}{number of CPU cores used.}

\item{verbose}{verbose level. By default, CrossSCC will output all logs as well as progress bars.}

\item{show.progress.bar}{Set to FALSE if you don't want to see progress bar.}

\item{log.file}{a test/debug use option. Messages will be directed to this new text file instead of stderr().
Note that verbose threshold will be set to -1 with this parameter.}
}
\value{
a data.tree object.
}
\description{
CrossSCC: An Single-Cell Clustering analysis framework Crossing batches of data

CrossSCC
}
\examples{
library(data.table)
markers <- fread(system.file("extdata", "markers.csv", package = "CrossSCC"))
data(immu)
handsome.zuo <- CrossSCC(immu, markers, ncores = 16)
}
