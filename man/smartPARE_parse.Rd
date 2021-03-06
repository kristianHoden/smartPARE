% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleavageModelAndPrediction.R
\name{smartPARE_parse}
\alias{smartPARE_parse}
\title{Parses the smartPARE output file into a dataset}
\usage{
smartPARE_parse(
  smartPAREresultFile = paste0(homePath1, "example_results/true.txt")
)
}
\arguments{
\item{smartPAREresultFile}{smartPARE output filepath}
}
\description{
Parses the smartPARE output file into a dataset based on the path of each true image
}
\examples{
smartPARE_parse(smartPAREresultFile = paste0(homePath1,"example_results/true.txt"))

}
\keyword{smartPARE_parse}
