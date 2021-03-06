% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleavageModelAndPrediction.R
\name{kerasCreateDataset_2d}
\alias{kerasCreateDataset_2d}
\title{Function to import data to train the model}
\usage{
kerasCreateDataset_2d(homePath = "example/", pixels = 28)
}
\arguments{
\item{homePath}{Path to your working directory}

\item{pixels}{Desired pixel size of the images}
}
\description{
Make sure to have a dir called train and your picures in dirs in this dir.
The pictures must be sorted so that a category of pics is in each subdir of the train dir.
Ex, all cat pics in a sub dir called cats and all dog pics in a subdir called dogs.
}
\examples{
kerasCreateDataset_2d(homePath = "keras/3categories/", pixels = 28)
}
\keyword{kerasCreateDataset}
