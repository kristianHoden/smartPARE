% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleavageModelAndPrediction.R
\name{smartPARE_examineCleavages}
\alias{smartPARE_examineCleavages}
\title{Uses the trained model to assess if the images are true or false cleavages}
\usage{
smartPARE_examineCleavages(
  examinePath,
  model,
  hitsNrs = c(1, 2),
  falseNrs = c(0),
  pixels = 28,
  picsPerloop = 100,
  avoidDir = vector(),
  reshape1 = F
)
}
\arguments{
\item{examinePath}{Path to the cleavage images (windows)}

\item{model}{Trained model}

\item{hitsNrs}{Directory number of the true cleavages. Default c(1,2)}

\item{falseNrs}{Directory number of the false cleavages. Default 0}

\item{pixels}{Pixel size of the images, same as assigned in smartPARE_train}

\item{picsPerloop}{Number of images per loop withing the function. Default 100}

\item{avoidDir}{directory that should not be assessed within the recursive examinePath}
}
\description{
Uses the trained model to assess if the images are true or false cleavages
}
\examples{
smartPARE_examineCleavages(examinePath,
model,
hitsNrs = c(1,2),
falseNrs = c(0),
pixels = 28,
picsPerloop = 100,
avoidDir = vector(),
reshape1 =F
)

}
\keyword{smartPARE_examineCleavages}
