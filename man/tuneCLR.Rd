% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleavageModelAndPrediction.R
\name{tuneCLR}
\alias{tuneCLR}
\title{Tuning the cyclical learning rate (CLR)}
\usage{
tuneCLR(
  batch_size2 = 64,
  epochs_find_LR = 20,
  lr_max = 0.1,
  optimizer2 = keras::optimizer_rmsprop(lr = lr_max, decay = 0),
  validation_split2 = 0.2,
  rollmeanSplit = 3
)
}
\arguments{
\item{batch_size2}{Batch size}

\item{epochs_find_LR}{Epochs to run to find the optimal learnig rate, often a low number is needed}

\item{lr_max}{Maximum approved learning rate}

\item{optimizer2}{Optimizer}

\item{validation_split2}{Fraction of the data used to validate the model}

\item{rollmeanSplit}{smothens the curve by rollmean division of this number}
}
\description{
The function is run to find the minimum and maximum values of the cyclical learning rate
}
\examples{
tuneCLR(batch_size2 = 64,
epochs_find_LR = 20,
lr_max = 0.1,
optimizer2 = optimizer_rmsprop(lr=lr_max, decay=0),
validation_split2 =0.2,
rollmeanSplit = 3
)
}
\keyword{tuneCLR}
