# smartPARE

__smartPARE__ is an R package designed for sRNA cleavage confirmation based on degradome data. smartPARE utilizes deep sequencing convolutional neural networks (CNN) to segregate true and false predictions of sRNA cleavage data produced by any sRNA cleavage prediction tool. 


__Overview__ - smartPARE consists of the 3 following stages: 

1. Generation of cleavage windows 
1. Cleavage window training
1. Cleavage confirmations  

# Installation

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center", 
  warning=FALSE, 
  message=FALSE, 
  idy.opts=list(width.cutoff=120))
library(smartPARE)
```

```{r}
#Installation requires devtools
#install.packages("devtools")
devtools::install_github('kristianHoden/smartPARE')
#Load the smartPARE R package
library(smartPARE)
```

For further instructions please see the vignette.