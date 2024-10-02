\name{CIpostSelect-package}
\alias{CIpostSelect-package}
\title{CIpostSelect}
\description{
  This package calculates post-selection confidence intervals for variables.
  It uses repeated data splitting with a voting mechanism and offers two methods for post-selection:
  Lasso and BIC. For Lasso, cross-validation is used to find the best lambda that fits the model.
  For BIC, since it's not possible to test all models, a backward or forward elimination method is applied.
  The selection is done on one part of the data, followed by calibration on the other part, and this process is repeated multiple times.
}
\details{
  This package provides two main functions:
  \cr
  \cr
  - \strong{lmps} : This function provides the model selection matrices for the different data splits,
    as well as the matrix of coefficient estimates for the selected models.
    Its `summary` method gives important information about the appropriate voting type to use with the CIps function.
  \cr
  \cr
  - \strong{CIps} : This function takes an `lmps` object as a argument, along with other parameters that specify the type of vote
    and the confidence level for the confidence intervals (calculated empirically).
}
\section{Package Information}{
  \strong{Package:} CIpostSelect
  \cr
  \strong{Version:} 0.1.0
  \cr
  \strong{Date:} 2024-09-26
  \cr
  \strong{License:} MIT
}
\section{Author and Maintainer}{
  \strong{Author:} Boubacar DIALLO
  \cr
  \strong{Maintainer:} Boubacar Diallo <boubacar.diallo0@yahoo.com>
}
\examples{
\donttest{
library(mlbench)
data("BostonHousing")
# Create lmps object
model = lmps(medv ~ ., data = BostonHousing, method = "Lasso", N = 100)
# Summary of lmps
summary(model) # helps choose the appropriate vote type
# Create CIps object
cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
# Results
print(cips)
# Summary plot
plot(cips)
}
}
\docType{package}
