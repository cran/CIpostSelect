% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uneFonction.R
\name{CIps}
\alias{CIps}
\title{Creates an object of class CIps based on the provided parameters.}
\usage{
CIps(x, vote, alpha, s.vote_coef = 0.5)
}
\arguments{
\item{x}{An object of class lmps, which contains the selection and coefficient estimation matrices.}

\item{vote}{The type of vote to perform: "model" for selection based on the most frequent model,
or "coef" for variable selection (e.g., if a variable is selected more than 50 percent of the time).}

\item{alpha}{Specifies the confidence level for the confidence intervals.}

\item{s.vote_coef}{A parameter between 0 and 1 that, when using "coef" voting,
indicates the frequency threshold for selecting a variable.}
}
\value{
An object of class CIps.
}
\description{
Creates an object of class CIps based on the provided parameters.
}
\details{
After obtaining the lmps object, which provides the selection matrices (models and coefficients),
this function allows us to compute confidence intervals that are calculated empirically based on the chosen voting method and the desired level of certainty.
The confidence intervals are obtained through empirical calculation on each vector of estimates for the corresponding coefficients.

CIps also provides an intercept (test version) estimated as follows: in the case of a vote on models, it takes the average of the intercept vector for the rows where the most frequently selected model in the N splits is chosen. For the vote on coefficients, the idea is to select the coefficient that has been chosen the least number of times among those retained and then average the intercept only for the rows where this coefficient is selected.
}
\examples{
library(mlbench)
data("BostonHousing")
# lmps object
model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50)
# CIps object
cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)



\donttest{
# lmps object
model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
# CIps object
cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
}
}
