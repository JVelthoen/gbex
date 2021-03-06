% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dev_per_step.R
\name{dev_per_step}
\alias{dev_per_step}
\title{Deviance per step}
\usage{
dev_per_step(object, y = NULL, X = NULL)
}
\arguments{
\item{object}{A fitted gbex object}

\item{y}{vector with a test set of observations}

\item{X}{data.frame with a testset of covariates}
}
\value{
A vector with the deviance after adding a set of trees for sigma and gamma parameters
}
\description{
Compute the deviance per boosting iteration
}
\details{
Without specifying y and X the data in the gbex object are used
}
