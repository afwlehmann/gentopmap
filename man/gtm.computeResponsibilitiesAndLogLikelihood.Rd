\name{gtm.computeResponsibilitiesAndLogLikelihood}
\alias{gtm.computeResponsibilitiesAndLogLikelihood}
\title{Compute responsibilities and log-likelihood.}
\usage{
  gtm.computeResponsibilitiesAndLogLikelihood(T, Y, beta)
}
\arguments{
  \item{T}{a matrix whose rows contain the samples in
  data-space}

  \item{Y}{a matrix whose rows contain the samples of X
  transformed into data-space}

  \item{beta}{the reciprocal of the variance of the
  radially symmetric Gaussians}
}
\value{
  a list consisting of \item{\code{Rin}}{the responsibility
  of the Gaussian centered at \code{Y[i,]} for having
  generated \code{T[n,]}}

  \item{\code{logLikelikehood}}{the overall log-likelihood
  of all samples and all Gaussians, i.e. \deqn{ \sum_n \log
  \left( \frac{1}{K} \sum_i p(t_n | x_i, W, \beta)
  \right)}{% sum_n( log( 1/K * sum_i( p(t_n | x_i, W, beta)
  ) ) ) } where \eqn{K} is the number of samples in
  \eqn{Y}}
}
\description{
  For each sample t in T, calculate the reponsibility of
  each Gaussian centered at one of the samples y in Y for
  having generated t.
}

