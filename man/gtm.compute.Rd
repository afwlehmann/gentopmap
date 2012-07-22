\name{gtm.compute}
\alias{gtm.compute}
\title{Compute a Generative Topographic Mapping for a high-dimensional dataset.}
\usage{
  gtm.compute(T, grid, sigma, K, epsilon = 1e-05,
    maxIterations = 100, callback = NULL, ...)
}
\arguments{
  \item{T}{a matrix whose rows contain the samples of the
  high-dimensional dataset}

  \item{grid}{a (probably regular) two- or
  three-dimensional grid in latent-space}

  \item{sigma}{the common standard deviation of the
  radially symmetric Gaussians centered at the grid points
  in latent-space}

  \item{K}{the number of samples to be drawn from the
  Gaussians in latent-space}

  \item{epsilon}{the convergence criterion}

  \item{maxIterations}{the maximum number of iterations}

  \item{callback}{a \code{function(i, ll, ...} that is
  called once during iteration with the current iteration
  number \code{i} and log-likelihood \code{ll}, as well as
  any given \code{...}}

  \item{...}{additional parameters passed on to the
  callback function}
}
\value{
  a list consisting of \item{\code{X}}{a matrix whose rows
  contain the latent-space samples} \item{\code{Rin}}{the
  responsibilities (see
  \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}})}
  \item{\code{histLogLikelihood}}{a vector containing the
  log-likelihood for each iteration}
}
\description{
  Compute a Generative Topographic Mapping for a
  high-dimensional dataset.
}
\references{
  Bishop et al., ``GTM: The Generative Topographic
  Mapping'', Neural Computation 10, No. 1, p. 215-234, 1998
}
\seealso{
  \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}}
}

