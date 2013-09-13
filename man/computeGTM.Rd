\name{computeGTM}
\alias{computeGTM}
\title{Compute a Generative Topographic Mapping for a high-dimensional dataset.}
\usage{
  computeGTM(T, grid, sigma, K, epsilon = 0.1,
    maxIterations = 100, callback = NULL, ...)
}
\arguments{
  \item{T}{a matrix whose rows contain the samples of the
  high-dimensional dataset}

  \item{grid}{a matrix whose rows correspond to the centers
  of the basis functions (radially symmetric Gaussians),
  typically layed out on a [0,1]x[0,1] grid}

  \item{sigma}{the common standard deviation of the basis
  functions}

  \item{K}{the number of latent-space samples}

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
  \code{\link{computeResponsibilities}})}
  \item{\code{histLLH}}{a vector containing the
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
  \code{\link{computeResponsibilities}}
}

