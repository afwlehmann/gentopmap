\name{computeGTM}
\alias{computeGTM}
\title{Compute a Generative Topographic Mapping for a high-dimensional dataset.}
\usage{
  computeGTM(T, grid, M, sigma, epsilon = 0.1,
    maxIterations = 50, verb = FALSE)
}
\arguments{
  \item{T}{a matrix whose rows contain the samples of the
  high-dimensional dataset}

  \item{grid}{a vector specifying the number of points per
  axis along a grid in latent-space, e.g. \code{c(100,
  100)}}

  \item{M}{the number of basis functions}

  \item{sigma}{the common standard deviation of the basis
  functions}

  \item{epsilon}{the convergence criterion}

  \item{maxIterations}{the maximum number of iterations}

  \item{verb}{\code{TRUE} for verbose output, else
  \code{FALSE}}
}
\value{
  a list, as instance of \code{gtm}, consisting of
  \item{\code{X}}{a matrix whose rows contain the
  latent-space samples} \item{\code{Rin}}{the
  responsibilities (see
  \code{\link{computeResponsibilities}})}
  \item{\code{mu}}{the centers of the basis functions in
  latent-space} \item{\code{sigma}}{the common standard
  deviation of the basis functions} \item{\code{beta}}{the
  common inverse variance of the radially symmetric
  Gaussians in data-space} \item{\code{W}}{a mapping from
  latent- into data-space} \item{\code{histLLH}}{a vector
  containing the log-likelihood for each iteration}
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

