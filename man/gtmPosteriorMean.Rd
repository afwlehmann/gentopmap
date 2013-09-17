\name{gtmPosteriorMean}
\alias{gtmPosteriorMean}
\title{Posterior-mean projection of the samples used for training the given GTM.}
\usage{
  gtmPosteriorMean(model)
}
\arguments{
  \item{model}{an instance of \code{gtm}}
}
\value{
  a LxN matrix whose rows represent the projection of the N
  samples into the L-dimensional latent-space
}
\description{
  Posterior-mean projection of the samples used for
  training the given GTM.
}
\seealso{
  \code{\link{computeResponsibilities}}
}

