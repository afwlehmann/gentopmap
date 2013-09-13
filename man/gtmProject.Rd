\name{gtmProject}
\alias{gtmProject}
\title{Project latent-space samples using data-space responsibilities.}
\usage{
  gtmProject(X, Rin)
}
\arguments{
  \item{X}{a matrix whose rows represent the latent-space
  samples}

  \item{Rin}{the responsibilities (see
  \code{\link{computeResponsibilities}})}
}
\value{
  a matrix whose rows represent the projection of the given
  samples
}
\description{
  The given latent-space samples \eqn{X} are weighted by
  the data-space responsibilties, i.e.
  \deqn{x_{\mathit{proj}_i} = \sum_n R_{in}}{x_proj_i =
  sum_n R_in} for the \eqn{i}-th row of \eqn{X}.
}
\seealso{
  \code{\link{computeResponsibilities}}
}

