#' Compute responsibilities and log-likelihood.
#'
#' For each sample t in T, calculate the reponsibility of each Gaussian centered
#' at one of the samples y in Y for having generated t.
#' @param T a matrix whose rows contain the samples in data-space
#' @param Y a matrix whose rows contain the samples of X transformed into
#'  data-space
#' @param beta the reciprocal of the variance of the radially symmetric
#'  Gaussians
#' @return a list consisting of
#'  \item{\code{Rin}}{the responsibility of the Gaussian centered at
#'      \code{Y[i,]} for having generated \code{T[n,]}}
#' 
#'  \item{\code{logLikelikehood}}{the overall log-likelihood of all samples and
#'      all Gaussians, i.e. \deqn{
#'          \sum_n \log \left( \frac{1}{K} \sum_i p(t_n | x_i,
#'          W, \beta) \right)}{% sum_n( log( 1/K * sum_i( p(t_n | x_i, W, beta) ) ) )
#'      } where \eqn{K} is the number of samples in \eqn{Y}}
gtm.computeResponsibilitiesAndLogLikelihood <- function(T, Y, beta) {
  K <- nrow(Y)
  N <- nrow(T)
  D <- ncol(T)
  # Compute the responsibilities of every Gaussian centered at T[n,] for
  # sample Y[i,]
  tmp <- ncol(T)/2 * (log(beta) - log(2*pi))
  Rin <- exp(tmp + -beta/2 *
             sapply(seq(N), function(n) rowSums(sweep(Y, 2, T[n,], `-`)^2)))
  # Compute the sums of each column as they are used for both the calculation
  # of the log-likelihood and normalization.
  auxSums <- colSums(Rin)
  # Compute the log-likelihood.
  logLikelihood <- sum(log(auxSums) - log(K))
  # Done.
  structure(
      list(Rin=scale(Rin, center=F, scale=replace(auxSums, auxSums <= 0, 1)),
           logLikelihood=logLikelihood),
      class=c("gtmResponsibilities"))
}


#' Compute a Generative Topographic Mapping for a high-dimensional dataset.
#'
#' @param T a matrix whose rows contain the samples of the high-dimensional dataset
#' @param grid a (probably regular) two- or three-dimensional grid in latent-space
#' @param sigma the common standard deviation of the radially symmetric
#'  Gaussians centered at the grid points in latent-space
#' @param K the number of samples to be drawn from the Gaussians in latent-space
#' @param epsilon the convergence criterion
#' @param maxIterations the maximum number of iterations
#' @param callback a \code{function(i, ll, ...} that is called once during
#'  iteration with the current iteration number \code{i} and log-likelihood
#'  \code{ll}, as well as any given \code{...}
#' @param ... additional parameters passed on to the callback function
#' @return a list consisting of
#'  \item{\code{X}}{a matrix whose rows contain the latent-space samples}
#'  \item{\code{Rin}}{the responsibilities (see \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}})}
#'  \item{\code{histLogLikelihood}}{a vector containing the log-likelihood for each iteration}
#' @seealso \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}}
#' @references
#'      Bishop et al., ``GTM: The Generative Topographic Mapping'',
#'      Neural Computation 10, No. 1, p. 215-234, 1998
#' @importFrom mixtools rmvnorm
#' @export
gtm.compute <- function(T, grid, sigma, K, epsilon=1e-5, maxIterations=100, callback=NULL, ...) {
  stopifnot(nrow(T) > 0)
  stopifnot(nrow(grid) > 0)
  stopifnot(sigma>0)
  stopifnot(K > 0)    
  stopifnot(epsilon > 0)
  stopifnot(maxIterations > 0)

  N <- nrow(T)       # number of points in data space
  D <- ncol(T)       # dimensionality of the data space
  L <- ncol(grid)    # dimensionality of the latent space
  M <- nrow(grid)    # number of basis functions (Gaussians centered on the grid points)
  stopifnot(K %% M == 0)
  
  # Define an auxiliary function that computes the inverse of a given matrix
  # by means of Cholesky decomposition.
  matrixInverse <- function(X) {
    tmp <- chol(X, pivot=TRUE)
    oo <- order(attr(tmp, "pivot"))
    chol2inv(tmp)[oo,oo]
  }
  
  # Define another auxiliary function that computes the mapping of the samples
  # from latent- into data-space.
  computeY <- function(W, Phi) { Phi %*% t(W) }
  
  # Draw samples from latent-space.
  #
  # Note that whereas in Bishop's paper the samples are stored in the columns
  # of X, here the samples are stored per row.
  numSamplesPerNode <- K / M
  X <- do.call(rbind, lapply(seq(nrow(grid)), function(i) {
    rmvnorm(n=numSamplesPerNode, mu=grid[i,], sigma=diag(rep(sigma^2, L)))
  }))
  
  # Compute Phi where Phi_ij is the probability of sample X[i,] to be drawn
  # from a Gaussian with mean grid[j,...] and variance sigma
  tmp <- -(L/2) * (log(2) + log(pi) + 2*log(sigma))
  Phi <- exp(tmp - 1/(2*sigma^2) *
             sapply(seq(M), function(i) rowSums(sweep(X, 2, grid[i,], `-`)^2)))

  # Initialize W and beta by means of the least squares solution to
  # || W phi(x) - U x ||^2 where the columns of U are the principal
  # eigenvectors of T.
  #
  # Alternatively, use D^{-1} U x instead of U x where D^{-1} is a diagonal
  # matrix whose diagonal elements are the reciprocals of the square roots
  # of the corresponding eigenvalues (thus scaling unit variance).
  ##Tprime <- scale(T)
  ##aux <- svd(Tprime)
  ##W <- aux$v[,1:L] %*% t(X) %*% Phi %*% matrixInverse(t(Phi) %*% Phi)
  ##stopifnot(nrow(W) == D && ncol(W) == M)
  ##beta <- 1 / (aux$d[L+1]^2 / (nrow(T)-1))
  eV <- eigen(cov(T))
  U <- eV$vectors[,1:L] #%*% diag(sqrt(eV$values[1:L]))
  W <- U %*% t(scale(X, scale=F)) %*% Phi %*% matrixInverse(t(Phi) %*% Phi)
  stopifnot(nrow(W) == D && ncol(W) == M)
  beta <- 1 / eV$values[L+1] # Alternatively see paper p. 8
  
  # Compute the projection of the samples from latent- into data space.
  Y <- computeY(W, Phi)
  
  # Shoot.
  histLogLikelihood <- NA
  for (iter in 1:maxIterations) {
    # Compute the responsibilities and the overall log-likelihood.
    auxRin <- gtm.computeResponsibilitiesAndLogLikelihood(T, Y, beta)
    Rin <- auxRin$Rin
    if (is.function(callback)) 
      callback(iter, auxRin$logLikelihood, ...)
    lastLLH <- tail(histLogLikelihood, 1)
    if (is.na(histLogLikelihood))
        histLogLikelihood <- auxRin$logLikelihood
    else {
      histLogLikelihood <- c(histLogLikelihood, auxRin$logLikelihood)
      # EM is guaranteed to converge, so the llh gets better at every iteration
      # step.
      if (lastLLH - auxRin$logLikelihood < log(epsilon))
          break
    }

    # Update W.
    #
    # Due to the way R handles element-wise multiplication of matrices and/or
    # vectors, Phi^T G Phi equals Phi^T * diag(G) * Phi where diag(G) is the
    # same as rowSums(Rin).
    PhiTGPhi <- t(Phi * rowSums(Rin)) %*% Phi
    W <- t(matrixInverse(PhiTGPhi) %*% t(Phi) %*% Rin %*% T)
    rm(PhiTGPhi)
    
    # Update the projection of the samples from latent- into data space (as
    # a necessary impliciation of updating beta).
    Y <- computeY(W, Phi)
    
    # Update beta accordingly.
    beta <- N * D / sum(sapply(1:N, function(n) {
      rowSums(sweep(Y, 2, T[n,], `-`)^2) * Rin[,n]
    }))
  }
  
  res <- list(X=X, Rin=Rin, histLogLikelihood=histLogLikelihood)
  class(res) <- "gtm"
  res
}


#' Formatted output of a given Generative Topographic Mapping.
#' @param x an instance of \code{gtm} (see \code{\link{gtm.compute}})
#' @param ... additional arguments to be passed on to a final print statement
#' @seealso \code{\link{gtm.compute}}
#' @method print gtm
#' @export print.gtm
print.gtm <- function(x, ...) {
  stopifnot(inherits(x, "gtm"))
  printSubset <- function(M, nx, ny) {
    fOut <- function(...) sprintf(paste(rep("%.6e", nx)), ...)
    for (i in 1:ny)
      print(paste(fOut(x$X[i,1:nx]), "...\n"))
  }
  print("Generative Topographic Mapping\n\n")
  print("Latent-space samples:\n")
  printSubset(x$X, 3, 3)
  print("\n\n")
  print("Responsibilities:\n")
  printSubset(x$Rin, 3, 3)
  print("\n\n")
  print("Log-likelihood:\n")
  printSubset(matrix(x$histLogLikelihood, byrow=T), 3, 1)
  print(...)
}


#' Project latent-space samples using data-space responsibilities.
#'
#' The given latent-space samples \eqn{X} are weighted by the
#' data-space responsibilties, i.e.
#' \deqn{x_{\mathit{proj}_i} = \sum_n R_{in}}{x_proj_i = sum_n R_in}
#' for the \eqn{i}-th row of \eqn{X}.
#' 
#' @param X a matrix whose rows represent the latent-space samples
#' @param Rin the responsibilities (see \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}})
#' @return a matrix whose rows represent the projection of the given samples
#' @seealso \code{\link{gtm.computeResponsibilitiesAndLogLikelihood}}
#' @export
gtm.project <- function(X, Rin) {
  N <- ncol(Rin)
  t(sapply(1:N, function(n) colSums(Rin[,n] * X)))
}
