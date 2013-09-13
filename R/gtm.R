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
computeResponsibilities <- function(T, Y, beta) {
  K <- nrow(Y)
  N <- nrow(T)
  D <- ncol(T)
  # Compute the responsibilities of every Gaussian centered at T[n,] for
  # sample Y[i,]
  Rin <- exp(D/2 * (log(beta) - log(2*pi)) - beta/2 *
             sapply(seq(N), function(n) rowSums(sweep(Y, 2, T[n,], `-`)^2)))
  # Compute the sums of each column as they are used for both the calculation
  # of the log-likelihood and normalization.
  auxSums <- colSums(Rin)
  # Compute the log-likelihood.
  llh <- sum(log(auxSums) - log(K))
  # Done.
  structure(
      list(Rin=scale(Rin, center=F, scale=replace(auxSums, auxSums <= 0, 1)),
           llh=llh),
      class=c("gtmResponsibilities"))
}


#' Compute a Generative Topographic Mapping for a high-dimensional dataset.
#'
#' @param T a matrix whose rows contain the samples of the high-dimensional dataset
#' @param grid a matrix whose rows correspond to the centers of the basis
#'  functions (radially symmetric Gaussians), typically layed out on a
#'  [0,1]x[0,1] grid
#' @param sigma the common standard deviation of the basis functions
#' @param K the number of latent-space samples
#' @param epsilon the convergence criterion
#' @param maxIterations the maximum number of iterations
#' @param callback a \code{function(i, ll, ...} that is called once during
#'  iteration with the current iteration number \code{i} and log-likelihood
#'  \code{ll}, as well as any given \code{...}
#' @param ... additional parameters passed on to the callback function
#' @return a list consisting of
#'  \item{\code{X}}{a matrix whose rows contain the latent-space samples}
#'  \item{\code{Rin}}{the responsibilities (see \code{\link{computeResponsibilities}})}
#'  \item{\code{histLLH}}{a vector containing the log-likelihood for each iteration}
#' @seealso \code{\link{computeResponsibilities}}
#' @references
#'      Bishop et al., ``GTM: The Generative Topographic Mapping'',
#'      Neural Computation 10, No. 1, p. 215-234, 1998
#' @importFrom mixtools rmvnorm
#' @export
computeGTM <- function(T, grid, sigma, K, epsilon=0.1, maxIterations=100, callback=NULL, ...) {
  stopifnot(is.matrix(T) && nrow(T) > 0)
  stopifnot(is.matrix(grid) && nrow(grid) > 0)
  stopifnot(sigma > 0)
  stopifnot(K > 0)    
  stopifnot(epsilon >= 0 && epsilon < 1)
  stopifnot(maxIterations > 0)

  N <- nrow(T)       # number of points in data space
  D <- ncol(T)       # dimensionality of the data space
  L <- ncol(grid)    # dimensionality of the latent space
  M <- nrow(grid)    # number of basis functions (Gaussians centered on the grid points)
  stopifnot(K %% M == 0)
  
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

  # Initialize W and beta with the least squares solution of
  # || W phi(x) - U x ||^2
  # where the columns of U are the first L principal eigenvectors of T.
  aux <- local({
    Tprime <- sweep(T, 2, colMeans(T), `-`)
    tmp <- svd(Tprime)
    U <- tmp$v[,1:L]
    ev <- tmp$d^2 / (nrow(Tprime)-1)
    pseudoInverse <- qr.solve(qr(t(Phi) %*% Phi)) %*% t(Phi)
    W <- t(pseudoInverse %*% X %*% t(U))
    beta <- 1 / ev[L+1]
    list(W=W, beta=beta)
  })
  W <- aux$W
  beta <- aux$beta
  
  # Auxiliary method for the mapping from latent- into data-space.
  computeY <- function(W, Phi) { Phi %*% t(W) }

  # Compute the projection of the samples from latent- into data space.
  Y <- computeY(W, Phi)
  
  # Shoot.
  histLLH <- NULL
  for (iter in 1:maxIterations) {
    # Compute the responsibilities and the overall log-likelihood.
    auxRin <- computeResponsibilities(T, Y, beta)
    Rin <- auxRin$Rin
    if (is.function(callback)) 
      callback(iter, auxRin$llh, ...)
    lastLLH <- tail(histLLH, 1)
    histLLH <- c(histLLH, auxRin$llh)
    if (!is.null(lastLLH) && auxRin$llh - lastLLH < log(epsilon))
        break

    # Update W.
    PhiTGPhi <- t(Phi * rowSums(Rin)) %*% Phi # same as t(Phi) %*% G %*% Phi
    invPhiTGPhi <- qr.solve(qr(PhiTGPhi))
    W <- t(invPhiTGPhi %*% t(Phi) %*% Rin %*% T)
    rm(PhiTGPhi)
    
    # Update the projection of the samples from latent- into data space (as
    # a necessary impliciation of updating beta).
    Y <- computeY(W, Phi)
    
    # Update beta accordingly.
    beta <- N * D / sum(sapply(1:N, function(n) {
      rowSums(sweep(Y, 2, T[n,], `-`)^2) * Rin[,n]
    }))
  }

  structure(list(X=X,
                 Rin=Rin,
                 histLLH=histLLH),
            class=c("gtm"))
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
  printSubset(matrix(x$histLLH, byrow=T), 3, 1)
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
#' @param Rin the responsibilities (see \code{\link{computeResponsibilities}})
#' @return a matrix whose rows represent the projection of the given samples
#' @seealso \code{\link{computeResponsibilities}}
#' @export
gtmProject <- function(X, Rin) {
  N <- ncol(Rin)
  t(sapply(1:N, function(n) colSums(Rin[,n] * X)))
}
