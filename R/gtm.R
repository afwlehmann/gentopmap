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
#'    \code{Y[i,]} for having generated \code{T[n,]}}
#'  \item{\code{llh}}{the log-likelihood of the model given the data, i.e.
#'    \deqn{
#'        \sum_n \log \left( \frac{1}{K} \sum_i p(t_n | x_i,
#'        W, \beta) \right)}{% sum_n( log( 1/K * sum_i( p(t_n | x_i, W, beta) ) ) )
#'    } where \eqn{K} is the number of samples in \eqn{Y}}
computeResponsibilities <- function(T, Y, beta) {
  K <- nrow(Y)
  N <- nrow(T)
  D <- ncol(T)
  # Compute the responsibilities of every Gaussian centered at Y[i,] for T[n,]
  #Rin <- exp(D/2 * (log(beta) - log(2*pi)) - beta/2 *
  #           sapply(seq(N), function(n) rowSums(sweep(Y, 2, T[n,], `-`)^2)))
  # Subtracting the mean vector by using matrix multiplication is way faster than
  # by using `sweep` or the like.
  Rin <- local({
    .Y <- cbind(Y, -1)
    .I <- diag(D)
    exp(D/2 * (log(beta) - log(2*pi)) - beta/2 *
        sapply(1:N, function(n) rowSums((.Y %*% rbind(.I, T[n,]))^2)))
  })
  # The sums of the columns are used for both the calculation of the
  # log-likelihood as well as normalization.
  auxSums <- colSums(Rin)
  llh <- sum(log(auxSums) - log(K))
  list(Rin=scale(Rin, center=F, scale=replace(auxSums, auxSums <= 0, 1)),
       llh=llh)
}


#' Compute a Generative Topographic Mapping for a high-dimensional dataset.
#'
#' @param T a matrix whose rows contain the samples of the high-dimensional dataset
#' @param grid a vector specifying the number of points per axis along a grid in
#'  latent-space, e.g. \code{c(100, 100)}
#' @param M the number of basis functions
#' @param sigma the common standard deviation of the basis functions
#' @param epsilon the convergence criterion
#' @param maxIterations the maximum number of iterations
#' @param verb \code{TRUE} for verbose output, else \code{FALSE}
#' @return a list, as instance of \code{gtm}, consisting of
#'  \item{\code{X}}{a matrix whose rows contain the latent-space samples}
#'  \item{\code{Rin}}{the responsibilities (see \code{\link{computeResponsibilities}})}
#'  \item{\code{mu}}{the centers of the basis functions in latent-space}
#'  \item{\code{sigma}}{the common standard deviation of the basis functions}
#'  \item{\code{beta}}{the common inverse variance of the radially symmetric
#'    Gaussians in data-space}
#'  \item{\code{W}}{a mapping from latent- into data-space}
#'  \item{\code{histLLH}}{a vector containing the log-likelihood for each iteration}
#' @seealso \code{\link{computeResponsibilities}}
#' @references
#'      Bishop et al., ``GTM: The Generative Topographic Mapping'',
#'      Neural Computation 10, No. 1, p. 215-234, 1998
#' @export
computeGTM <- function(T, grid, M, sigma, epsilon=0.1, maxIterations=50, verb=FALSE) {
  stopifnot(is.matrix(T) && nrow(T) > 0)
  stopifnot(is.vector(grid) && length(grid) > 0 && all(grid > 0))
  stopifnot(0 < M)
  stopifnot(0 < sigma)
  stopifnot(0 <= epsilon && epsilon < 1)
  stopifnot(0 < maxIterations)

  T <- cbind(T, 1)   # bias
  
  N <- nrow(T)       # number of points in data space
  D <- ncol(T)       # dimensionality of the data space
  L <- length(grid)  # dimensionality of the latent space

  # Compute the grid in latent-space.
  X <- as.matrix(do.call(expand.grid, lapply(grid, function(n) seq(-1, 1, length.out=n))))

  # Determine the centers of the basis functions.
  mu <- local({
    tmp <- rep(floor(M^(1/L)), L-1)
    numBFPerDim <- c(tmp, prod(tmp) %% M)
    as.matrix(do.call(expand.grid, lapply(numBFPerDim, function(n) seq(-1, 1, length.out=n))))
  })
  
  # Compute Phi where Phi_ij is the probability of sample X[i,] to be drawn
  # from a Gaussian with mean mu[j,...] and variance sigma
  tmp <- -(L/2) * (log(2) + log(pi) + 2*log(sigma))
  Phi <- exp(tmp - 1/(2*sigma^2) *
             sapply(seq(M), function(i) rowSums(sweep(X, 2, mu[i,], `-`)^2)))

  # Initialize W and beta with the least squares solution of
  # || W phi(x) - U x ||^2
  # where the columns of U are the first L principal eigenvectors of T.
  aux <- local({
    Tprime <- sweep(T, 2, colMeans(T), `-`)
    tmp <- svd(Tprime)
    U <- tmp$v[,1:L]
    ev <- tmp$d^2 / (nrow(Tprime)-1)
    pseudoInverse <- qr.solve(qr(t(Phi) %*% Phi)) %*% t(Phi)
    W <- pseudoInverse %*% X %*% t(U)
    beta <- 1 / ev[L+1]
    list(W=W, beta=beta)
  })
  W <- aux$W
  beta <- aux$beta
  rm(aux)
  
  # Auxiliary method for the mapping from latent- into data-space.
  computeY <- function(W, Phi) { Phi %*% W }

  # Compute the projection of the samples from latent- into data space.
  Y <- computeY(W, Phi)
  
  # Shoot.
  histLLH <- NULL
  for (iter in 1:maxIterations) {
    # Compute the responsibilities and the overall log-likelihood.
    auxRin <- computeResponsibilities(T, Y, beta)
    Rin <- auxRin$Rin

    if (verb)
        cat(sprintf("%.18e\n", auxRin$llh))

    # See if we have already converged.
    lastLLH <- tail(histLLH, 1)
    histLLH <- c(histLLH, auxRin$llh)
    if (!is.null(lastLLH) && auxRin$llh - lastLLH < log(epsilon))
        break

    # Update W.
    PhiTGPhi <- t(Phi * rowSums(Rin)) %*% Phi # same as t(Phi) %*% G %*% Phi
    invPhiTGPhi <- qr.solve(qr(PhiTGPhi))
    W <- invPhiTGPhi %*% t(Phi) %*% Rin %*% T
    # Clean up after ourselves.
    rm(PhiTGPhi, invPhiTGPhi)
    
    # Update the mapping's parameter set before updating beta.
    Y <- computeY(W, Phi)
    beta <- N * D / sum(sapply(1:N, function(n) {
      rowSums(sweep(Y, 2, T[n,], `-`)^2) * Rin[,n]
    }))
  }
  if (length(histLLH) < 2 || diff(tail(histLLH, 2)) >= log(epsilon))
      warning("GTM did not converge.")

  structure(list(X=X,
                 Rin=Rin,
                 mu=mu,
                 sigma=sigma,
                 beta=beta,
                 W=W,
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


#' Posterior-mean projection of the samples used for training the given GTM.
#' @param model an instance of \code{gtm}
#' @return a LxN matrix whose rows represent the projection of the N samples
#'   into the L-dimensional latent-space
#' @seealso \code{\link{computeResponsibilities}}
#' @export
gtmPosteriorMean <- function(model) {
  stopifnot(inherits(model, "gtm"))
  N <- ncol(model$Rin)
  t(sapply(1:N, function(n) colSums(model$Rin[,n] * model$X)))
}


#' Posterior-mode projection of the samples used for training the given GTM.
#' @param model an instance of \code{gtm}
#' @return a LxN matrix whose rows represent the projection of the N samples
#'   into the L-dimensional latent-space
#' @seealso \code{\link{computeResponsibilities}}
#' @export
gtmPosteriorMode <- function(model) {
  stopifnot(inherits(model, "gtm"))
  N <- ncol(model$Rin)
  t(sapply(1:N, function(n) model$X[which.max(model$Rin[,n]),]))
}
