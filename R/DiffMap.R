#' e-value estimation function
#'
#' Computes best e-value according to Lafon criterion
#' @param data raw data for which diffusion coordinates are calculated
#' @export

lafon <- function(data) {
  # Calculates a best e-value according to the Lafon method

  t <- as.matrix(dist(data))
  t[which(t == min(t))] = max(t)
  e <- sum(apply(t,2,mean)) / ncol(t)
  return(e)
}

#' e-value estimation function
#'
#' Computes best e-value according to Singer criterion
#' @param data raw data for which diffusion coordinates are calculated
#' @export

singer <- function(data, start, end) {
  # Calculate a best e-value according to the Singer method

  eo <- exp(seq(log(start), log(end), length.out = 100))
  e2 <- array(dim = length(eo))
  idx <- 1
  for (i in eo) {
    L <- exp(-t / (2 * i^2))
    e2[idx] <- sum(L)
    idx <- idx + 1
  }
  e3 <- e2 / max(e2)
  eb <- eo[which(abs(e3-0.5) == min(abs(e3-0.5)))]
  return(eb)
}

#' Compute diffusion coordinates
#'
#' Computes diffusion coordinates
#' @param data raw data for which diffusion coordinates are calculated
#' @param e gaussian kernel width. Optimize with either LAFON, SINGER, or custom function
#' @param a number of steps in the diffusion Markov chain
#' @export

diffmap <- function(data, e = lafon(data), a = 1) {
  # Calculate Diffusion Map Coordinates
  # Inputs
  # data  - n x G data frame of rows of cells and columns of genes
  # e - gaussian kernel width
  # a - number of diffusion steps
  # Output
  # DM - Diffusion Map Coordinates

  # Calculate Guassian distance matrix
  W <- exp(-(as.matrix(dist(data))) / e)
  D <- diag(colSums(W))
  L <- solve(D) %*% W
  for (i in 1:a) {
    W = W %*% W
    D = D %*% D
    L = solve(D) %*% L %*% solve(D)
    M = solve(D) %*% L
  }
  # Calculate Eigen vectors and values of the normalized distance matrix
  Eg = eigen(M)
  Em = t(t(data) %*% Eg[[2]])
  return(Em)
}

