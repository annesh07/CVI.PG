#' Title Expected lower bound
#'
#' @param N number of observed data
#' @param D dimension of each observed data
#' @param T0 total clusters fixed for the variational distribution
#' @param s1 shape parameter for the prior distribution of alpha
#' @param s2 rate parameter for the prior distribution of alpha
#' @param A0 shape parameter matrix for the prior distribution of lambda
#' @param B0 rate parameter matrix for the prior distribution of lambda
#' @param X observed data, N x D matrix
#' @param W1 shape parameter for the posterior distribution of alpha
#' @param W2 rate parameter for the posterior distribution of alpha
#' @param A shape parameter matrix for the posterior distribution of lambda
#' @param B rate parameter matrix for the posterior distribution of lambda
#' @param Plog logarithm of latent allocation probability matrix
#'
#' @return value of the ELBO function
#' @export
#'
#' @examples ELBO_PG(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01,
#'              A0 = matrix(0.001, nrow = 10, ncol = 2),
#'              B0 = matrix(0.001, nrow = 10, ncol = 2),
#'              X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)),
#'                  nrow = 100, ncol = 2),
#'              W1 = 0.01, W2 = 0.01, A = matrix(0.001, nrow = 10, ncol = 2),
#'              B = matrix(0.001, nrow = 10, ncol = 2),
#'              Plog = matrix(-3, nrow = 100, ncol = 10))
ELBO_PG <- function(N, D, T0, s1, s2, A0, B0, X, W1, W2, A, B, Plog){
  P0 <- exp(Plog)
  #the alpha
  e0 <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1)-log(W2)) - s2*(W1/W2)

  #the z's
  eni <- colSums(P0)
  vni <- colSums(P0*(1 - P0))
  enj <- rowSums(apply(P0, 1, f0))
  vnj <- rowSums(apply(P0, 1, f1))
  e10 <- lgamma(1 + eni) + 0.5*trigamma(1 + eni)*vni +
    lgamma((W1/W2) + enj) + 0.5*trigamma((W1/W2) + enj)*((W1/(W2^2)) + vnj) -
    lgamma(1 + (W1/W2)+eni+enj) -
    0.5*trigamma(1 + (W1/W2)+eni+enj)*((W1/(W2^2))+vni+vnj)
  e1 <- T0*(digamma(W1) - log(W2)) + sum(e10)

  #the Lamda_ij's
  e20 <- A0*log(B0) - lgamma(A0) + (A0 - 1)*(digamma(A) - log(B)) - B0*(A/B)
  e2 <- sum(e20)

  #the x_n's
  e30 <- P0*(X %*% t(digamma(A) - log(B)))
  e31 <- - P0 %*% matrix(rowSums(A/B), nrow = T0)
  e32 <- sweep(P0, 1, - rowSums(lfactorial(X)), "*")
  e3 <- sum(e30) + sum(e31) + sum(e32)

  #the variational distributions
  e40 <- W1*log(W2) - lgamma(W1) + (W1 - 1)*(-log(W2) + digamma(W1)) - W1
  e41 <- sum(exp(Plog)*Plog)
  e42 <- A*log(B) - lgamma(A) + (A - 1)*(digamma(A) - log(B)) - A
  e4 <- e40 + e41 + sum(e42)

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
}
