#' Title Collapsed Variational Inference (CVI) for Poisson-Gamma model
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
#' @param maxit maximum number of iterations the updates are run
#'
#' @return the posterior value of alpha, the number of clusters from the
#' latent allocation probability, the proportions of each cluster and the
#' latent allocation probability matrix
#' @export
#'
#' @examples CVI_PG(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01,
#'              A0 = matrix(0.001, nrow = 10, ncol = 2),
#'              B0 = matrix(0.001, nrow = 10, ncol = 2),
#'              X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)),
#'                  nrow = 100, ncol = 2),
#'              W1 = 0.01, W2 = 0.01, A = matrix(0.001, nrow = 10, ncol = 2),
#'              B = matrix(0.001, nrow = 10, ncol = 2),
#'              Plog = matrix(-3, nrow = 100, ncol = 10), maxit = 1000)
CVI_PG <- function(N, D, T0, s1, s2,A0, B0, X, W1, W2, A, B, Plog, maxit){
  f <- list()
  f[[1]] <- ELBO_PG(N, D, T0, s1, s2,A0, B0, X, W1, W2, A, B, Plog)

  maxit <- 1000

  for (m in 1:maxit){
    P0 <- exp(Plog)
    #different updates for i = 1, i = {2, ..., T0-1} and i = T0
    P230 <- X %*% t(digamma(A) - log(B))
    P231 <- rowSums(A/B)
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]

      p_eni <- colSums(P1)
      p_vni <- colSums(P1*(1-P1))
      p_enj <- rowSums(apply(P1, 1, f0))
      p_vnj <- rowSums(apply(P1, 1, f1))

      P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
                                                            (W1/W2)) + (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj +
                                                                                                         (W1/W2))^2)

      P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
      P22 <- c(0, cumsum(P21)[1:(T0-1)])

      # P231 <- rep(NA, T0)
      # for (i in 1:T0){
      #   P231[i] <- -0.5*L21[i,, drop=FALSE] %*% t(L21[i,, drop=FALSE])
      # }
      P2 <- P20 + P22 + P230[n, ] - P231
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }

    #labelling the probability matrix so that non-zero cluster allocations
    #are present in the beginning of the matrix
    P00 <- exp(Plog)
    Csum <- colSums(P00)
    index <- which(Csum > 0.0001)
    l0 <- length(index)
    for (l in 1:l0){
      Plog[, c(l, index[l])] <- Plog[, c(index[l], l)]
    }
    #final updated and labelled probability matrix
    Pf0 <- exp(Plog)

    A <- A0 + t(Pf0) %*% X
    B <- sweep(B0, 1, colSums(Pf0), "+")

    #update of the shape parameter of alpha
    W1 <- s1 + l0 - 1
    #update of the rate parameter of alpha
    a0 <- l0/log(N)
    a_eni <- colSums(Pf0[,1:l0])
    a_vni <- colSums(Pf0[,1:l0]*(1-Pf0[,1:l0]))
    a_enj <- rowSums(apply(Pf0[,1:l0], 1, f0))
    a_vnj <- rowSums(apply(Pf0[,1:l0], 1, f1))
    W20 <- log(a0 + a_eni[1:(l0 - 1)] + a_enj[1:(l0 - 1)]) -
      0.5*(a_vni[1:(l0 - 1)] + a_vnj[1:(l0 - 1)])/((a0 + a_eni[1:(l0 - 1)]
                                                    + a_enj[1:(l0 - 1)])^2) - log(a0 + a_enj[1:(l0 - 1)]) +
      0.5*a_vnj[1:(l0 - 1)]/((a0 + a_enj[1:(l0 - 1)])^2)
    W21 <- log(a0 + a_eni[l0]) - 0.5*a_vni[l0]/((a0 + a_eni[l0])^2) -
      log(a0)
    W2 <- s2  + sum(W20) + W21

    f[[m+1]] <- ELBO_PG(N, D, T0, s1, s2,A0, B0, X, W1, W2, A, B, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }
  alpha0 <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))

  posterior <- list("alpha"=alpha0, "Clusters"=clustnum,
                    "Proportions"=clust, "Clustering" = Plog)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
