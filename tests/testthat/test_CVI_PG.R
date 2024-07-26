test_that("Number of clusters from latent allocation:",{
  expect_equal(CVI_PG(N = 100, D = 2, T0 = 10, s1 = 0.01, s2 = 0.01,
                      A0 = matrix(0.01, nrow = 10, ncol = 2),
                      B0 = matrix(0.01, nrow = 10, ncol = 2),
                      X = matrix(c(rep(0, 50), rep(10, 50), rep(0, 50), rep(10, 50)),
                                 nrow = 100, ncol = 2),
                      W1 = 0.01, W2 = 0.01,
                      A = matrix(0.001, nrow = 10, ncol = 2),
                      B = matrix(0.001, nrow = 10, ncol = 2),
                      Plog = matrix(-3, nrow = 100, ncol = 10),
                      maxit = 1000)$posterior$Clusters, 2)})
