test_that("qrcpp: decomposition works", {
  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
  h9 <- hilbert(9)
  qr1 <- qrcpp(h9, tol = 1e-10)
  expect_lt(max(abs(qr1$Q %*% qr1$R - h9[,qr1$pivot])), 1e-07)
})
