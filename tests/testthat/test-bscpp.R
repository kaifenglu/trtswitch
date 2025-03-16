library(splines)

testthat::test_that("bscpp: basis and attributes", {
  bs1 <- bscpp(women$weight, df = 5)
  bs2 <- bs(women$weight, df = 5)
  testthat::expect_equal(`attributes<-`(bs1, NULL), 
                         `attributes<-`(bs2, NULL))
  testthat::expect_equal(attr(bs1, "degree"), attr(bs2, "degree"))
  testthat::expect_equal(attr(bs1, "knots"), attr(bs2, "knots"))
  testthat::expect_equal(attr(bs1, "boundary_knots"), 
                         attr(bs2, "Boundary.knots"))
  testthat::expect_equal(attr(bs1, "intercept"), attr(bs2, "intercept"))
})


