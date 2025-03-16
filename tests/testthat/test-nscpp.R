library(splines)

testthat::test_that("nscpp: basis and attributes", {
  ns1 <- nscpp(women$weight, df = 5)
  ns2 <- ns(women$weight, df = 5)
  testthat::expect_equal(`attributes<-`(ns1, NULL), 
                         `attributes<-`(ns2, NULL))
  testthat::expect_equal(attr(ns1, "degree"), attr(ns2, "degree"))
  testthat::expect_equal(attr(ns1, "knots"), attr(ns2, "knots"))
  testthat::expect_equal(attr(ns1, "boundary_knots"), 
                         attr(ns2, "Boundary.knots"))
  testthat::expect_equal(attr(ns1, "intercept"), attr(ns2, "intercept"))
})


testthat::test_that("nscpp: prediction", {
  ht <- seq(57, 73, length.out = 200)
  fm1 <- lm(weight ~ nscpp(height, df = 5), data = women)
  pr1 <- predict(fm1, data.frame(height = ht))
  fm2 <- lm(weight ~ ns(height, df = 5), data = women)
  pr2 <- predict(fm2, data.frame(height = ht))
  testthat::expect_equal(pr1, pr2)
})
