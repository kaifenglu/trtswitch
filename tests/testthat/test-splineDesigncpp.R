library(splines)

testthat::test_that("splineDesigncpp: default values", {
  testthat::expect_equal(splineDesigncpp(knots = 1:10, x = 4:7),
                         splineDesign(knots = 1:10, x = 4:7))
})

testthat::test_that("splineDesigncpp: derivatives", {
  testthat::expect_equal(splineDesigncpp(knots = 1:10, x = 4:7, 
                                         derivs = 0:2),
                         splineDesign(knots = 1:10, x = 4:7, 
                                      derivs = 0:2))
})
