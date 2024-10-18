test_that("splineDesigncpp: default values", {
  expect_equal(splineDesigncpp(knots = 1:10, x = 4:7),
               splines::splineDesign(knots = 1:10, x = 4:7))
})

test_that("splineDesigncpp: derivatives", {
  expect_equal(splineDesigncpp(knots = 1:10, x = 4:7, derivs = 0:2),
               splines::splineDesign(knots = 1:10, x = 4:7, derivs = 0:2))
})
