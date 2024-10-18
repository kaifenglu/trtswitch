test_that("bscpp: basis and attributes", {
  bs1 <- bscpp(women$weight, df = 5)
  bs2 <- splines::bs(women$weight, df = 5)
  expect_equal(`attributes<-`(bs1, NULL), `attributes<-`(bs2, NULL))
  expect_equal(attr(bs1, "degree"), attr(bs2, "degree"))
  expect_equal(attr(bs1, "knots"), attr(bs2, "knots"))
  expect_equal(attr(bs1, "boundary_knots"), attr(bs2, "Boundary.knots"))
  expect_equal(attr(bs1, "intercept"), attr(bs2, "intercept"))
})


