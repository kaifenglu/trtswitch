library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("zph_phregr: right-censored data", {
  fit1 <- phregr(data = liver, time = "Time", event = "Status", 
                 covariates = c("log(Bilirubin)", "log(Protime)", 
                                "log(Albumin)", "Age", "Edema"))
  
  fit2 <-  coxph(Surv(Time, Status) ~ log(Bilirubin) + log(Protime) + 
                   log(Albumin) + Age + Edema, data = liver)
  
  zph1 <- zph_phregr(fit1, transform = "log")
  zph2 <- cox.zph(fit2, transform = "log")
  
  rownames(zph2$y) <- NULL
  zph2$call <- NULL
  
  testthat::expect_equal(zph1, zph2)
})

