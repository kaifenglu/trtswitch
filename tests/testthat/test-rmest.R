library(dplyr, warn.conflicts = FALSE)

test_that("rmest: valid milestone time works", {
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time", 
               event="status", milestone=161)
  
  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  rmst = 52.64545
  stderr = 19.8286
  
  expect_equal(round(df1$rmst, 5), rmst)
  expect_equal(round(df1$stderr, 4), stderr)
})


test_that("rmest: bias correction for variance estimate", {
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time", 
               event="status", milestone=161, biascorrection=TRUE)
  
  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(bc tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  stderr = 21.4173
  
  expect_equal(round(df1$stderr, 4), stderr)
})


test_that("rmest: milestone should be <= largest observed time", {
  expect_error(rmest(aml %>% filter(x == "Nonmaintained"), time="time", 
               event="status", milestone=161))
})
