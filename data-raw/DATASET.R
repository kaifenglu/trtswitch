aml = lrstat::aml
heart = lrstat::heart
immdef = lrstat::immdef
ingots = lrstat::ingots
rawdata = lrstat::rawdata
shilong = lrstat::shilong
six = lrstat::six
tobin = lrstat::tobin
sexagg = logistf::sexagg

usethis::use_data(aml, overwrite = TRUE)
usethis::use_data(heart, overwrite = TRUE)
usethis::use_data(immdef, overwrite = TRUE)
usethis::use_data(ingots, overwrite = TRUE)
usethis::use_data(rawdata, overwrite = TRUE)
usethis::use_data(shilong, overwrite = TRUE)
usethis::use_data(six, overwrite = TRUE)
usethis::use_data(tobin, overwrite = TRUE)
usethis::use_data(sexagg, overwrite = TRUE)

adsl <- readxl::read_excel("notes/adsl_sampled.xlsx")
adtdc <-readxl::read_excel("notes/adtdc_sampled.xlsx")

usethis::use_data(adsl, overwrite = TRUE)
usethis::use_data(adtdc, overwrite = TRUE)
