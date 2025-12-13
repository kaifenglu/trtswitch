#ifndef __LOGISTIC_REGRESSION__
#define __LOGISTIC_REGRESSION__

#include <Rcpp.h>

struct DataFrameCpp; 
struct ListCpp; 

ListCpp logisregcpp(const DataFrameCpp& data,
                    const std::string& event = "event",
                    const std::vector<std::string>& covariates = {},
                    const std::string& freq = "",
                    const std::string& weight = "",
                    const std::string& offset = "",
                    const std::string& id = "",
                    const std::string& link = "logit",
                    const std::vector<double>& init = {},
                    const bool robust = false,
                    const bool firth = false,
                    const bool flic = false,
                    const bool plci = false,
                    const double alpha = 0.05,
                    const int maxiter = 50,
                    const double eps = 1.0e-9);

Rcpp::List logisregRcpp(SEXP data,
                        const std::string& event,
                        const std::vector<std::string>& covariates,
                        const std::string& freq = "",
                        const std::string& weight = "",
                        const std::string& offset = "",
                        const std::string& id = "",
                        const std::string& link = "logit",
                        const std::vector<double>& init = {},
                        const bool robust = false,
                        const bool firth = false,
                        const bool flic = false,
                        const bool plci = false,
                        const double alpha = 0.05,
                        const int maxiter = 50,
                        const double eps = 1.0e-9);

#endif // __LOGISTIC_REGRESSION__