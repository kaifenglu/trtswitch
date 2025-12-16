#ifndef __SURVIVAL_ANALYSIS__
#define __SURVIVAL_ANALYSIS__

#include <Rcpp.h>

struct DataFrameCpp; 
struct ListCpp; 


DataFrameCpp survQuantilecpp(const std::vector<double>& time,
                             const std::vector<int>& event,
                             const double cilevel = 0.95,
                             const std::string& transform = "loglog",
                             const std::vector<double>& probs = {0.25, 0.5, 0.75});


DataFrameCpp kmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time = "time",
                      const std::string& time2 = "",
                      const std::string& event = "event",
                      const std::string& weight = "",
                      const std::string& conftype = "log-log",
                      const double conflev = 0.95,
                      const bool keep_censor = false);


#endif // __SURVIVAL_ANALYSIS__
