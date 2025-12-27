#ifndef __SURVIVAL_ANALYSIS__
#define __SURVIVAL_ANALYSIS__

#include <Rcpp.h>

struct DataFrameCpp; 
struct ListCpp; 

std::vector<double> fsurvci(double surv, double sesurv, std::string& ct, double z);

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

DataFrameCpp kmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const double milestone = 0,
                       const double survDiffH0 = 0,
                       const double conflev = 0.95);

DataFrameCpp lrtestcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const bool weight_readj = false,
                       const double rho1 = 0,
                       const double rho2 = 0);

DataFrameCpp rmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time,
                      const std::string& event,
                      const double milestone,
                      const double conflev,
                      const bool biascorrection);

DataFrameCpp rmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& event,
                       const double milestone,
                       const double rmstDiffH0,
                       const double conflev,
                       const bool biascorrection);

#endif // __SURVIVAL_ANALYSIS__
