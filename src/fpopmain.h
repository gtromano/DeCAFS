/* 
 This code and all the code in this package is publicly released on CRAN 
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3 
 See DESCRIPTION for further information about the authors.
 */

#ifndef ___FPOPGAETA_H___
#define ___FPOPGAETA_H___
#include "quadratic.h"


std::tuple<std::vector<int>, std::list<double>, std::vector<DeCAFS::quad>> FPOPmain (std::vector<double>&, double&, double&, double&, double&, std::string);

#endif
