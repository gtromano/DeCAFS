/* 
 This code and all the code in this package is publicly released on CRAN 
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3 
 See DESCRIPTION for further information about the authors.
 */

#ifndef ___ALGORITHMS_H___
#define ___ALGORITHMS_H___

#include <tuple>
#include <vector>
#include <list>
#include "quadratic.h"

std::vector<DeCAFS::quad> getMinOfTwoQuads(const std::vector<DeCAFS::quad>&, const std::vector<DeCAFS::quad>&);
std::vector<int> backtracking(std::vector<int>&);
std::tuple<std::vector<int>, std::list<double>> sigBacktrackingRWAR(std::list<std::vector<DeCAFS::quad>>, std::vector<double>&, double &, double&, double&, double&);
std::tuple<std::vector<int>, std::list<double>> sigBacktrackingAR(std::list<std::vector<DeCAFS::quad>>, std::vector<double>&, double &, double&, double&);
std::tuple<double, double> getGlobalMinimum(std::vector<DeCAFS::quad>&);

std::vector<DeCAFS::quad> recomputeIntervals(const std::vector<DeCAFS::quad>&, const double&, const double&, const double&);
std::vector<DeCAFS::quad> infConv(std::vector<DeCAFS::quad>, const double&, const std::vector<double>&);
std::vector<double> generateAutoRegressive(const double&, const double&, const std::vector<double>&, const std::vector<double>&);
std::vector<DeCAFS::quad>  reverseCost (std::vector<DeCAFS::quad>);
  
std::vector<DeCAFS::quad> getQtil(std::vector<DeCAFS::quad>, const double&,  const double&, const double&);
std::vector<DeCAFS::quad> addNewPoint(std::vector<DeCAFS::quad>, const double&,  const double&, const double&);

#endif
