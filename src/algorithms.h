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

std::list<DeCAFS::quad> getMinOfTwoQuads(const std::list<DeCAFS::quad>&, const std::list<DeCAFS::quad>&);
std::vector<int> backtracking(std::vector<int>&);
std::tuple<std::vector<int>, std::list<double>> sigBacktrackingRWAR(std::list<std::list<DeCAFS::quad>>, std::vector<double>&, double &, double&, double&, double&);
std::tuple<std::vector<int>, std::list<double>> sigBacktrackingAR(std::list<std::list<DeCAFS::quad>>, std::vector<double>&, double &, double&, double&);
std::tuple<double, double> getGlobalMinimum(std::list<DeCAFS::quad>&);

std::list<DeCAFS::quad> recomputeIntervals(const std::list<DeCAFS::quad>&, const double&, const double&, const double&);
std::list<DeCAFS::quad> infConv(std::list<DeCAFS::quad>, const double&, const std::vector<double>&);
std::vector<double> generateAutoRegressive(const double&, const double&, const std::vector<double>&, const std::vector<double>&);
std::list<DeCAFS::quad>  reverseCost (std::list<DeCAFS::quad>);
  
std::list<DeCAFS::quad> getQtil(std::list<DeCAFS::quad>, const double&,  const double&, const double&);
std::list<DeCAFS::quad> addNewPoint(std::list<DeCAFS::quad>, const double&,  const double&, const double&);

#endif
