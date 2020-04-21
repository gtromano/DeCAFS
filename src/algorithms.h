
#ifndef ___ALGORITHMS_H___
#define ___ALGORITHMS_H___

#include <tuple>
#include <vector>
#include <list>
#include "quadratic.h"

std::vector<quad> getMinOfTwoQuads(const std::vector<quad>&, const std::vector<quad>&);
std::vector<int> backtracking(std::vector<int>&);
std::list<double> sigBacktrackingRWAR(std::list<std::vector<quad>>, std::vector<double>&, double &, double&, double&, double&);
std::list<double> sigBacktrackingAR(std::list<std::vector<quad>>, std::vector<double>&, double &, double&, double&);
std::tuple<double, double> getGlobalMinimum(std::vector<quad>&);

std::vector<quad> recomputeIntervals(const std::vector<quad>&, const double&, const double&, const double&);
std::vector<quad> infConv(std::vector<quad>, const double&, const std::vector<double>&);
std::vector<double> generateAutoRegressive(const double&, const double&, const std::vector<double>&, const std::vector<double>&);

std::vector<quad> getQtil(std::vector<quad>, const double&,  const double&, const double&);
std::vector<quad> addNewPoint(std::vector<quad>, const double&,  const double&, const double&);

#endif
