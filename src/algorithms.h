
#ifndef ___ALGORITHMS_H___
#define ___ALGORITHMS_H___

#include <tuple>
#include <vector>
#include "quadratic.h"

std::vector<quad> getMinOfTwoQuads(std::vector<quad>&, std::vector<quad>&);
std::vector<int> backtracking(std::vector<int>&);
std::vector<quad> recomputeIntervals(std::vector<quad>&, const double&, const double&, const double&);
std::vector<quad> getCostLeq(std::vector<quad>&, const int&);
std::vector<quad> applyl2Penalty(std::vector<quad>&, const double&, const std::vector<double>&);

#endif
