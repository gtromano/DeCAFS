
#ifndef ___QUADRATIC_H___
#define ___QUADRATIC_H___

#include <tuple>
#include <vector>

typedef std::tuple<int, double, double, double, double, double> quad;

int tau(const quad&);
double l(const quad&);
double u(const quad&);
double a(const quad&);
double b(const quad&);
double c(const quad&);

std::tuple<double, double> getminimum(const quad& );
std::tuple<double, double> getintersections(const quad&, const quad&);
void print_costf (std::vector<quad>&);

#endif
