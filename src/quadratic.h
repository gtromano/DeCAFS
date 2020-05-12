
#ifndef ___QUADRATIC_H___
#define ___QUADRATIC_H___

#include <tuple>
#include <vector>

namespace DeCAFS {
typedef std::tuple<int, double, double, double, double, double> quad;
}

int tau(const DeCAFS::quad&);
double l(const DeCAFS::quad&);
double u(const DeCAFS::quad&);
double a(const DeCAFS::quad&);
double b(const DeCAFS::quad&);
double c(const DeCAFS::quad&);

std::tuple<double, double> getminimum(const DeCAFS::quad& );
std::tuple<double, double> getintersections(const DeCAFS::quad&, const DeCAFS::quad&);
void print_costf (std::vector<DeCAFS::quad>&);

#endif
