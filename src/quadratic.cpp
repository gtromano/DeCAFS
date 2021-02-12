/* 
 This code and all the code in this package is publicly released on CRAN 
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3 
 See DESCRIPTION for further information about the authors.
 */

/*
In this file the quadratic and the operations necessary to it are defined
*/
#include <tuple>
#include <iostream>
#include <vector>
#include <cmath>
#include "quadratic.h"
using namespace std;

//#define const double infinity =  1.0/0.0;

// this quadratic polynomial is defined in the following way:
//                 tau, from l, to u,   a X^2,  b X,    c
//typedef std::tuple<int, double, double, double, double, double> DeCAFS::quad;

// to build the quadratic simpy do DeCAFS::quad q(0, 1, 2, 3, 4, 5)

// functins for unpacking the values of tau, l, u, a, b ,c
int tau(const DeCAFS::quad &q) {
  return std::get<0>(q);
}

double l(const DeCAFS::quad &q) {
  return std::get<1>(q);
}

double u(const DeCAFS::quad &q) {
  return std::get<2>(q);
}

double a(const DeCAFS::quad &q) {
  return std::get<3>(q);
}

double b(const DeCAFS::quad &q) {
  return std::get<4>(q);
}

double c(const DeCAFS::quad &q) {
  return std::get<5>(q);
}



//////////////////////////////
//// GETTING THE MINIMUM ////
/////////////////////////////
// returns the minimum of a quadratic and where it is reached
tuple<double, double> getminimum(const DeCAFS::quad& q) {
  if (a(q) == 0) {
    return(make_tuple(c(q), l(q)));
  } else {
    auto at = (- b(q)) / (2 * a(q));
    if (at < l(q)) {
      at = l(q);
    } else if (at > u(q)) {
      at = u(q);
    }
    double minim = a(q) * (at * at) + b(q) * at + c(q);
    if (minim == -INFINITY) {
      //cout << "THIS MIN SHOULD NOT HAPPEN: " << tau(q) << "\n";
      };
    return(make_tuple(minim, at));
  }
}


/////////////////////////////////
//// GETTING INTERSECTIONS /////
////////////////////////////////

// returns the right and left intersection of a quadtratic with a line
tuple<double, double> getintersections(const DeCAFS::quad& q1, const DeCAFS::quad& q2) {
	auto ac = a(q1) - a(q2);
	auto bc = b(q1) - b(q2);
	auto cc = c(q1) - c(q2);

  double z = (bc * bc) - (4 * ac * cc);
  double left = (- bc - sqrt(z)) / (2 * ac);
  double right = (- bc + sqrt(z)) / (2 * ac);
  return(make_tuple(left, right));
}


void print_costf (std::vector<DeCAFS::quad>& costS) {
	for (auto& q: costS) {
		// cout << get<0>(q) << " "
		// 		 << get<1>(q) << " "
		// 		 << get<2>(q) << " "
		// 		 << get<3>(q) << " "
		// 		 << get<4>(q) << " "
		// 		 << get<5>(q) << endl;
	}
}

//////////////////////////////
////// main for testing /////
////////////////////////////
/*
int main(void) {
  vector<DeCAFS::quad> cost;
  cost.push_back(DeCAFS::quad(0,-INFINITY,INFINITY,1,2,4));
  cost.push_back(DeCAFS::quad(0, -INFINITY, INFINITY, 0.0, 0, 2));

  cout << l(addNewPoint(cost[0], 4.0)) << endl;
  cout << get<0>(getminimum(cost[0])) << endl;
  cout << get<1>(getminimum(cost[0])) << endl;

  cout << get<1>(getminimum(cost[1])) << endl;

  cout << get<0>(getintersections(cost[0], 6.0)) << endl;
  return 0;
}
*/
