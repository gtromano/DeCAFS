#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>      // std::iota
#include <vector>
#include <cmath>
#include <tuple>
//#include "quadratic.h"
#include "algorithms.h"
using namespace std;


template < typename Type >
int whichMin(const std::vector<Type>& v) {
  return distance(v.begin(), min_element(v.begin(), v.end()));
}



// main GFPOP function, takes a vector of data by reference and a penalty as a double
vector<int> FPOPmain (vector<double> &y, double &l0penalty, double &l2penalty, std::string type) {
  int N = y.size();
  vector<quad> Q = {addNewPoint(quad(0, -INFINITY, INFINITY, 0, 0, 0), y[0])};
  vector<int> taus;

  for (size_t t = 1; t < N; t++) {
    // applying the l2transformation
    if (l2penalty != INFINITY) { Q = applyl2Penalty(Q, l2penalty, y); }

    //getting the minimum in Q and the relative tau
    std::vector<double> mins(Q.size());
    transform(Q.begin(), Q.end(), mins.begin(), [](quad& q){return get<0>(getminimum(q));});
    auto tau_ind = whichMin(mins); //auto tau_ind = distance(mins.begin(), min_element(mins.begin(), mins.end()));
    taus.push_back(tau(Q[tau_ind]));

    vector<quad> constraint;

    if (type == "std") {
      // performing the minimization (pruning with a line at height F + beta)
      //vector<quad> constraint (1, quad(t + 1, -INFINITY, INFINITY, 0, 0, mins[tau_ind] + l0penalty));
      constraint = {quad(t + 1, -INFINITY, INFINITY, 0, 0, mins[tau_ind] + l0penalty)};

    } else if (type == "isotonic") {

      // getting the min less operator
      constraint = getCostLeq(Q, t);
      // pushing up the piecewise quadtratic of the l0penalty
      for_each(constraint.begin(), constraint.end(), [&l0penalty](quad& q) {
        get<5>(q) += l0penalty;
      });
      //cout << "costraint" << endl; print_costf(constraint);

    } else {
      cout << "Please specify a type of change, either 'std' or 'isotonic'. Exiting." << endl;
      break;
    }

    Q = getMinOfTwoQuads(Q, constraint);
    //cout << "after minimization" << endl; print_costf(Q);


    // updating the values in each piecewise quadtratic in Q
    transform(Q.begin(), Q.end(), Q.begin(), [&y, &t](quad& q){return addNewPoint(q, y[t]);}); // RENAME TO ADD NEW POINT
    //cout << "Cost function after the UPDATE" << endl; print_costf(Q);
    //cout << "------------------------------\n" << endl;
    //cout << "Press enter to continue" << endl; cin.get();
  }

  std::vector<double> mins(Q.size());
  transform(Q.begin(), Q.end(), mins.begin(), [](quad& q){return get<0>(getminimum(q));});
  auto tau_ind = whichMin(mins);
  taus.push_back(tau(Q[tau_ind]));

  auto cp = backtracking(taus);

  return cp;
}
