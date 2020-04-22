#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>      // std::iota
#include <vector>
#include <list>
#include <cmath>
#include <tuple>
//#include "quadratic.h"
#include "algorithms.h"
using namespace std;

template < typename Type >
int whichMin(const std::vector<Type>& v) {
  return distance(v.begin(), min_element(v.begin(), v.end()));
}



// main l2FPOP function, takes a vector of data by reference and a penalty as a double
std::tuple<vector<int>, std::list<double>, vector<quad>> FPOPmain (vector<double> &y, double &beta, double &lambda, double &gamma, double& phi, std::string type) {
  
  int N = y.size();
  
  vector<quad> Q = {quad(1, -INFINITY, INFINITY,
                         gamma / (1 - phi * phi),
                         -2 * y[0] * gamma / (1 - phi * phi),
                         y[0] * y[0] * gamma / (1 - phi * phi))}; // adding the first point

  vector<int> taus; // initializing the taus list
  list<vector<quad>> QStorage {Q}; // initializing the cost list
  //list<double> signal;
  
  for (size_t t = 1; t < N; t++) {
    
    //getting the minimum in Q and the relative tau
    std::vector<double> mins(Q.size());
    transform(Q.begin(), Q.end(), mins.begin(), [](quad& q){return get<0>(getminimum(q));});
    auto tau_ind = whichMin(mins);
    taus.push_back(tau(Q[tau_ind]));
    //signal.push_back(get<1>(getGlobalMinimum(Q)));
    
    
    if (type == "isotonic") {
      //cout << "Please specify 'std' as the type argument as isotonic change is yet to be implemented." << endl;
      break;
    }
    
    // computing the increment
    auto zt = y[t] - phi * y[t - 1];
    
    // getting the Qtilde
    auto Qtil = getQtil(Q, gamma, phi, zt);
    //QStorage.push_front(Q); // saving the piecewise quadratic list at time t
    
    // getting the cost for no change
    vector<quad> Qeq;
    if (lambda != 0 && lambda != INFINITY) {
      Qeq = infConv(Qtil, gamma * phi + lambda, y);
      Qeq = addNewPoint(move(Qeq), gamma, phi, zt);
    } else {
      Qeq = addNewPoint(Qtil, gamma, phi, zt);
    }
    

    // getting the cost for the change
    auto Qneq = infConv(Qtil, gamma * phi, y);
    Qneq = addNewPoint(move(Qneq), gamma, phi, zt);
    
    for (auto& q: Qneq) {
      get<0>(q) = t + 1;
      get<5>(q) = c(q) + beta;
    } // adding the beta penalty and updating the tau
    
    Q = getMinOfTwoQuads(Qeq, Qneq);
    QStorage.push_front(Q); // saving the piecewise quadratic list at time t
    
    //cout << "------------------------------\n" << endl;
    //cout << "Press enter to continue" << endl; cin.get();
    
  }
  
  //cout << "Q" << endl; print_costf(Q);
  //for (auto& p : taus) cout << p << endl;
  std::vector<double> mins(Q.size());
  transform(Q.begin(), Q.end(), mins.begin(), [](quad& q){return get<0>(getminimum(q));});
  auto tau_ind = whichMin(mins);
  taus.push_back(tau(Q[tau_ind]));
  
  auto cp = backtracking(taus);
  
  if (lambda != 0 && lambda != INFINITY) {
    auto signal = sigBacktrackingRWAR(move(QStorage), y, beta, lambda, gamma, phi);
    return std::make_tuple(cp, signal, Q);
  } else {
    auto signal = sigBacktrackingAR(move(QStorage), y, beta, gamma, phi);
    return std::make_tuple(cp, signal, Q);
  }
  
}
