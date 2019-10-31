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



// main l2FPOP function, takes a vector of data by reference and a penalty as a double
vector<int> FPOPmain (vector<double> &y, double &beta, double &lambda, double &gamma, double& phi, std::string type) {
  
  int N = y.size();
  vector<quad> Q = {quad(1, -INFINITY, INFINITY,
                         gamma / (1 - phi * phi),
                         -2 * y[0] * gamma / (1 - phi * phi),
                         y[0] * y[0] * gamma / (1 - phi * phi))}; // adding the first point

  vector<int> taus; // initializing the taus list

  for (size_t t = 1; t < N; t++) {
    
    //getting the minimum in Q and the relative tau
    std::vector<double> mins(Q.size());
    transform(Q.begin(), Q.end(), mins.begin(), [](quad& q){return get<0>(getminimum(q));});
    auto tau_ind = whichMin(mins); //auto tau_ind = distance(mins.begin(), min_element(mins.begin(), mins.end()));
    taus.push_back(tau(Q[tau_ind]));
    
    
    if (type == "isotonic") {
      cout << "Please specify 'std' as the type argument as isotonic change is yet to be implemented." << endl;
      break;
    }
    
    // computing the increment
    auto zt = y[t] - phi * y[t - 1];
    
    // getting the Qtilde
    auto Qtil = getQtil(move(Q), gamma, phi, zt);
    //cout << "Qtil" << endl; print_costf(Qtil);
    
    
    // getting the cost for no change
    auto Qeq = infConv(Qtil, gamma * phi + lambda, y);
    Qeq = addNewPoint(move(Qeq), gamma, phi, zt);
    //cout << "Qeq" << endl; print_costf(Qeq);
    
    
    // getting the cost for the change
    auto Qneq = infConv(move(Qtil), gamma * phi, y);
    Qneq = addNewPoint(move(Qneq), gamma, phi, zt);
    
    for (auto& q: Qneq) {
      get<0>(q) = t + 1;
      get<5>(q) = c(q) + beta;
    } // adding the beta penalty and updating the tau
    
    //cout << "Qneq" << endl; print_costf(Qneq);
    
    Q = getMinOfTwoQuads(Qeq, Qneq);
    
    //cout << "------------------------------\n" << endl;
    //cout << "Press enter to continue" << endl; cin.get();
    
  }
  
  //cout << "Q" << endl; print_costf(Q);
  //for (auto& p : taus) cout << p << endl;
  auto cp = backtracking(taus);

  return cp;
}
