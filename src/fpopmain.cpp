/* 
 This code and all the code in this package is publicly released on CRAN 
 under the license GPL 3.0 available at https://cran.r-project.org/web/licenses/GPL-3 
 See DESCRIPTION for further information about the authors.
*/


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
std::tuple<vector<int>, std::list<double>, vector<DeCAFS::quad>> FPOPmain (vector<double> &y, double &beta, double &lambda, double &gamma, double& phi, std::string type) {
  
  int N = y.size();
  
  // this extends the recursion to the negative phi case
  // with negative_phi = true, we reverse the order of the quadratics before crucial operations
   bool negative_phi = false;
   if (phi < 0) {
     phi = -phi;
     negative_phi = true;
   }
  
  vector<DeCAFS::quad> Q = {DeCAFS::quad(1, -INFINITY, INFINITY,
                         gamma / (1 - phi * phi),
                         -2 * y[0] * gamma / (1 - phi * phi),
                         y[0] * y[0] * gamma / (1 - phi * phi))}; // adding the first point

  vector<int> taus; // initializing the taus list
  list<vector<DeCAFS::quad>> QStorage {Q}; // initializing the cost list
  //list<double> signal;
  auto y_min = *min_element(y.begin(), y.end());
  auto y_max = *max_element(y.begin(), y.end());
  
  for (size_t t = 1; t < N; t++) {
    if (type == "isotonic") {
      //cout << "Please specify 'std' as the type argument as isotonic change is yet to be implemented." << endl;
      break;
    }
    
    // computing the increment
    auto zt = y[t] - phi * y[t - 1];
    
    // IF NEGATIVE PHI reverse the cost
    if (negative_phi) {
      Q = reverseCost(std::move(Q));
    }
    
    // getting the Qtilde
    auto Qtil = getQtil(Q, gamma, phi, zt);

    // getting the cost for no change
    vector<DeCAFS::quad> Qeq;
    if (lambda != 0 && lambda != INFINITY) {
      Qeq = infConv(Qtil, gamma * phi + lambda, y, y_min, y_max);
      Qeq = addNewPoint(move(Qeq), gamma, phi, zt);
    } else {
      Qeq = addNewPoint(Qtil, gamma, phi, zt);
    }
    

    // getting the cost for the change
    auto Qneq = infConv(Qtil, gamma * phi, y, y_min, y_max);
    Qneq = addNewPoint(move(Qneq), gamma, phi, zt);
    
    for (auto& q: Qneq) {
      get<0>(q) = t + 1;
      get<5>(q) = c(q) + beta;
    } // adding the beta penalty and updating the tau
    
    Q = getMinOfTwoQuads(Qeq, Qneq);
    
    QStorage.push_front(Q); // saving the piecewise quadratic list at time t
    
    if (negative_phi) {
      Q = reverseCost(std::move(Q));
    }
    
    //cout << "------------------------------\n" << endl;
    //cout << "Press enter to continue" << endl; cin.get();
  }
  
  if(negative_phi)
    phi = -phi;
  
  // backtracking procedures
  if (lambda != 0 && lambda != INFINITY) {
    auto signal = sigBacktrackingRWAR(move(QStorage), y, beta, lambda, gamma, phi);
    return std::make_tuple(get<0>(signal), get<1>(signal), Q);
  } else {
    auto signal = sigBacktrackingAR(move(QStorage), y, beta, gamma, phi);
    return std::make_tuple(get<0>(signal), get<1>(signal), Q);
  }
  
}
