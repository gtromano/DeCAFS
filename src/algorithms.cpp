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
#include "quadratic.h"
using namespace std;

template < typename Type >
int whichMin(const std::vector<Type>& v) {
  return distance(v.begin(), min_element(v.begin(), v.end()));
}

template < typename Type >
int whichMin(const std::list<Type>& v) {
  return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}


//////////////////////////////////////////////
///// COMBINE TWO PIECEWISE QUADRATICS ///////
//////////////////////////////////////////////

// this function performs the minimization between a cost function and a costraint
std::list<DeCAFS::quad> getMinOfTwoQuads(const std::list<DeCAFS::quad>& costS, const std::list<DeCAFS::quad>& costR) {
  std::list<DeCAFS::quad> outcost;
  std::list<DeCAFS::quad>::const_iterator i = costS.begin(); // iterator for the previous cost function
  std::list<DeCAFS::quad>::const_iterator j = costR.begin(); // iterator for the contraint function (this will always be equal to 0 in FPOP)
  double low = -INFINITY; // index for the linesearch
  double upp;
  
  while (low != INFINITY) {
    upp = std::min(u(*i), u(*j));
    
    // checking a piecewise quadratic with a line
    // we will always have to check a piecewise quadtratic with a line in case of FPOP
    //auto minS = get<0>(getminimum(costS[i]));
    //auto minCr = get<0>(getminimum(costR[j]));
    auto minS = get<0>(getminimum(DeCAFS::quad(tau(*i), low, upp, a(*i), b(*i), c(*i))));
    auto minCr = get<0>(getminimum(DeCAFS::quad(tau(*j), low, upp, a(*j), b(*j), c(*j))));
    
    auto inters = getintersections(*i, *j);
    auto left_in_range = (get<0>(inters) > low) && (get<0>(inters) <= upp);
    auto right_in_range = (get<1>(inters) > low) && (get<1>(inters) <= upp);
    auto in_range = left_in_range + right_in_range;
    
    if (in_range == 0) {
      // we do not have any interaction in range, so a check is needed to see if
      // the cost or the constraint are greater
      if (minCr <= minS) {
        outcost.push_back(DeCAFS::quad(tau(*j), low, upp, a(*j), b(*j), c(*j)));
      } else {
        outcost.push_back(DeCAFS::quad(tau(*i), low, upp, a(*i), b(*i), c(*i)));
      }
      
    } else if (in_range == 1) {
      // in this case we have only one interaction and we have to understand if it's
      // the right one or the left one
      if (right_in_range) {
        // in this case it's the right
        outcost.push_back(DeCAFS::quad(tau(*i), low,            get<1>(inters), a(*i), b(*i), c(*i)));
        outcost.push_back(DeCAFS::quad(tau(*j), get<1>(inters), upp,            a(*j), b(*j), c(*j)));
      } else {
        // in this case is the left
        outcost.push_back(DeCAFS::quad(tau(*j), low,            get<0>(inters), a(*j), b(*j), c(*j)));
        outcost.push_back(DeCAFS::quad(tau(*i), get<0>(inters), upp,            a(*i), b(*i), c(*i)));
      }
      
    } else if (in_range == 2) {
      // in this case we have two interaction points
      outcost.push_back(DeCAFS::quad(tau(*j), low,            get<0>(inters), a(*j), b(*j), c(*j)));
      outcost.push_back(DeCAFS::quad(tau(*i), get<0>(inters), get<1>(inters), a(*i), b(*i), c(*i)));
      outcost.push_back(DeCAFS::quad(tau(*j), get<1>(inters), upp,            a(*j), b(*j), c(*j)));
    }
    
    // updating the lower bound
    low = upp;
    // updating the iterators
    if (upp == u(*i)) {
      i++;
    }
    if (upp == u(*j)) {
      j++;
    }
  }
  
  return outcost;
}

      

///////////////////////////////////////////
///// RECONSTRUCT THE CHANGEPOINTS  ///////
//////////////////////////// //////////////

std::vector<int> backtracking(std::vector<int>& taus) {
  std::vector<int> cp;
  //int s = taus.size() + 1;
  int s = taus.size();
  //cp.push_back(s);
  
  // for (auto& k:taus) { std::cout << k << ' '; }
  // std::cout << "\n";
  
  // s - 2 since 0 indexing
  while (s != 1) {
    cp.push_back(s);
    s = taus[s - 2];
    /*
    cout << "Current tau " << s << endl;
    cin.get();
    */
  }
  reverse(cp.begin(), cp.end());
  return cp;
}



////////////////////////////////////////////////////////////
///// RECOMPUTE THE INTERVALS OF DISJOINT QUADRATICS ///////
////////////////////////////////////////////////////////////

// This is a O(n) version. It has been tested to work with the simple FPOP
// It might work with the GFPOP, however a proof is needed since at the current state
// the O(n^2) algorithm (still to code) is the best one
// std::vector<DeCAFS::quad> recomputeIntervals(const std::vector<DeCAFS::quad>& cost, const double& lower, const double& upper, const double& sigma) {
// 
//   // initializing an empty vector
//   std::vector<DeCAFS::quad> outcost (3 * cost.size(), DeCAFS::quad(0, 0, 0, 0, 0, 0));
// 
//   std::vector<int> index(cost.size()); // vector with n ints needed as a reference index for pruning
//   iota (std::begin(index), std::end(index), 0); // Fill with 0, 1, ..., n
// 
//   int i = 0; // index for the current cost function
//   int k = 0; // index for the outcost
//   double prevInter = -INFINITY; // previous interaction
// 
//   while (true) {
// 
//     // if we are at the end of our list
//     if ((i + 1) == index.size()) {
//       outcost[k] = cost[index[i]];
//       get<1>(outcost[k]) = prevInter; // from to the previous interaction
//       get<2>(outcost[k]) = INFINITY; // to beyond (actually infinity)
//       break; // terminate if we are at the end
// 
//     } else {
//       // getting (the bigger) intersection between the current and the next one
//       auto inter = get<1>(getintersections(cost[index[i]], cost[index[i+1]]));
// 
//       ////// JUST SOME SAFETY CHECKS ON THE INTERSECTION /////
//       if (inter < lower || inter > upper) {inter = NAN;} // if out of range make it NAN
//       // if there is no valid intersection then quit SHOULD NOT HAPPEN
//       if (isnan(inter)) {
//         //cout << "Detected no intersection in Range, stopping cicle" << endl;
//         outcost[k] = cost[index[i]];
//         get<1>(outcost[k]) = prevInter; // from to the previous interaction
//         get<2>(outcost[k]) = INFINITY; // to beyond (actually infinity)
//         break; // terminate if we are at the end
//       }
//       ////////////////////////////////////////////////
// 
//       if (inter <= prevInter) {
// 
//         index.erase(index.begin() + i); // erasing the previous quadratic from the index
//         // reverting the changes made in the previous step
//         i--;
//         k--;
//         prevInter = l(outcost[k]); // restoring the previous intersection
// 
//       } else {
//         outcost[k] = cost[index[i]];
//         get<1>(outcost[k]) = prevInter; // from to the previous interaction
//         get<2>(outcost[k]) = inter; // to beyond (actually infinity)
//         k++;
//         i++;
//         prevInter = inter;
//       }
//     }
//   }
//   outcost.resize(k + 1);
//   return outcost;
// }

// std::list<DeCAFS::quad> recomputeIntervals(std::list<DeCAFS::quad>& cost, const double& lower, const double& upper, const double& sigma) {
//   
//   // initializing an empty list
//   std::list<DeCAFS::quad> outcost;
//   
//   std::list<DeCAFS::quad>::iterator i = cost.begin(); // iterator for the current cost function
//   double prevInter = -INFINITY; // previous interaction
//   
//   while (true) {
//     
//     // if we are at the end of our list
//     if (std::next(i) == cost.end()) {
//       outcost.push_back(*i);
//       get<1>(outcost.back()) = prevInter; // from to the previous interaction
//       get<2>(outcost.back()) = INFINITY; // to beyond (actually infinity)
//       break; // terminate if we are at the end
//       
//     } else {
//       // getting (the bigger) intersection between the current and the next one
//       auto inter = get<1>(getintersections(*i, *std::next(i)));
//       
//       ////// JUST SOME SAFETY CHECKS ON THE INTERSECTION /////
//       if (inter < lower || inter > upper) {inter = NAN;} // if out of range make it NAN
//       // if there is no valid intersection then quit SHOULD NOT HAPPEN
//       if (isnan(inter)) {
//         //cout << "Detected no intersection in Range, stopping cicle" << endl;
//         outcost.push_back(*i);
//         get<1>(outcost.back()) = prevInter; // from to the previous interaction
//         get<2>(outcost.back()) = INFINITY; // to beyond (actually infinity)
//         break; // terminate if we are at the end
//       }
//       ////////////////////////////////////////////////
//       
//       if (inter <= prevInter) {
//         
//         cost.erase(i); // erasing the previous quadratic from the list
//         // reverting the changes made in the previous step
//         i = cost.begin();
//         outcost.pop_back();
//         prevInter = l(outcost.back()); // restoring the previous intersection
//         
//       } else {
//         outcost.push_back(*i);
//         get<1>(outcost.back()) = prevInter; // from to the previous interaction
//         get<2>(outcost.back()) = inter; // to beyond (actually infinity)
//         i++;
//         prevInter = inter;
//       }
//     }
//   }
//   return outcost;
// }

// std::list<DeCAFS::quad> recomputeIntervals(std::list<DeCAFS::quad>& cost, const double& lower, const double& upper, const double& sigma) {
//   
//   // initializing an empty list
//   std::list<DeCAFS::quad> outcost;
//   
//   std::list<DeCAFS::quad>::iterator i = cost.begin(); // iterator for the current cost function
//   double prevInter = -INFINITY; // previous interaction
//   
//   while (true) {
//     
//     // if we are at the end of our list
//     if (std::next(i) == cost.end()) {
//       outcost.push_back(*i);
//       get<1>(outcost.back()) = prevInter; // from to the previous interaction
//       get<2>(outcost.back()) = INFINITY; // to beyond (actually infinity)
//       break; // terminate if we are at the end
//       
//     } else {
//       // getting (the bigger) intersection between the current and the next one
//       auto inter = get<1>(getintersections(*i, *std::next(i)));
//       
//       ////// JUST SOME SAFETY CHECKS ON THE INTERSECTION /////
//       if (inter < lower || inter > upper) {inter = NAN;} // if out of range make it NAN
//       // if there is no valid intersection then quit SHOULD NOT HAPPEN
//       if (isnan(inter)) {
//         //cout << "Detected no intersection in Range, stopping cicle" << endl;
//         outcost.push_back(*i);
//         get<1>(outcost.back()) = prevInter; // from to the previous interaction
//         get<2>(outcost.back()) = INFINITY; // to beyond (actually infinity)
//         break; // terminate if we are at the end
//       }
//       ////////////////////////////////////////////////
//       
//       if (inter <= prevInter) {
//         
//         cost.erase(i); // erasing the previous quadratic from the list
//         // and then we restart! 
//         i = cost.begin();
//         
//         outcost.clear();
// 
//         prevInter = -INFINITY;
//         
//       } else {
//         outcost.push_back(*i);
//         get<1>(outcost.back()) = prevInter; // from to the previous interaction
//         get<2>(outcost.back()) = inter; // to beyond (actually infinity)
//         i++;
//         prevInter = inter;
//       }
//     }
//   }
//   return outcost;
// }

void recomputeIntervals(std::list<DeCAFS::quad>& cost, const double& lower, const double& upper, const double& sigma) {
  
  // initializing an empty list
  // std::list<DeCAFS::quad> outcost;
  
  std::list<DeCAFS::quad>::iterator i = cost.begin(); // iterator for the current cost function
  double prevInter = -INFINITY; // previous interaction
  
  while (true) {
    
    // if we are at the end of our list
    if (std::next(i) == cost.end()) {
      //outcost.push_back(*i);
      get<1>(*i) = prevInter; // from to the previous interaction
      get<2>(*i) = INFINITY; // to beyond (actually infinity)
      break; // terminate if we are at the end
      
    } else {
      // getting (the bigger) intersection between the current and the next one
      auto inter = get<1>(getintersections(*i, *std::next(i)));
      
      ////// JUST SOME SAFETY CHECKS ON THE INTERSECTION /////
      if (inter < lower || inter > upper) {inter = NAN;} // if out of range make it NAN
      // if there is no valid intersection then quit SHOULD NOT HAPPEN
      if (isnan(inter)) {
        //cout << "Detected no intersection in Range, stopping cicle" << endl;
        //outcost.push_back(*i);
        get<1>(*i) = prevInter; // from to the previous interaction
        get<2>(*i) = INFINITY; // to beyond (actually infinity)
        break; // terminate if we are at the end
      }
      ////////////////////////////////////////////////
      
      if (inter <= prevInter) {
        
        cost.erase(i); // erasing the previous quadratic from the list
        // and then we restart! 
        i = cost.begin();
        
        //outcost.clear();
        
        prevInter = -INFINITY;
        
      } else {
        //outcost.push_back(*i);
        get<1>(*i) = prevInter; // from to the previous interaction
        get<2>(*i) = inter; // to beyond (actually infinity)
        i++;
        prevInter = inter;
      }
    }
  }
//  return outcost;
}


/*
// this function gets what in Rigaill & al. is known as the min less then equal constraint
std::vector<DeCAFS::quad> getCostLeq(std::vector<DeCAFS::quad>& cost, const int& t) {
  
  // initializing an empty vector
  std::vector<DeCAFS::quad> outcost (2 * cost.size(), DeCAFS::quad(0, 0, 0, 0, 0, 0));
  
  
  int i = 0; // index for the current quadratic
  int j = 1; // index for the quadratic we look the interaction with
  int k = 0; // index for the outcost
  
  while (true) {
    
    auto MINandAT = getminimum(cost[i]);
    
    // if this is true it means we're at the end of the list
    if ((i + j) >= cost.size()) {
      
      // adding the current quadtratic to the output matrix
      outcost[k] = cost[i];
      get<2>(outcost[k]) = get<1>(MINandAT); // to the point where the min is reached
      k++;
      // adding the constraint line
      outcost[k] = DeCAFS::quad(t + 1, get<1>(MINandAT), INFINITY, 0, 0, get<0>(MINandAT));
      break; // terminate
      
    } else {
      
      // this is the candidate line for quadratic i, starting from its minimum
      DeCAFS::quad candidateLine = DeCAFS::quad(t + 1, get<1>(MINandAT), INFINITY, 0, 0, get<0>(MINandAT));
      
      // checking if there is an intersection with the line and with the next quadtratic
      auto inter = get<1>(getintersections(candidateLine, cost[i + j]));
      auto in_range = (inter > l(cost[i + j])) && (inter < u(cost[i + j]));
      
      // if so we have found an interaction
      if (in_range) {
        
        // add the current DeCAFS::quad to the output
        outcost[k] = cost[i];
        get<2>(outcost[k]) = get<1>(MINandAT);
        k++;
        // add the line to the output till inter
        outcost[k] = candidateLine;
        get<2>(outcost[k]) = inter;
        k++;
        // set the lower of the next quadratic to the inter
        get<1>(cost[i + j]) = inter;
        
        
        i = i + j; // we jump at the quadtratic we found the interaction with
        j = 1; // we reset the j
        
      } else {
        // here we jump to the second next quadtratic to see if there is an interaction
        j = j + 1;
      } // end if - else
    } // end if - else
  } // end while
  
  outcost.resize(k + 1);  // we resize and return the cost
  return outcost;
}
*/


/////////////////////////////////////////////////////////////
///// INFIMUM CONVOLUTION APPLIED TO A PIECEWISE QUAD ///////
/////////////////////////////////////////////////////////////

// std::vector<DeCAFS::quad> infConv(std::vector<DeCAFS::quad> cost, const double& omega, const std::vector<double>& y) {
//   
//   //cout << "Applying the transormation!" << endl;
//   //cout << "cost before transformation" << endl; print_costf(cost);
//   
//   
//   // applying the tranformation and obtaining a set of disjoint quadratics
//   // please pay attention at the order of operation to avoid overwriting a, b or c before they're used!
//   for_each(cost.begin(), cost.end(), [&omega](auto& q){
//     auto minAt = get<1>(getminimum(q));
//     get<1>(q) = minAt + (l(q) - minAt) * (a(q) / omega + 1); // not necessary in case 1
//     get<2>(q) = minAt + (u(q) - minAt) * (a(q) / omega + 1); // not necessary in case 1
//     get<5>(q) =  c(q) - ((b(q) * b(q)) / (4 * (a(q) + omega)));
//     get<4>(q) = b(q) * (omega / (a(q) + omega));
//     get<3>(q) = a(q) * (omega / (a(q) + omega));
//   });
//   
//   //cout << "cost before erase" << endl; print_costf(cost);
//   
//   
//   // REMOVE OBSERVATIONS OUTISE THE RANGE OF Y
//   auto y_min = *min_element(y.begin(), y.end());
//   auto y_max = *max_element(y.begin(), y.end());
//   
//   cost.erase(remove_if(cost.begin(), cost.end(), [&y_min, &y_max] (auto& q){
//     return(((u(q) < y_min) || (l(q) > y_max)));
//     }), cost.end());
//   
//   get<1>(cost[0]) = -INFINITY;
//   get<2>(cost[cost.size() - 1]) = INFINITY;
//   ////////////////////////////////////////////
//   
//   
//   //cout << "cost after erase" << endl; print_costf(cost);
//   
//   // if we're comparing just lines we pick the smallest one
//   auto allLines = true;
//   for (const auto& q : cost){
//     if (a(q) != 0) allLines = false;
//   }
//   if (allLines) {
//     //cout << "running" << endl;
//     std::vector<double> mins(cost.size());
//     transform(cost.begin(), cost.end(), mins.begin(), [](DeCAFS::quad& q){return c(q);});
//     cost = {cost[whichMin(mins)]};
//   }
//   /////////////////////////////////////////////
//   
//   
//   
//   
//   // recomputing the intervals
//   cost = recomputeIntervals(cost, -INFINITY, INFINITY, 0.0);
//   
//   return cost;
// }
void infConv(std::list<DeCAFS::quad>& cost, const double& omega, const std::vector<double>& y, const double& y_min, const double& y_max) {
  // applying the transformation and obtaining a set of disjoint quadratics
  // please pay attention at the order of operation to avoid overwriting a, b or c before they're used!
  for_each(cost.begin(), cost.end(), [&omega](auto& q){
    auto minAt = get<1>(getminimum(q));
    get<1>(q) = minAt + (l(q) - minAt) * (a(q) / omega + 1); // not necessary in case 1
    get<2>(q) = minAt + (u(q) - minAt) * (a(q) / omega + 1); // not necessary in case 1
    get<5>(q) =  c(q) - ((b(q) * b(q)) / (4 * (a(q) + omega)));
    get<4>(q) = b(q) * (omega / (a(q) + omega));
    get<3>(q) = a(q) * (omega / (a(q) + omega));
  });
  
  // REMOVE OBSERVATIONS OUTISE THE RANGE OF Y
  cost.remove_if([&y_min, &y_max] (auto& q){
    return(((u(q) < y_min) || (l(q) > y_max)));
  });
  
  get<1>(cost.front()) = -INFINITY;
  get<2>(cost.back()) = INFINITY;
  
  // if we're comparing just lines we pick the smallest one
  auto allLines = true;
  for (const auto& q : cost){
    if (a(q) != 0) allLines = false;
  }
  if (allLines) {
    std::list<double> mins;
    transform(cost.begin(), cost.end(), std::back_inserter(mins), [](DeCAFS::quad& q){return c(q);});
    //cost = {cost[whichMin(mins)]};
    auto minElement = whichMin(mins);
    
    auto it = cost.begin();
    advance(it, minElement);
    cost = {*it};
  }
  
  // recomputing the intervals
  //cost = recomputeIntervals(cost, -INFINITY, INFINITY, 0.0);
  recomputeIntervals(cost, -INFINITY, INFINITY, 0.0);  
  //return cost;
}



//////////////////////////////
///// UPDATE FUNCTIONS ///////
//////////////////////////////

/*
DeCAFS::quad addNewPointNOPENALTY(const DeCAFS::quad& q, const double& yi) {
  DeCAFS::quad newq(tau(q), l(q), u(q), a(q) + 1, b(q) - 2 * yi, c(q) + yi * yi);
  return(newq);
}
*/

// 
// std::vector<DeCAFS::quad> addNewPoint(std::vector<DeCAFS::quad> cost, const double& gamma,  const double& phi, const double& zt) {
//   
//   for_each(cost.begin(), cost.end(), [&gamma, &zt, &phi](DeCAFS::quad& q){
//     get<3>(q) = a(q) - gamma * phi + gamma;
//     get<4>(q) = b(q) - 2 * zt * gamma;
//     get<5>(q) = c(q) - gamma * (zt * zt) / (phi - 1);
//   });
//   
//   return cost;
// }

void addNewPoint(std::list<DeCAFS::quad>& cost, const double& gamma,  const double& phi, const double& zt) {
  std::for_each(cost.begin(), cost.end(), [&gamma, &zt, &phi](DeCAFS::quad& q){
    get<3>(q) = a(q) - gamma * phi + gamma;
    get<4>(q) = b(q) - 2 * zt * gamma;
    get<5>(q) = c(q) - gamma * (zt * zt) / (phi - 1);
  });
  
  //return cost;
}

// ////////////////////////////////
// ///// GET Q TIL FUNCTION ///////
// ////////////////////////////////
// 
// std::vector<DeCAFS::quad> getQtil(std::vector<DeCAFS::quad> cost, const double& gamma,  const double& phi, const double& zt) {
//   
//   for_each(cost.begin(), cost.end(), [&gamma, &zt, &phi](DeCAFS::quad& q){
//     get<3>(q) = a(q) - gamma * phi * (1 - phi);
//     get<4>(q) = b(q) + 2 * gamma * phi * zt;
//     get<5>(q) = c(q) - gamma * phi * (zt * zt) / (1 - phi);
//   });
//   
//   return cost;
// }

std::list<DeCAFS::quad> getQtil(std::list<DeCAFS::quad> cost, const double& gamma,  const double& phi, const double& zt) {
  std::for_each(cost.begin(), cost.end(), [&gamma, &zt, &phi](DeCAFS::quad& q){
    get<3>(q) = a(q) - gamma * phi * (1 - phi);
    get<4>(q) = b(q) + 2 * gamma * phi * zt;
    get<5>(q) = c(q) - gamma * phi * (zt * zt) / (1 - phi);
  });
  
  return cost;
}
////////////////////////////////////////
///// reverse the mu axis /////////////
//////////////////////////////////////

// std::vector<DeCAFS::quad> reverseCost (std::vector<DeCAFS::quad> cost) {
//   
//   for_each(cost.begin(), cost.end(), [](DeCAFS::quad& q){
//     get<4>(q) = -b(q); // inverting the b coefficient
//     
//     auto temp = u(q);
//     get<2>(q) = -l(q);
//     get<1>(q) = -temp;
//   });
//   
//   std::reverse(std::begin(cost), std::end(cost));
//   
//   return(cost); // that's all folks! 
// }
std::list<DeCAFS::quad> reverseCost (std::list<DeCAFS::quad> cost) {
  
  for_each(cost.begin(), cost.end(), [](DeCAFS::quad& q){
    get<4>(q) = -b(q); // inverting the b coefficient
    
    auto temp = u(q);
    get<2>(q) = -l(q);
    get<1>(q) = -temp;
  });
  
  std::reverse(std::begin(cost), std::end(cost));
  
  return(cost); // that's all folks! 
}

///////////////////////////////////
////// get global minimum ////////
/////////////////////////////////
std::tuple<double, double> getGlobalMinimum(std::list<DeCAFS::quad>& Q) {
  
  std::vector<double> mins(Q.size());
  transform(Q.begin(), Q.end(), mins.begin(), [](DeCAFS::quad& q){return get<0>(getminimum(q));});
  
  std::vector<double> ats(Q.size());
  transform(Q.begin(), Q.end(), ats.begin(), [](DeCAFS::quad& q){return get<1>(getminimum(q));});
  
  
  auto minElement = whichMin(mins);
  auto min = mins[minElement];
  auto at = ats[minElement];
  return std::make_tuple(min, at);
}


/////////////////////////////////
////// evaluate cost func ///////
/////////////////////////////////

double evalCost(const std::list<DeCAFS::quad>& Q, const double& at) {
  for (auto& q : Q) {
    if (l(q) < at && u(q) > at) return a(q) * (at * at) + b(q) * at + c(q);
  }
  return nan("");
}

//////////////////////////////
///// SIGNAL BACKTRACKING ////
//////////////////////////////

std::tuple<std::vector<int>, std::list<double>> sigBacktrackingRWAR(std::list<std::list<DeCAFS::quad>> QStorage, vector<double>& y, double &beta, double& lambda, double& gamma, double& phi) {
  int N = y.size();
  std::list<double> muHatStorage;
  double muHat;
  
  // initialization
  muHat = get<1>(getGlobalMinimum(QStorage.front()));
  std::list<int> tauHatStorage;
  muHatStorage.push_front(muHat);
  QStorage.pop_front();
  
  auto t = N - 2;
  
  for (auto& Qt : QStorage) {
    
    // making the first b piecewise function (in case of a change)
    std::list<DeCAFS::quad> B1(Qt.size());
    transform(Qt.begin(), Qt.end(), B1.begin(), [&t, &muHat, &y, &beta, &gamma, &phi](DeCAFS::quad& q){
      auto zt = y[t + 1] - phi *  y[t];
      DeCAFS::quad newq(tau(q),
                l(q),
                u(q),
                a(q) + gamma * (phi * phi),
                b(q) + 2 * gamma * (zt - muHat) * phi,
                c(q) + beta + gamma * (zt - muHat) * (zt - muHat));
      return newq;
    });
    auto B1Min = getGlobalMinimum(B1);
    
    
    // making the second b piecewise function (in case of no change)
    std::list<DeCAFS::quad> B2(Qt.size());
    transform(Qt.begin(), Qt.end(), B2.begin(), [&t, &muHat, &y, &lambda, &gamma, &phi](DeCAFS::quad& q){
      auto zt = y[t + 1] - phi *  y[t];
      DeCAFS::quad newq(tau(q),
                l(q),
                u(q),
                a(q) + gamma * (phi * phi) + lambda,
                b(q) - 2 * lambda * muHat + 2 * gamma * (zt - muHat) * phi,
                c(q) + lambda * muHat * muHat + gamma * (zt - muHat) * (zt - muHat));
      return newq;
    });
    
    
    auto B2Min = getGlobalMinimum(B2);
    //cout<< get<1>(B1Min) << "   " << get<1>(B2Min) << endl;
    // if the minimum of the first cost function is smaller than the second take its argmin
    
    if (get<0>(B1Min) >= get<0>(B2Min)) {
      muHat = get<1>(B2Min);
    } else {
      muHat = get<1>(B1Min);
      tauHatStorage.push_front(t + 1);
    }        

    muHatStorage.push_front(muHat);
    t -= 1;
  } // end for
  
  
  std::vector<int> cp_n{ std::make_move_iterator(tauHatStorage.begin()), 
                         std::make_move_iterator(tauHatStorage.end())};
  
  return std::make_tuple(cp_n, muHatStorage);
}



std::tuple<std::vector<int>, std::list<double>> sigBacktrackingAR(std::list<std::list<DeCAFS::quad>> QStorage, vector<double>& y, double &beta, double& gamma, double& phi) {
  int N = y.size();
  std::list<int> tauHatStorage;
  std::list<double> muHatStorage;
  double muHat;
  
  //cout << QStorage.size() << endl;
  //print_costf(QStorage.front());
  
  // initialization
  muHat = get<1>(getGlobalMinimum(QStorage.front()));
  muHatStorage.push_front(muHat);
  QStorage.pop_front();
  
  auto t = N - 2;
  
  for (auto& Qt : QStorage) {
    
    // making the first b piecewise function (in case of a change)
    std::list<DeCAFS::quad> B1(Qt.size());
    transform(Qt.begin(), Qt.end(), B1.begin(), [&t, &muHat, &y, &beta, &gamma, &phi](DeCAFS::quad& q){
      auto zt = y[t + 1] - phi *  y[t];
      DeCAFS::quad newq(tau(q), l(q), u(q),
                a(q) + gamma * (phi * phi),
                b(q) + 2 * gamma * (zt - muHat) * phi,
                c(q) + beta + gamma * (zt - muHat) * (zt - muHat));
      return newq;
    });
    auto B1Min = getGlobalMinimum(B1);
    
    
    // making the second b piecewise function (in case of no change)
    std::list<DeCAFS::quad> B2(Qt.size());
    transform(Qt.begin(), Qt.end(), B2.begin(), [&t, &muHat, &y, &gamma, &phi](DeCAFS::quad& q){
      auto zt = y[t + 1] - phi *  y[t];
      DeCAFS::quad newq(tau(q), l(q), u(q),
                a(q) + gamma * (phi * phi),
                b(q) + 2 * gamma * (zt - muHat) * phi,
                c(q) + gamma * (zt - muHat) * (zt - muHat));
      return newq;
    });
    auto Bstar = evalCost(B2, muHat);

    //cout<< get<0>(B1Min) << "   " << Bstar << endl;
    if (get<0>(B1Min) < Bstar) {
      muHat = get<1>(B1Min);
      tauHatStorage.push_front(t + 1);
    }
    
    muHatStorage.push_front(muHat);
    t -= 1;
  } // end for
  
  std::vector<int> cp_n{ std::make_move_iterator(tauHatStorage.begin()), 
                         std::make_move_iterator(tauHatStorage.end())};
  
  return std::make_tuple(cp_n, muHatStorage);
}


/////////////////////////
///// GENERATE AR ///////
/////////////////////////

// this function generates a ar process [to be intended to call from R function dataRW]
std::vector<double> generateAutoRegressive(const double& phi, const double& y0, const std::vector<double>& mu, const std::vector<double>& ynoise) {
  std::vector<double> y(ynoise.size());
  //y[0] = y0 * phi + mu[0] + ynoise[0];
  y[0] += mu[0];
  for (size_t t = 1; t < (ynoise.size()); t++) {
    y[t] = phi * (y[t-1] - mu[t - 1]) + mu[t] + ynoise[t];
  }
  return y;
}
