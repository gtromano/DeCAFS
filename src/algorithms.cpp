#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>      // std::iota
#include <vector>
#include <cmath>
#include <tuple>
#include "quadratic.h"
using namespace std;


// this function performs the minimization between a cost function and a costraint
std::vector<quad> getMinOfTwoQuads(std::vector<quad>& costS, std::vector<quad>& costR) {
  std::vector<quad> outcost (3 * costS.size(), quad(0, 0, 0, 0, 0, 0));
  int i = 0; // index for the previous cost function
  int j = 0; // index for the contraint function (this will always be equal to 0 in FPOP)
  int k = 0; // index for the outcost
  double low = -INFINITY; // index for the linesearch
  double upp;
  
  while (low != INFINITY) {
    upp = min(u(costS[i]), u(costR[j]));
    
    if (a(costR[j]) != 0) {
      // checking a piecewise quadtratic with another quadtratic
      // in this case we are comparing two quadratics, and costS is always the best by a factor of l0penalty
      outcost[k] = costS[i];
      get<1>(outcost[k]) = low;
      get<2>(outcost[k]) = upp;
      k++;
      
    } else {
      
      // checking a piecewise quadratic with a line
      // we will always have to check a piecewise quadtratic with a line in case of FPOP
      auto minS = get<0>(getminimum(costS[i]));
      auto minCr = get<0>(getminimum(costR[j]));
      
      auto inters = getintersections(costS[i], quad(0, 0, 0, 0, 0, minCr));
      auto left_in_range = (get<0>(inters) > low) && (get<0>(inters) < upp);
      auto right_in_range = (get<1>(inters) > low) && (get<1>(inters) < upp);
      auto in_range = left_in_range + right_in_range;
      
      if (in_range == 0) {
        // we do not have any interaction in range, so a check is needed to see if
        // the cost or the constraint are greater
        if (minCr <= minS) {
          outcost[k] = quad(tau(costR[j]), low, upp, a(costR[j]), b(costR[j]), c(costR[j])); k++;
        } else {
          outcost[k] = quad(tau(costS[i]), low, upp, a(costS[i]), b(costS[i]), c(costS[i])); k++;
        }
        
      } else if (in_range == 1) {
        // in this case we have only one interaction and we have to understand if it's
        // the right one or the left one
        if (right_in_range) {
          // in this case it's the right
          outcost[k] = quad(tau(costS[i]), low,            get<1>(inters), a(costS[i]), b(costS[i]), c(costS[i])); k++;
          outcost[k] = quad(tau(costR[j]), get<1>(inters), upp,            a(costR[j]), b(costR[j]), c(costR[j])); k++;
        } else {
          // in this case is the left
          outcost[k] = quad(tau(costR[j]), low,            get<0>(inters), a(costR[j]), b(costR[j]), c(costR[j])); k++;
          outcost[k] = quad(tau(costS[i]), get<0>(inters), upp,            a(costS[i]), b(costS[i]), c(costS[i])); k++;
        }
        
      } else if (in_range == 2) {
        // in this case we have two interactions, so first the constraint, then the costf, then the constraint
        outcost[k] = quad(tau(costR[j]), low,            get<0>(inters), a(costR[j]), b(costR[j]), c(costR[j])); k++;
        outcost[k] = quad(tau(costS[i]), get<0>(inters), get<1>(inters), a(costS[i]), b(costS[i]), c(costS[i])); k++;
        outcost[k] = quad(tau(costR[j]), get<1>(inters), upp,            a(costR[j]), b(costR[j]), c(costR[j])); k++;
      } else {
        std::cout << "THIS SHOULD NOT HAPPEN" << '\n';
      }
    }
    
    // updating i, j and l
    if (u(costS[i]) == upp) {
      i++;
    } else {
      j++;
    }
    
    // sliding to the next window
    low = upp;
  }
  
  outcost.resize(k);
  
  // picking the unique values
  std::vector<int> to_select (outcost.size(), true);
  
  for (int k = 1; k < outcost.size(); k++) {
    if (a(outcost[k]) == a(outcost[k - 1])) {
      get<1>(outcost[k]) = get<1>(outcost[k - 1]);
      to_select[k - 1] = false;
    }
  }
  
  std::vector<quad> output;
  for (int i=0; i < outcost.size(); i++) {
    if (to_select[i]) {
      output.push_back(outcost[i]);
    }
  }
  return (output);
}


std::vector<int> backtracking(std::vector<int>& taus) {
  std::vector<int> cp;
  int s = taus.size() + 1;
  //cp.push_back(s);
  
  // for (auto& k:taus) { std::cout << k << ' '; }
  // std::cout << "\n";
  
  // s - 2 since 0 indexing
  while (s != 0) {
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


// This is a O(n) version. It has been tested to work with the simple FPOP
// It might work with the GFPOP, however a proof is needed since at the current state
// the O(n^2) algorithm (still to code) is the best one
std::vector<quad> recomputeIntervals(std::vector<quad>& cost, const double& lower, const double& upper, const double& sigma) {
  
  // initializing an empty vector
  std::vector<quad> outcost (2 * cost.size(), quad(0, 0, 0, 0, 0, 0));
  
  std::vector<int> index(cost.size()); // vector with n ints needed as a reference index for pruning
  iota (std::begin(index), std::end(index), 0); // Fill with 0, 1, ..., n
  
  int i = 0; // index for the current cost function
  int k = 0; // index for the outcost
  double prevInter = -INFINITY; // previous interaction
  
  while (true) {
    
    // if we are at the end of our list
    if ((i + 1) == index.size()) {
      outcost[k] = cost[index[i]];
      get<1>(outcost[k]) = prevInter; // from to the previous interaction
      get<2>(outcost[k]) = INFINITY; // to beyond (actually infinity)
      break; // terminate if we are at the end
      
    } else {
      // getting (the bigger) intersection between the current and the next one
      auto inter = get<1>(getintersections(cost[index[i]], cost[index[i+1]]));
      
      ////// JUST SOME SAFETY CHECKS ON THE INTERSECTION /////
      if (inter < lower || inter > upper) {inter = NAN;} // if out of range make it NAN
      // if there is no valid intersection then quit SHOULD NOT HAPPEN
      if (isnan(inter)) {
        outcost[k] = cost[index[i]];
        get<1>(outcost[k]) = prevInter; // from to the previous interaction
        get<2>(outcost[k]) = INFINITY; // to beyond (actually infinity)
        break; // terminate if we are at the end
      }
      ////////////////////////////////////////////////
      
      if (inter <= prevInter) {
        
        index.erase(index.begin() + i); // erasing the previous quadratic from the index
        // reverting the changes made in the previous step
        i--;
        k--;
        prevInter = l(outcost[k]); // restoring the previous intersection
        
      } else {
        outcost[k] = cost[index[i]];
        get<1>(outcost[k]) = prevInter; // from to the previous interaction
        get<2>(outcost[k]) = inter; // to beyond (actually infinity)
        k++;
        i++;
        prevInter = inter;
      }
    }
  }
  outcost.resize(k + 1);
  
  return outcost;
}



// this function gets what in Rigaill & al. is known as the min less then equal constraint
std::vector<quad> getCostLeq(std::vector<quad>& cost, const int& t) {
  
  // initializing an empty vector
  std::vector<quad> outcost (2 * cost.size(), quad(0, 0, 0, 0, 0, 0));
  
  
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
      outcost[k] = quad(t + 1, get<1>(MINandAT), INFINITY, 0, 0, get<0>(MINandAT));
      break; // terminate
      
    } else {
      
      // this is the candidate line for quadratic i, starting from its minimum
      quad candidateLine = quad(t + 1, get<1>(MINandAT), INFINITY, 0, 0, get<0>(MINandAT));
      
      // checking if there is an intersection with the line and with the next quadtratic
      auto inter = get<1>(getintersections(candidateLine, cost[i + j]));
      auto in_range = (inter > l(cost[i + j])) && (inter < u(cost[i + j]));
      
      // if so we have found an interaction
      if (in_range) {
        
        // add the current quad to the output
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


std::vector<quad> applyl2Penalty(std::vector<quad>& cost, const double& l2penalty, const std::vector<double>& y) {
  // cout << "Applying the transormation!" << endl;
  // applying the tranformation and obtaining a set of disjoint quadratics
  for_each(cost.begin(), cost.end(), [&l2penalty](quad& q){
    auto qold = q; // this is needed to avoid overwriting the current quad
    auto minAt = get<1>(getminimum(qold));
    get<1>(q) = minAt + (l(qold) - minAt) * (a(qold) / l2penalty + 1); // not necessary in case 1
    get<2>(q) = minAt + (u(qold) - minAt) * (a(qold) / l2penalty + 1); // not necessary in case 1
    get<3>(q) = a(qold) * (l2penalty / (a(qold) + l2penalty));
    get<4>(q) = b(qold) * (l2penalty / (a(qold) + l2penalty));
    get<5>(q) =  c(qold) - ( (b(qold) * b(qold)) / (4 * (a(qold) + l2penalty)));
  });
  
  
  // REMOVE OBSERVATIONS OUTISE THE RANGE OF Y
  auto y_min = *min_element(y.begin(), y.end());
  for (int k = 0; k < cost.size(); k++) {
    if (u(cost[k]) < y_min) {
      cost.erase(cost.begin() + k);
    }
  }
  get<1>(cost[0]) = -INFINITY;
  ////////////////////////////////////////////
  
  // recomputing the intervals
  cost = recomputeIntervals(cost, -INFINITY, INFINITY, 0.0);
  
  return cost;
}


// this function generates a Random Walk [to be intended to call from R function dataRW]
std::vector<double> generateAutoRegressive(const double& gamma, const double& y0, const std::vector<double>& mu, const std::vector<double>& ynoise) {
  std::vector<double> y(ynoise.size());
  y[0] = y0 * gamma + mu[0] + ynoise[0];
  for (size_t t = 1; t < (ynoise.size()); t++) {
    y[t] = gamma * (y[t-1] - mu[t - 1]) + mu[t] + ynoise[t];
  }
  return y;
}