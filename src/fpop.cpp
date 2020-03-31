#include <math.h>
#include<Rcpp.h>
#include <cmath>

#include "algorithms.h"
#include "fpopmain.h"

using namespace Rcpp;
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(.DeCAFS)]]
List DeCAFS(std::vector<double> vectData, double beta, double lambda, double gamma, double phi=0, std::string type="std")
{
  
  auto fpopOut = FPOPmain (vectData, beta, lambda, gamma, phi, type);

  //////// unpacking cost function /////////
  std::list<int> tau;
  std::list<double> l, u, a, b, c;
  for (auto& q:get<2>(fpopOut)) {
    tau.push_back(get<0>(q));
    l.push_back(get<1>(q));
    u.push_back(get<2>(q));
    a.push_back(get<3>(q));
    b.push_back(get<4>(q));
    c.push_back(get<5>(q));
  }

  DataFrame outCost = DataFrame::create(Named("tau") = tau,
                                        Named("l") = l,
                                        Named("u") = u,
                                        Named("a") = a,
                                        Named("b") = b,
                                        Named("c") = c);
  //////////////////////////////////////////


  // returning list
  List res = List::create(
    _["changepoints"] = get<0>(fpopOut),
    _["signal"] = get<1>(fpopOut),
    _["costFunction"] = outCost
  );

  return res;
}

// [[Rcpp::export(.dataAR_c)]]
List dataAR_c(const double& phi, const double& epsilon0, const std::vector<double>& mu, const std::vector<double>& ynoise)
{

  List res = List::create(
    _["z"] = generateAutoRegressive(phi, epsilon0, mu, ynoise)
  );

  return res;
}
