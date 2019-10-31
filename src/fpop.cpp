#include <math.h>
#include<Rcpp.h>
#include <cmath>

#include "algorithms.h"
#include "fpopmain.h"

using namespace Rcpp;
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List l2fpop(std::vector<double> vectData, double beta, double lambda, double gamma, double phi=0, std::string type="std")
{

  List res = List::create(
    _["changepoints"] = FPOPmain (vectData, beta, lambda, gamma, phi, type)
  );

  return res;
}

// [[Rcpp::export]]
List dataAR_c(const double& gamma, const double& y0, const std::vector<double>& mu, const std::vector<double>& ynoise)
{
  
  List res = List::create(
    _["z"] = generateAutoRegressive(gamma, y0, mu, ynoise)
  );
  
  return res;
}