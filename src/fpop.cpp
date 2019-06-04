#include <math.h>
#include<Rcpp.h>
#include <cmath>

#include "fpopmain.h"

using namespace Rcpp;
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List l2fpop(std::vector<double> vectData, double l0penalty, double l2penalty, std::string type="std")
{


  /// RETURN
  List res = List::create(
    _["changepoints"] = FPOPmain(vectData, l0penalty, l2penalty, type)
  );

  return res;
}
