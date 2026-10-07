#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double h_AMISE_cpp(NumericVector x) {
  int n = x.size();
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  double a = 0.25 * (n - 1), b = 0.75 * (n - 1);
  int i = a, j = b;
  a = (1 - (a - i)) * y[i] + (a - i) * y[i + 1];
  b = (1 - (b - j)) * y[j] + (b - j) * y[j + 1];
  double s = sd(x);
  s = std::min(s, (b - a) / 1.34);
  return std::pow(40 * std::sqrt(std::acos(-1.0)), 0.2) * s * std::pow(n, -0.2);
}

NumericVector cv_cpp(double h, Function f, IntegerVector g, List z, NumericVector x) {
  int n = x.size();
  NumericVector v(n);
  for (int i = 0; i < n; i++) {
    v[i] = as<double>(f(x[i], z[g[i] - 1], h));
  }
  return v;
}

// [[Rcpp::export]]
double nll_cpp(double h, Function f, IntegerVector g, List z, NumericVector x) {
  NumericVector v = cv_cpp(h, f, g, z, x);
  long double s = 0;
  for (int i = 0; i < v.size(); i++) {
    s += std::log(v[i]);
  }
  return -s;
}

// [[Rcpp::export]]
double ucv_cpp(double h, Function f, Function s, IntegerVector g, List z, NumericVector x) {
  NumericVector v = cv_cpp(h, f, g, z, x);
  long double t = 0;
  for (int i = 0; i < v.size(); i++) {
    t += v[i];
  }
  return as<double>(s(_["x"] = x, _["h"] = h)) - 2.0 / x.size() * static_cast<double>(t);
}
