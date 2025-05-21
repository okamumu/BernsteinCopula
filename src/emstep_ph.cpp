#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List emstep_ph_build_interval_data(NumericVector x, NumericVector w1, NumericVector w2) {
  int n = x.size();
  List data(3 * n);
  NumericVector weights(3 * n);

  for (int i = 0; i < n; i++) {
    double xi = x[i];

    // left-censored: [0, xi)
    data[3 * i] = NumericVector::create(0.0, xi);
    weights[3 * i] = w1[i];

    // observed: (xi, xi)
    // Note: This is a point mass at xi, so we can represent it as a closed interval
    // with the same start and end point.
    data[3 * i + 1] = xi;
    weights[3 * i + 1] = 1.0;

    // right-censored: (xi, Inf)
    // Note: This is represented as an open interval starting from xi to infinity.
    // In R, we can use Inf to represent infinity.
    // This is a half-open interval, so we can represent it as [xi, ∞)
    // where the left endpoint is included and the right endpoint is excluded.
    data[3 * i + 2] = NumericVector::create(xi, R_PosInf);
    weights[3 * i + 2] = w2[i];
  }

  return List::create(
    Named("data") = data,
    Named("weights") = weights
  );
}

// [[Rcpp::export]]
List emstep_ph_build_group_data(NumericVector x, NumericVector w1, NumericVector w2) {
  int n = x.size();
  List data(2 * n);
  NumericVector weights(2 * n);

  for (int i = 0; i < n; i++) {
    double xi = x[i];

    // left-censored: [0, xi)
    data[2 * i] = NumericVector::create(0.0, xi);
    weights[2 * i] = w1[i];

    // right-censored: (xi, Inf)
    // Note: This is represented as an open interval starting from xi to infinity.
    // In R, we can use Inf to represent infinity.
    // This is a half-open interval, so we can represent it as [xi, ∞)
    // where the left endpoint is included and the right endpoint is excluded.
    data[2 * i + 1] = NumericVector::create(xi, R_PosInf);
    weights[2 * i + 1] = w2[i];
  }

  return List::create(
    Named("data") = data,
    Named("weights") = weights
  );
}
