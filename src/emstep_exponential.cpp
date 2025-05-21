#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List emstep_exponential(const NumericVector& x,
    const NumericVector& w1,
    const NumericVector& w2,
    const List& params) {

    double lambda = as<double>(params["rate"]);

    int K = x.size();
    if (w1.size() != K || w2.size() != K) {
        stop("x, w1, and w2 must be the same length.");
    }
    
    double ex = 0.0;
    double en = 0.0;
    
    for (int k = 0; k < K; k++) {
        double E_le = 1.0 / lambda - x(k) * exp(-lambda * x(k)) / (1 - exp(-lambda * x(k)));
        double E_ge = x(k) + 1.0 / lambda;
        ex += x(k) + w1(k) * E_le + w2(k) * E_ge;
        en += 1 + w1(k) + w2(k);
    }

    double mu_new = ex / en;
    
    return List::create(
        Named("rate") = mu_new
    );
}

// [[Rcpp::export]]
List emstep_exponential_group(const NumericVector& x,
    const NumericVector& w1,
    const NumericVector& w2,
    const List& params) {

    double lambda = as<double>(params["rate"]);

    int K = x.size();
    if (w1.size() != K || w2.size() != K) {
        stop("x, w1, and w2 must be the same length.");
    }
    
    double ex = 0.0;
    double en = 0.0;
    
    for (int k = 0; k < K; k++) {
        double E_le = 1.0 / lambda - x(k) * exp(-lambda * x(k)) / (1 - exp(-lambda * x(k)));
        double E_ge = x(k) + 1.0 / lambda;
        ex += w1(k) * E_le + w2(k) * E_ge;
        en += w1(k) + w2(k);
    }

    double mu_new = ex / en;
    
    return List::create(
        Named("rate") = mu_new
    );
}