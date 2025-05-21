#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List emstep_normal(const NumericVector& x,
    const NumericVector& w1,
    const NumericVector& w2,
    const List& params) {

    double mu = as<double>(params["mean"]);
    double sigma = as<double>(params["sd"]);

    int K = x.size();
    if (w1.size() != K || w2.size() != K) {
        stop("x, w1, and w2 must be the same length.");
    }
    
    double ex = 0.0;
    double ex2 = 0.0;
    double en = 0.0;
    
    for (int k = 0; k < K; k++) {
        double z = (x(k) - mu) / sigma;
        double phi = R::dnorm(z, 0.0, 1.0, false);
        double Phi = R::pnorm(z, 0.0, 1.0, true, false);
        double Phi_c = R::pnorm(z, 0.0, 1.0, false, false);
        
        double E_le = mu - sigma * phi / Phi;
        double E_ge = mu + sigma * phi / Phi_c;
        double E_le2 = mu * mu + sigma * sigma - (sigma * sigma * z + 2 * mu * sigma) * phi / Phi;
        double E_ge2 = mu * mu + sigma * sigma + (sigma * sigma * z + 2 * mu * sigma) * phi / Phi_c;
        
        ex += x(k) + w1(k) * E_le + w2(k) * E_ge;
        ex2 += x(k)*x(k) + w1(k) * E_le2 + w2(k) * E_ge2;
        en += 1 + w1(k) + w2(k);
    }

    double mu_new = ex / en;
    double sigma_new = std::sqrt(ex2 / en - mu_new * mu_new);
    
    return List::create(
        Named("mean") = mu_new,
        Named("sd") = sigma_new
    );
}

// [[Rcpp::export]]
List emstep_normal_group(const NumericVector& x,
    const NumericVector& w1,
    const NumericVector& w2,
    const List& params) {

    double mu = as<double>(params["mean"]);
    double sigma = as<double>(params["sd"]);

    int K = x.size();
    if (w1.size() != K || w2.size() != K) {
        stop("x, w1, and w2 must be the same length.");
    }
    
    double ex = 0.0;
    double ex2 = 0.0;
    double en = 0.0;
    
    for (int k = 0; k < K; k++) {
        double z = (x(k) - mu) / sigma;
        double phi = R::dnorm(z, 0.0, 1.0, false);
        double Phi = R::pnorm(z, 0.0, 1.0, true, false);
        double Phi_c = R::pnorm(z, 0.0, 1.0, false, false);
        
        double E_le = mu - sigma * phi / Phi;
        double E_ge = mu + sigma * phi / Phi_c;
        double E_le2 = mu * mu + sigma * sigma - (sigma * sigma * z + 2 * mu * sigma) * phi / Phi;
        double E_ge2 = mu * mu + sigma * sigma + (sigma * sigma * z + 2 * mu * sigma) * phi / Phi_c;
        
        ex += w1(k) * E_le + w2(k) * E_ge;
        ex2 += w1(k) * E_le2 + w2(k) * E_ge2;
        en += w1(k) + w2(k);
    }

    double mu_new = ex / en;
    double sigma_new = std::sqrt(ex2 / en - mu_new * mu_new);
    
    return List::create(
        Named("mean") = mu_new,
        Named("sd") = sigma_new
    );
}