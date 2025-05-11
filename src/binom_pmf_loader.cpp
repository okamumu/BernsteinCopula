#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dbinom_loader(int n, double p) {
    NumericVector logP(n + 1);
    NumericVector P(n + 1);
    
    if (p <= 0.0) {
        P(0) = 1.0;
        return P;
    }
    if (p >= 1.0) {
        P(n) = 1.0;
        return P;
    }
    
    int k0 = std::floor((n + 1) * p);
    double q = 1.0 - p;
    
    double mu = n * p;
    double sigma = std::sqrt(n * p * q);
    double z = (k0 - mu) / sigma;
    double log_pdf = -0.5 * z * z - std::log(sigma) - 0.5 * std::log(2 * M_PI);
    logP(k0) = log_pdf;
    
    for (int k = k0 + 1; k <= n; ++k) {
        double delta = std::log((n - k + 1.0) / k) + std::log(p / q);
        logP(k) = logP(k - 1) + delta;
    }
    
    for (int k = k0 - 1; k >= 0; --k) {
        double delta = std::log((k + 1.0) / (n - k)) + std::log(q / p);
        logP(k) = logP(k + 1) + delta;
    }
    
    double maxLog = max(logP);
    double sumExp = 0.0;
    for (int k = 0; k <= n; ++k) {
        sumExp += std::exp(logP(k) - maxLog);
    }
    double logZ = maxLog + std::log(sumExp);
    
    for (int k = 0; k <= n; ++k) {
        P(k) = std::exp(logP(k) - logZ);
    }
    
    return P;
}

// [[Rcpp::export]]
void dbinom_loader_inplace(int n, double p, NumericVector& P) {
    P.fill(0.0);
    NumericVector logP(n + 1);
    
    if (p <= 0.0) {
        P(0) = 1.0;
        return;
    }
    if (p >= 1.0) {
        P(n) = 1.0;
        return;
    }
    
    int k0 = std::floor((n + 1) * p);
    double q = 1.0 - p;
    
    double mu = n * p;
    double sigma = std::sqrt(n * p * q);
    double z = (k0 - mu) / sigma;
    double log_pdf = -0.5 * z * z - std::log(sigma) - 0.5 * std::log(2 * M_PI);
    logP(k0) = log_pdf;
    
    for (int k = k0 + 1; k <= n; ++k) {
        double delta = std::log((n - k + 1.0) / k) + std::log(p / q);
        logP(k) = logP(k - 1) + delta;
    }
    
    for (int k = k0 - 1; k >= 0; --k) {
        double delta = std::log((k + 1.0) / (n - k)) + std::log(q / p);
        logP(k) = logP(k + 1) + delta;
    }
    
    double maxLog = max(logP);
    double sumExp = 0.0;
    for (int k = 0; k <= n; ++k) {
        sumExp += std::exp(logP(k) - maxLog);
    }
    double logZ = maxLog + std::log(sumExp);
    
    for (int k = 0; k <= n; ++k) {
        P(k) = std::exp(logP(k) - logZ);
    }
}
