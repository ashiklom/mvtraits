#include "mvtraits_zigg.h"

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
void zsetseed(unsigned long int s) {
    zigg.setSeed(s);
    return;
}

// [[Rcpp::export]]
unsigned long int zgetseed() {
    return zigg.getSeed();
}

// [[Rcpp::export]]
double zrnorm() {
    return zigg.norm();
}
