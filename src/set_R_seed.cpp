#include "mvtraits_zigg.h"

Ziggurat::Ziggurat::Ziggurat zigg;

//[[Rcpp::export]]
void set_R_seed() {
    Rcpp::Environment glob = Rcpp::Environment::global_env();
    Rcpp::NumericVector rseed_vec = glob[".Random.seed"];
    unsigned int index = 2 + R::runif(0, 1) * (rseed_vec.size() - 2);
    unsigned long rseed = rseed_vec[index];
    Rcpp::Rcout << "Using seed: " << rseed << std::endl;
    zigg.setSeed(rseed);
    return;
}
