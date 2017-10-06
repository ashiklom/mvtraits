#ifndef _DISTRIBUTIONS_
#define _DISTRIBUTIONS_
#include "mvtraits_common.h"

arma::mat c_random_mvnorm(unsigned int n, arma::mat mu, arma::mat Sigma);
arma::mat rwishart(double df, arma::mat S);

#endif
