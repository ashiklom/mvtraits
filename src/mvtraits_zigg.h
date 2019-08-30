#ifndef _ZIGG_
#define _ZIGG_
#include "mvtraits_common.h"
#include <Ziggurat.h>

// [[Rcpp::depends(RcppZiggurat)]]

void zset_seed(unsigned long int s);
unsigned long int zget_seed();
double zrnorm();

#endif
