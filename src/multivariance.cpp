#include <Rcpp.h>
using namespace Rcpp;

//' fast euclidean distance matrix computation
//'
//' @param x the vector for which the distanc matrix is computed
//' @examples
//' #require(microbenchmark)
//' #x = rnorm(100)
//' #microbenchmark(fastdist(as.matrix(x)),as.matrix(dist(x)))
//' @export
// [[Rcpp::export]]
NumericMatrix fastdist (const NumericMatrix & x){
  unsigned int outrows = x.nrow(), i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(outrows,outrows);

  for (i = 0; i < outrows - 1; i++){
    Rcpp::NumericVector v1 = x.row(i);
    for (j = i + 1; j < outrows ; j ++){
      d = sqrt(sum(pow(v1-x.row(j), 2.0)));
      out(j,i)=d;
      out(i,j)=d;
    }
  }

  return out;
}

/*** R
*/
