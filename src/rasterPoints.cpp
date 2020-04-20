#include <Rcpp.h>
#include <omp.h>
#include <math.h>

using namespace Rcpp;

//' data2raster
//' 
//' This function convert a 2 column matrix of data in a img matrix.
//'
//' @param x matrix of dim = height X 2.
//' @param col An character vector with color codes with length = nrow(x).
//' @param usr An double vector. See par("usr").
//' @param width An integer.
//' @param height An integer.
//' @param cex An integer.
//' @param ncores Number of cores to use. Default to 0 (max OpenMP avaiable cores).
// [[Rcpp::export]]
CharacterMatrix data2raster(NumericMatrix x, CharacterVector col, NumericVector usr, int width, int height, int cex=1, int ncores=0) {
  NumericMatrix data=clone(x);
  data(_,0)=((data(_,0)-usr[0])/(usr[1]-usr[0]))*width;
  data(_,1)=((data(_,1)-usr[2])/(usr[3]-usr[2]))*height;
  int n=data.nrow();
  CharacterMatrix res(height, width);
  std::fill(res.begin(), res.end(), NA_STRING) ;
  int nthreads=omp_get_max_threads();
  if (cex < 1) cex = 1;
  if ((ncores<=0) || (ncores>nthreads)) ncores=nthreads;
  int minidx = -(int)(((float)(cex)-0.5)/2);
  int maxidx = (int)(cex/2);
  float center = 0.5;
  if (cex & 1) center = 0;
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < n; i++) {
    int r=(height-1)-(int)data(i,1);
    int c=(int)data(i,0);
    for (int j = minidx; j <= maxidx; j++) {
      for (int k = minidx; k <= maxidx; k++) {
        int r2 = r+j;
        int c2 = c+k;
        if ((r2<0) | (c2<0) | (r2>=height) | (c2>=width)) continue;
        float distj = fabsf((float)j-center) + 0.5;
        float distk = fabsf((float)k-center) + 0.5;
        if (sqrt(pow(distj,2) + pow(distk,2)) > ((float)(cex)/2)+0.5) continue;
        res(r2, c2)=(char*)col(i);
      }
    }
  }
  return(res);
}
