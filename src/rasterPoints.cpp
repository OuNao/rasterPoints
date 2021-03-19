#include <Rcpp.h>
#include <omp.h>
#include <math.h>

using namespace Rcpp;

//' data2raster
//' 
//' This function convert a 2 column matrix of data in a img matrix.
//'
//' @param x Data matrix (ncol = 2).
//' @param col An character vector with color codes with length = nrow(x).
//' @param usr An double vector. See par("usr").
//' @param width An integer.
//' @param height An integer.
//' @param cex An integer.
//' @param ncores Number of cores to use. Default to 0 (use max OpenMP avaiable cores).
//' @param colorder Named numeric vector color priority..
// [[Rcpp::export]]
CharacterMatrix data2raster(NumericMatrix x, CharacterVector col, NumericVector colorder, NumericVector usr, int width, int height, int cex=1, int ncores=0) {
  NumericMatrix data=clone(x);
  data(_,0)=((data(_,0)-usr[0])/(usr[1]-usr[0]))*(width-1);
  data(_,1)=((data(_,1)-usr[2])/(usr[3]-usr[2]))*(height-1);
  int n=data.nrow();
  CharacterVector colv(n);
  if (col.length() == 1) colv.fill(col(0)); else colv=col;
  CharacterMatrix res(height, width);
  std::fill(res.begin(), res.end(), NA_STRING) ;
  if (cex < 1) cex = 1;
  int minidx = -(int)(((float)(cex)-0.5)/2);
  int maxidx = (int)(cex/2);
  float center = 0.5;
  if (cex & 1) center = 0;
  int nthreads=omp_get_max_threads();
  if ((ncores<=0) || (ncores>nthreads)) ncores=nthreads;
  if (colorder(0) != 0) ncores=1;
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < n; i++) {
    int r=(height-2)-(int)round(data(i,1));
    int c=(int)round(data(i,0))-1;
    for (int j = minidx; j <= maxidx; j++) {
      for (int k = minidx; k <= maxidx; k++) {
        int r2 = r+j;
        int c2 = c+k;
        if ((r2<0) | (c2<0) | (r2>=height) | (c2>=width)) continue;
        float distj = fabsf((float)j-center) + 0.5;
        float distk = fabsf((float)k-center) + 0.5;
        if (sqrt(pow(distj,2) + pow(distk,2)) > ((float)(cex)/2)+0.5) continue;
        if (colorder(0) == 0 && res(r2, c2) != NA_STRING) continue;
        if (colorder(0) != 0 && res(r2, c2) != NA_STRING && (int)colorder[(char*)res(r2, c2)] >= (int)colorder[(char*)colv(i)]) continue;
        res(r2, c2)=(char*)colv(i);
      }
    }
  }
  return(res);
}
