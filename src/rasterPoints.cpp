#include <Rcpp.h>
#include <omp.h>

using namespace Rcpp;

//' data2raster
//' 
//' This function convert a 2 column matrix of data in a img matrix.
//'
//' @param x Data matrix (ncol = 2).
//' @param col_idx Integer vector with color indices (1-based from R). With length = nrow(x).
//' @param colorder Numeric vector where index is color_id and value is priority.
//' @param usr An double vector. See par("usr").
//' @param width Output raster width.
//' @param height Output raster height.
//' @param cex Point size vector.
//' @param ncores Number of cores to use. Default to 0 (use max OpenMP avaiable cores).
// [[Rcpp::export]]
IntegerMatrix data2raster(const NumericMatrix& x, 
                              const IntegerVector& col_idx, 
                              const IntegerVector& colorder, 
                              NumericVector usr, 
                              int width, int height, 
                              const IntegerVector cex, 
                              int ncores = 0) {
  
  int n = x.nrow();
  IntegerMatrix res(height, width); 
  IntegerMatrix p_map(height, width); 
  
  // Initialize matrices: 0 for empty color, -1 for lowest priority
  std::fill(res.begin(), res.end(), 0);
  std::fill(p_map.begin(), p_map.end(), -1);

  int nthreads = omp_get_max_threads();
  if (ncores <= 0 || ncores > nthreads) ncores = nthreads;

  // Pre-calculate scales to move divisions out of the 5M loop
  double x_range = usr[1] - usr[0];
  double y_range = usr[3] - usr[2];
  double scale_x = (width - 1) / x_range;
  double scale_y = (height - 1) / y_range;
  
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < n; i++) {
    int px = (int)((x(i, 0) - usr[0]) * scale_x);
    int py = (int)((x(i, 1) - usr[2]) * scale_y);
    
    int r = (height - 2) - py;
    int c = px - 1;
    
    int cexp = cex[i];
    float radius_sq = pow((float)cexp / 2.0 + 0.5, 2);
    float center = (cexp & 1) ? 0.0 : 0.5;
    
    int minidx = -(int)(((float)cexp - 0.5) / 2.0);
    int maxidx = (int)(cexp / 2);
    
    int current_cor = col_idx[i];
    int current_prio = colorder[current_cor - 1]; // Maps color ID to its priority
    
    for (int j = minidx; j <= maxidx; j++) {
      for (int k = minidx; k <= maxidx; k++) {
        int r2 = r + j;
        int c2 = c + k;
        
        if (r2 < 0 || c2 < 0 || r2 >= height || c2 >= width) continue;
        
        // Distance check using squared values (O(1) vs O(N) of sqrt)
        float dist_sq = pow(fabsf((float)j - center) + 0.5, 2) + pow(fabsf((float)k - center) + 0.5, 2);
        
        if (dist_sq > radius_sq) continue;
        
        // Priority logic: only overwrite if current point has higher or equal priority.
        // Using integers makes this comparison extremely fast and safer for concurrent writes.
        if (current_prio >= p_map(r2, c2)) {
          p_map(r2, c2) = current_prio;
          res(r2, c2) = current_cor;
        }
      }
    }
  }
  return res;
}
