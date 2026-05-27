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

#include <Rcpp.h>
#include <omp.h>
#include <math.h>

using namespace Rcpp;

//' data2raster_density
//'
//' High-performance density aggregation for large scale flow cytometry data.
//' This function uses a binning approach (Datashader-style) instead of KDE
//' to achieve near-instant rendering of millions of events.
//'
//' @param x NumericMatrix with 2 columns (fluorescence channels).
//' @param usr NumericVector with plot limits c(x1, x2, y1, y2).
//' @param width Integer, output image width in pixels.
//' @param height Integer, output image height in pixels.
//' @param n_bins Integer, number of color steps (e.g., 256).
//' @param ncores Integer, number of OpenMP threads.
// [[Rcpp::export]]
IntegerMatrix data2raster_density(const NumericMatrix& x, 
                                  NumericVector usr, 
                                  int width, int height, 
                                  int n_bins = 256,
                                  int ncores = 0) {
  
  int n = x.nrow();
  
  // Internal buffer to store the hit count per pixel (density map)
  IntegerMatrix dens_map(height, width);
  std::fill(dens_map.begin(), dens_map.end(), 0);
  
  // Pre-calculate projection scales to move divisions out of the 5M loop
  double x_range = usr[1] - usr[0];
  double y_range = usr[3] - usr[2];
  double scale_x = (width - 1) / x_range;
  double scale_y = (height - 1) / y_range;
  
  int nthreads = omp_get_max_threads();
  if (ncores <= 0 || ncores > nthreads) ncores = nthreads;
  
  // LOOP 1: Aggregation (Binning)
  // We process 5M points and increment the count at the corresponding pixel.
  // Using #pragma omp atomic for thread-safe increments on the shared matrix.
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < n; i++) {
    int px = (int)((x(i, 0) - usr[0]) * scale_x);
    int py = (int)((x(i, 1) - usr[2]) * scale_y);
    
    // Invert Y axis to match R's raster coordinate system (top-down)
    int r = (height - 1) - py;
    int c = px;
    
    // Bounds check before writing to memory
    if (r >= 0 && c >= 0 && r < height && c < width) {
#pragma omp atomic
      dens_map(r, c)++;
    }
  }
  
  // Find the global maximum density to normalize colors
  int max_d = 0;
  for (int i = 0; i < dens_map.length(); i++) {
    if (dens_map[i] > max_d) max_d = dens_map[i];
  }
  
  // Final raster matrix containing color indices
  IntegerMatrix res(height, width);
  if (max_d == 0) return res;
  
  // LOOP 2: Normalization and Colormap mapping
  // We only iterate over the screen pixels (e.g., 800x600), making this very fast.
#pragma omp parallel for num_threads(ncores)
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      int count = dens_map(r, c);
      if (count > 0) {
        // Log-normalization: Essential for flow cytometry to visualize 
        // low-density populations alongside high-density cores.
        double norm = log1p(count) / log1p(max_d);
        
        // Map to 1-based index for R palette accessibility
        res(r, c) = (int)(norm * (n_bins - 1)) + 1;
      } else {
        // Zero represents transparency (NA) in the R wrapper
        res(r, c) = 0;
      }
    }
  }
  
  return res;
}
