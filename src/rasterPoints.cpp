#include <Rcpp.h>
#include <omp.h>
#include <cmath>

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
//' @param smooth Logical, Apply Gaussian neighborhood kernel smoothing in density mode?
//' @param smooth_radius Integer, Kernel neighborhood radius in pixels. Default is 4L.
//' @param smooth_sigma Double, Gaussian standard deviation. Default is 2.0.
//' @param margin_pct Double, interior bounding box margin ratio (0.05). Prevents axis saturation artifacts.
//' @param cex Double, point size scaling parameter.
//' @param n_bins Integer, number of color steps (e.g., 256).
//' @param ncores Integer, number of OpenMP threads.
// [[Rcpp::export]]
IntegerMatrix data2raster_density(const NumericMatrix& x, 
                                  NumericVector usr, 
                                  int width, int height, 
                                  int n_bins = 256,
                                  bool smooth = false,
                                  int smooth_radius = 4,
                                  double smooth_sigma = 2.0,
                                  double margin_pct = 0.05,
                                  double cex = 1.0,
                                  int ncores = 0) {
  int n = x.nrow();
  std::vector<int> raw_counts(width * height, 0);
  
  double scale_x = (double)width / (usr[1] - usr[0]);
  double scale_y = (double)height / (usr[3] - usr[2]);
  
  // ---------------------------------------------------------------------------
  // STEP 1: 2D Binning / Raw Point Aggregation (Exact 1-pixel mapping)
  // ---------------------------------------------------------------------------
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < n; i++) {
    int px = (int)((x(i, 0) - usr[0]) * scale_x);
    int py = (int)((x(i, 1) - usr[2]) * scale_y);
    
    int r = (height - 1) - py;
    int c = px;
    
    if (r >= 0 && c >= 0 && r < height && c < width) {
#pragma omp atomic
      raw_counts[r * width + c]++;
    }
  }
  
  // ---------------------------------------------------------------------------
  // STEP 2: Optional Gaussian Neighborhood Density Smoothing
  // ---------------------------------------------------------------------------
  std::vector<double> neighborhood_density(width * height, 0.0);
  
  if (!smooth) {
    for (size_t i = 0; i < raw_counts.size(); i++) {
      neighborhood_density[i] = (double)raw_counts[i];
    }
  } else {
    std::vector<double> kernel;
    int k_size = 2 * smooth_radius + 1;
    kernel.reserve(k_size * k_size);
    double two_sigma2 = 2.0 * smooth_sigma * smooth_sigma;
    
    for (int dr = -smooth_radius; dr <= smooth_radius; dr++) {
      for (int dc = -smooth_radius; dc <= smooth_radius; dc++) {
        double dist2 = (double)(dr * dr + dc * dc);
        kernel.push_back(std::exp(-dist2 / two_sigma2));
      }
    }
    
#pragma omp parallel for num_threads(ncores)
    for (int r = 0; r < height; r++) {
      for (int c = 0; c < width; c++) {
        if (raw_counts[r * width + c] == 0) continue;
        
        double sum_density = 0.0;
        int k_idx = 0;
        for (int dr = -smooth_radius; dr <= smooth_radius; dr++) {
          int nr = r + dr;
          if (nr < 0 || nr >= height) {
            k_idx += k_size;
            continue;
          }
          for (int dc = -smooth_radius; dc <= smooth_radius; dc++) {
            int nc = c + dc;
            if (nc >= 0 && nc < width) {
              sum_density += raw_counts[nr * width + nc] * kernel[k_idx];
            }
            k_idx++;
          }
        }
        neighborhood_density[r * width + c] = sum_density;
      }
    }
  }
  
  // ---------------------------------------------------------------------------
  // STEP 3: Robust Maximum Search & Linear Mapping with Centered Circle Dilators
  // ---------------------------------------------------------------------------
  int margin_x = (int)(width * margin_pct);
  int margin_y = (int)(height * margin_pct);
  
  if (margin_x >= width / 2) margin_x = 0;
  if (margin_y >= height / 2) margin_y = 0;
  
  double max_d = 0.0;
  
  for (int r = margin_y; r < height - margin_y; r++) {
    for (int c = margin_x; c < width - margin_x; c++) {
      double d = neighborhood_density[r * width + c];
      if (d > max_d) max_d = d;
    }
  }
  
  if (max_d == 0.0) {
    for (size_t i = 0; i < neighborhood_density.size(); i++) {
      if (neighborhood_density[i] > max_d) max_d = neighborhood_density[i];
    }
  }
  
  IntegerMatrix res(height, width);
  if (max_d == 0.0) return res;
  
  int cex_int = (int)std::round(cex);
  if (cex_int < 1) cex_int = 1;
  
#pragma omp parallel for num_threads(ncores)
  for (int r = 0; r < height; r++) {
    for (int c = 0; c < width; c++) {
      if (raw_counts[r * width + c] > 0) {
        double d_val = neighborhood_density[r * width + c];
        if (d_val > max_d) d_val = max_d;
        
        double norm = d_val / max_d;
        int idx = (int)(norm * (n_bins - 1)) + 1;
        if (idx > n_bins) idx = n_bins;
        
        if (cex_int == 1) {
          res(r, c) = idx;
        } else {
          // Symmetric pixel offsets for exact centering
          int r_offset_back = (cex_int - 1) / 2;
          int r_offset_fwd  = cex_int / 2;
          
          int r_min = std::max(0, r - r_offset_back);
          int r_max = std::min(height - 1, r + r_offset_fwd);
          int c_min = std::max(0, c - r_offset_back);
          int c_max = std::min(width - 1, c + r_offset_fwd);
          
          // Squared radius threshold for circular masking
          double max_r2 = (cex_int / 2.0) * (cex_int / 2.0);
          
          for (int dr = r_min; dr <= r_max; dr++) {
            for (int dc = c_min; dc <= c_max; dc++) {
              
              // Circular condition active for all cex >= 2
              double dist_r = dr - r;
              double dist_c = dc - c;
              if ((dist_r * dist_r + dist_c * dist_c) > max_r2) continue;
              
              // Retain higher density color index on visual overlap
              if (idx > res(dr, dc)) {
                res(dr, dc) = idx;
              }
            }
          }
        }
      }
    }
  }
  return res;
}
