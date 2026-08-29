# rasterPoints

[![R-CMD-check](https://github.com/OuNao/rasterPoints/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/OuNao/rasterPoints/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`rasterPoints` provides high-performance, parallelized C++ rasterization engines for massive single-cell and flow cytometry datasets in R. It delivers fast rendering of up to 10 million events with priority-based categorical display and dynamically normalized density plotting.

## Features

- Priority-Based Categorical Rendering (`data2raster`): Uses the `colorder` priority system to prevent rare cellular subsets from being obscured by high-density background populations.
- Dynamic Density Rasterization (`data2raster_density`): Combines 2D binning, parallel Gaussian smoothing, margin-based axis saturation guard, and full Min-Max normalization to prevent color palette truncation.
- High Performance: Multi-threaded execution via OpenMP and Rcpp, offering up to 230x speedup over base R graphics with minimal memory footprint (~1.37 GB peak for 10M events).

## Installation

You can install the development version of rasterPoints from GitHub:

```r
# install.packages("devtools")
devtools::install_github("OuNao/rasterPoints")
```

## Quick Start

```r
library(rasterPoints)

# Generate synthetic high-density data (1M events)
N <- 1e6
x <- rnorm(N, mean = 5, sd = 1)
y <- rnorm(N, mean = 5, sd = 1)
mat <- cbind(x, y)
usr <- c(range(x), range(y))

# 1. Density Rasterization
img_density <- data2raster_density(
  x = mat, usr = usr, width = 800, height = 800,
  smooth = TRUE, smooth_radius = 4, smooth_sigma = 2.0,
  margin_pct = 0.05, cex = 1.0, n_bins = 256
)

# Render raster in base R
plot(NA, xlim = usr[1:2], ylim = usr[3:4], xlab = "FSC-A", ylab = "SSC-A")
rasterImage(as.raster(matrix(terrain.colors(256)[img_density], 800, 800)), 
            usr[1], usr[3], usr[2], usr[4])
```

## Citation

If you use rasterPoints in your research, please cite our paper:

- Fortier, S. C. (2026). rasterPoints: Fast Density-Based Bivariate Rasterization and Priority-Based Rendering for Large-Scale Cytometry and Single-Cell Data in R. Journal of Open Source Software (Pending).

## License

This project is licensed under the GPL-3 License.