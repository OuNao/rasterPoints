---
title: 'rasterPoints: Fast Density-Based Bivariate Rasterization and Priority-Based Rendering for Large-Scale Cytometry and Single-Cell Data in R'
tags:
  - R
  - C++
  - OpenMP
  - flow cytometry
  - single-cell
  - visualization
  - rasterization
authors:
  - name: Sérgio Costa Fortier
    orcid: 0000-0002-4539-6825
    affiliation: '1'
affiliations:
  - index: 1
    name: D’Or Institute for Research and Education (IDOR), São Paulo, Brazil
date: 29 August 2026
bibliography: paper.bib
---

# Summary

Modern single-cell technologies, including multi-parameter flow cytometry, mass cytometry (CyTOF), and single-cell RNA sequencing (scRNA-seq), routinely generate datasets with millions of observations. `rasterPoints` is an R package designed to solve high-throughput bivariate visualization challenges for biological datasets containing up to 10 million events. Powered by an optimized C++ backend (`Rcpp` and OpenMP multithreading), `rasterPoints` provides memory-lean, parallelized engines for dynamic density estimation and priority-based categorical display, enabling rapid data exploration directly within base R graphics pipelines.

# Statement of need

Bivariate scatter plots are fundamental for exploratory data analysis, population identification, and gating quality control in single-cell biology. However, displaying multi-million-point scatter plots using standard R graphics devices introduces severe operational bottlenecks:

1. **Overplotting and Memory Latency:** Standard vector-based plotting functions (e.g., `graphics::points()`) render individual points as discrete vector objects, leading to high CPU execution latency and excessive memory usage.
2. **Computational Overhead of Density Estimators:** Traditional density-colored scatter plots rely on two-dimensional kernel density estimation (e.g., `grDevices::densCols()`), which scales poorly as sample sizes exceed $10^6$ events.
3. **Obscuration of Rare Populations:** In pre-categorized or gated cytometry plots, standard rasterization engines overwrite pixels chronologically, often masking rare, high-value cell subsets (e.g., antigen-specific T cells or stem cells comprising $<1\%$ of total events) under dominant background populations.

`rasterPoints` targets computational biologists, cytometrists, and single-cell bioinformaticians requiring fast, high-fidelity visualization pipelines without heavy framework dependencies.

# State of the field

Existing high-throughput plotting frameworks in R, such as `scattermore` [@scattermore] or `datashader` implementations, significantly accelerate point rendering. However, existing packages focus predominantly on standard pixel blending or simple point density aggregation without addressing the domain-specific challenges of single-cell cytometry.

Rather than contributing incremental patches to general-purpose scatter libraries, a dedicated build was justified to implement two domain-tailored rasterization algorithms:
- **Dynamic Density Rasterization (`data2raster_density`):** Integrates 2D spatial binning, parallel Gaussian neighborhood smoothing, axis saturation prevention via margin bounding, and centered circular pixel dilation into a single C++ pass with dynamic Min-Max density normalization across populated grid cells.
- **Priority-Based Categorical Rendering (`colorder`):** Evaluates pre-categorized populations by explicit priority tiers rather than chronological drawing order, ensuring that high-value rare subsets are preserved during pixel projection (\autoref{fig:colorder}).

![Visual comparison of standard rasterization versus priority-based rendering (`colorder`) on $10^7$ flow cytometry events. **A)** Standard chronological rasterization allows dominant background events (grey) to obscure a rare population (red). **B)** Priority-based rendering ensures that high-priority target events remain crisp and visible regardless of render order or density overlap. \label{fig:colorder}](Figure1_colorder_comparison.png)

# Software design

`rasterPoints` prioritizes computational performance, low memory latency, and API simplicity. The package architecture delegates heavy matrix transformations and pixel mapping directly to C++ via `Rcpp` [@Rcpp], leveraging OpenMP multithreading [@OpenMP] for parallelized binning and Gaussian smoothing.

Design trade-offs were made to balance resolution against rendering speed:
1. **Grid Bounding vs. Vector Coordinates:** The engine maps continuous float coordinate matrices into discrete $N \times N$ integer grids (default 256–1024 bins), decoupling rendering latency from raw event counts.
2. **In-Memory Bitmaps vs. Scene Graphs:** Rather than producing intermediate ggplot/graphics objects, the C++ core outputs native integer matrix bitmaps directly compatible with base R `rasterImage()`.
3. **Priority Evaluation Overhead:** Enforcing `colorder` pixel checks requires conditional array evaluation per point, but multi-threaded C++ execution keeps this overhead negligible compared to standard non-prioritized rasterization.

# Research impact statement

The efficiency of `rasterPoints` has been demonstrated through performance benchmarks and real-world deployment on an Intel Core i7-10700 CPU (8 physical cores, 16 threads, 32 GB RAM) running R 4.3.2 on Windows 11 x64.

At $10^7$ events, `rasterPoints` completes full-frame rendering in **2.10 s** (Default mode), **2.01 s** (`colorder` mode), and **696.16 ms** (Density mode). This represents up to a **202.5x speedup over base R graphics** (141 seconds) and up to a **10.5x speedup over `scattermore` integrated with `ggplot2`** (7.29 seconds). While `scattermore` via `ggplot2` accumulated a peak memory allocation of **12.74 GB** at $10^7$ points, `rasterPoints` capped memory usage at **1.37 GB**—delivering an **~89% reduction in memory overhead**. Microbenchmarks isolating the `data2raster_density` C++ engine against the standard `densCols()` pipeline on $5 \times 10^6$ events demonstrated a **~9.4x net speedup** (315.37 ms vs 2959.71 ms).
Furthermore, `rasterPoints` powers the visualization engine behind [FlowDraw](https://www.flowdraw.com.br), an online platform for interactive single-cell and flow cytometry data analysis, serving as concrete proof of production-grade stability and real-world research utility.

```r
library(rasterPoints)

# Generate synthetic flow cytometry data (10M events)
N <- 1e7
x <- c(rnorm(N * 0.95, mean = 5, sd = 1), rnorm(N * 0.05, mean = 8, sd = 0.5))
y <- c(rnorm(N * 0.95, mean = 5, sd = 1), rnorm(N * 0.05, mean = 8, sd = 0.5))
mat <- cbind(x, y)
usr <- c(range(x), range(y))

# 1. High-Performance Density Mode
rasterPoints(
  x = mat, usr = usr, width = 800, height = 800,
  type = "density", smooth = TRUE, smooth_radius = 4, smooth_sigma = 2.0,
  margin_pct = 0.05, cex = 1.0, n_bins = 256
)

# 2. Priority-Based Categorical Mode (Preserving Rare Populations)
clusters <- c(rep(1, N * 0.95), rep(2, N * 0.05)) # Cluster 2 is rare
col_priority <- c(1, 2) # Give Cluster 2 higher priority
rasterPoints(
  x = mat, col_idx = clusters, colorder = col_priority,
  usr = usr, width = 800, height = 800, cex = rep(1, N)
)
```

# Acknowledgements

We acknowledge the open-source R and C++ developer communities, particularly the authors and maintainers of `Rcpp`.

# References