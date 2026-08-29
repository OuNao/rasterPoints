\---

title: 'rasterPoints: Fast Density-Based Bivariate Rasterization and Priority-Based Rendering for Large-Scale Cytometry and Single-Cell Data in R'

tags:

&#x20; - R

&#x20; - C++

&#x20; - OpenMP

&#x20; - flow cytometry

&#x20; - single-cell

&#x20; - visualization

&#x20; - rasterization

authors:

&#x20; - name: Sérgio Costa Fortier

&#x20;   orcid: 0000-0002-4539-6825

&#x20;   affiliation: 1

affiliations:

&#x20; - index: 1

&#x20;   name: 'D''Or Institute for Research and Education (IDOR), São Paulo, Brazil'

date: 29 August 2026

bibliography: paper/paper.bib

\---



\# Summary



Modern single-cell technologies, including multi-parameter flow cytometry, mass cytometry (CyTOF), and single-cell RNA sequencing (scRNA-seq), routinely generate datasets with millions of observations. Bivariate scatter plots are fundamental for exploratory data analysis, population identification, and gating quality control. However, displaying multi-million-point scatter plots using standard R graphics devices introduces severe computational bottlenecks:



1\. \*\*Overplotting and Memory Latency:\*\* Standard vector-based plotting functions (e.g., `graphics::points()`) render individual points as discrete vector objects, leading to high CPU execution latency and excessive memory usage.

2\. \*\*Computational Overhead of Density Estimators:\*\* Traditional density-colored scatter plots rely on two-dimensional kernel density estimation (e.g., `grDevices::densCols()`), which scales poorly as sample sizes exceed $10^6$ events.

3\. \*\*Obscuration of Rare Populations:\*\* In pre-categorized or gated cytometry plots, standard rasterization engines overwrite pixels chronologically, often masking rare, high-value cell subsets (e.g., antigen-specific T cells or stem cells comprising $<1\\%$ of total events) under dominant background populations.



To resolve these limitations, `rasterPoints` provides a high-performance, lightweight R package powered by an optimized C++ backend using `Rcpp` \[@Rcpp] and OpenMP multithreading \[@OpenMP].



\# Statement of Need



Existing high-throughput plotting frameworks in R, such as `scattermore` \[@scattermore] or `datashader` implementations, significantly accelerate point rendering. However, `rasterPoints` fills a critical operational gap in single-cell analysis by offering two unified, memory-lean rasterization paradigms:



\- \*\*Integrated Dynamic Density Rasterization (`data2raster\_density`):\*\* Combines 2D spatial binning, parallel Gaussian neighborhood smoothing, axis saturation prevention via margin bounding, and centered circular pixel dilation into a single C++ pass. Crucially, it incorporates dynamic Min-Max density normalization across populated grid cells, ensuring that pre-gated or dense populations utilize the full color palette spectrum without premature midtone shifts.

\- \*\*Priority-Based Categorical Rendering (`colorder`):\*\* Introduces priority-based pixel evaluation (`colorder`) for pre-categorized populations. Instead of chronological overwriting, the C++ raster engine ensures that higher-priority color indices (representing rare cellular subsets) always overwrite lower-priority background pixels, guaranteeing high-fidelity population preservation (\\autoref{fig:colorder}).



!\[Visual comparison of standard rasterization versus priority-based rendering (`colorder`) on $10^7$ flow cytometry events. \*\*A)\*\* Standard chronological rasterization allows dominant background events (grey) to obscure a rare population (red). \*\*B)\*\* Priority-based rendering ensures that high-priority target events remain crisp and visible regardless of render order or density overlap. \\label{fig:colorder}](paper/Figure1\_colorder\_comparison.png)



\# Computational Performance and Benchmarks



All benchmarks were evaluated on an Intel Core i7-10700 CPU (8 physical cores, 16 threads, 32 GB RAM) running R 4.3.2 on Windows 11 x64.



At $10^7$ events, `rasterPoints` completes full-frame rendering in \*\*612.40 ms\*\* (Default mode), \*\*635.10 ms\*\* (`colorder` mode), and \*\*696.16 ms\*\* (Density mode). This represents up to a \*\*230.2x speedup over base R graphics\*\* (141 seconds) and a \*\*11.9x speedup over `scattermore` integrated with `ggplot2`\*\* (7.29 seconds). Furthermore, enforcing explicit priority ordering (`colorder`) introduces less than \*\*3.7% computational overhead\*\* compared to standard categorical rasterization.



While `scattermore` via `ggplot2` accumulated a peak memory allocation of \*\*12.74 GB\*\* at $10^7$ points, `rasterPoints` capped memory usage at \*\*1.37 GB\*\*—delivering an \*\*\~89% reduction in memory overhead\*\*. Microbenchmarks isolating the `data2raster\_density` C++ engine against the standard `densCols()` pipeline on $5 \\times 10^6$ events demonstrated a \*\*\~9.4x net speedup\*\* (315.37 ms vs 2959.71 ms).



\# Usage Example



The following example demonstrates categorical priority rendering using `colorder` and dynamic kernel density visualization:



```r

library(rasterPoints)



\# Generate synthetic flow cytometry data (10M events)

N <- 1e7

x <- c(rnorm(N \* 0.95, mean = 5, sd = 1), rnorm(N \* 0.05, mean = 8, sd = 0.5))

y <- c(rnorm(N \* 0.95, mean = 5, sd = 1), rnorm(N \* 0.05, mean = 8, sd = 0.5))

mat <- cbind(x, y)



\# 1. High-Performance Density Mode

usr <- c(range(x), range(y))

img\_density <- data2raster\_density(

&#x20; x = mat, usr = usr, width = 800, height = 800,

&#x20; smooth = TRUE, smooth\_radius = 4, smooth\_sigma = 2.0,

&#x20; margin\_pct = 0.05, cex = 1.0, n\_bins = 256

)



\# 2. Priority-Based Categorical Mode (Preserving Rare Populations)

clusters <- c(rep(1, N \* 0.95), rep(2, N \* 0.05)) # Cluster 2 is rare

col\_priority <- c(1, 2) # Give Cluster 2 higher priority

img\_cat <- data2raster(

&#x20; x = mat, col\_idx = clusters, colorder = col\_priority,

&#x20; usr = usr, width = 800, height = 800, cex = rep(1, N)

)

```



\# Acknowledgements



We acknowledge the open-source R and C++ developer communities, particularly the authors and maintainers of `Rcpp`.



\# References

