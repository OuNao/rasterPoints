# ==============================================================================
# COMPLETE BENCHMARK AND VALIDATION SCRIPT: RASTERPOINTS PACKAGE
# ==============================================================================
# This script evaluates the rendering performance (execution time and memory) 
# and OpenMP scalability of the rasterPoints package compared to traditional 
# methods (Base R) and the scattermore library (Base R and ggplot2).
# ==============================================================================

library(rasterPoints)
library(scattermore)
library(ggplot2)
library(bench)
library(parallel)

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTION: SYNTHETIC FCS DATA GENERATION (DENSITY + RARE POPULATION)
# ------------------------------------------------------------------------------
generate_fcs_data <- function(n_points = 1e6) {
  set.seed(42)
  
  # Main background population (95% of points)
  n_bg <- round(n_points * 0.95)
  x_bg <- rnorm(n_bg, mean = 5, sd = 1.5)
  y_bg <- rnorm(n_bg, mean = 5, sd = 1.5)
  
  # Rare target population of interest (5% of points)
  n_rare <- n_points - n_bg
  x_rare <- rnorm(n_rare, mean = 8, sd = 0.3)
  y_rare <- rnorm(n_rare, mean = 8, sd = 0.3)
  
  # Combine coordinates
  coords <- rbind(cbind(x_bg, y_bg), cbind(x_rare, y_rare))
  colnames(coords) <- c("V1", "V2")
  
  # Color vector (Gray for background, Red for rare population)
  colors <- c(rep("gray70", n_bg), rep("red3", n_rare))
  random_order <- sample(1:length(colors))
  
  return(list(coords = coords[random_order, ], colors = colors[random_order]))
}

# ------------------------------------------------------------------------------
# 2. ENVIRONMENT AND HARDWARE SPECIFICATION CAPTURE
# ------------------------------------------------------------------------------
get_system_specs <- function() {
  cat("==================================================================\n")
  cat("               ENVIRONMENT AND TEST SPECIFICATIONS                \n")
  cat("==================================================================\n")
  cat("OS:              ", sessionInfo()$running, "\n")
  cat("R Version:       ", R.version$version.string, "\n")
  cat("Platform:        ", R.version$platform, "\n")
  cat("Physical Cores:  ", parallel::detectCores(logical = FALSE), "\n")
  cat("Logical Cores:   ", parallel::detectCores(logical = TRUE), "\n")
  cat("==================================================================\n\n")
}

# ------------------------------------------------------------------------------
# 3. RENDERING AND MEMORY BENCHMARK
# ------------------------------------------------------------------------------
run_rendering_benchmark <- function(sample_sizes = c(1e5, 1e6, 5e6, 1e7)) {
  results_list <- list()
  
  for (n in sample_sizes) {
    message(sprintf("Running Benchmark for N = %s points...", format(n, scientific = FALSE)))
    fcs <- generate_fcs_data(n)
    data <- fcs$coords
    colors <- fcs$colors
    
    res <- bench::mark(
      # 1. Traditional Base R (Vector points)
      Base_R = {
        png(tempfile(), width = 800, height = 800)
        plot(data, type = "n")
        points(data, pch = ".", col = colors)
        dev.off()
      },
      
      # 2. Scattermore in Base R (Without ggplot2 overhead)
      Scattermore_Base = {
        png(tempfile(), width = 800, height = 800)
        scattermore::scattermoreplot(data[, 1], data[, 2], col = colors)
        dev.off()
      },
      
      # 3. Scattermore with ggplot2
      Scattermore_ggplot = {
        df <- as.data.frame(data)
        p <- ggplot(df, aes(x = V1, y = V2)) + 
          scattermore::geom_scattermore(color = colors)
        ggsave(tempfile(fileext = ".png"), plot = p, width = 800, height = 800, units = "px")
      },
      
      # 4. rasterPoints Default (No colorder, no density, no interpolation)
      Rasterpoints_Default = {
        png(tempfile(), width = 800, height = 800)
        rasterPoints(data, col = colors, colorder = "default", 
                     cex = 1, interpolate = FALSE, force_new = TRUE)
        dev.off()
      },
      
      # 5. rasterPoints with Priority-Based Rendering (Colorder without interpolation)
      Rasterpoints_Colorder = {
        png(tempfile(), width = 800, height = 800)
        rasterPoints(data, col = colors, colorder = c("gray70", "red3"), 
                     cex = 1, interpolate = FALSE, force_new = TRUE)
        dev.off()
      },
      
      # 6. rasterPoints Density Mode (Log-scale without interpolation)
      Rasterpoints_Density = {
        png(tempfile(), width = 800, height = 800)
        rasterPoints(data, colorder = "density", 
                     cex = 1, interpolate = FALSE, force_new = TRUE)
        dev.off()
      },
      
      iterations = 3,
      check = FALSE,
      memory = TRUE
    )
    
    res$n_points <- n
    results_list[[as.character(n)]] <- res
  }
  
  return(do.call(rbind, results_list))
}

# ------------------------------------------------------------------------------
# 4. OPENMP SCALABILITY BENCHMARK (MULTITHREADING)
# ------------------------------------------------------------------------------
run_openmp_scaling_benchmark <- function(n_points = 1e7, threads_vector = c(1, 2, 4, 8, 12, 16)) {
  max_avail <- parallel::detectCores()
  threads_vector <- threads_vector[threads_vector <= max_avail]
  
  message(sprintf("Evaluating OpenMP scalability with N = %s events...", format(n_points, scientific = FALSE)))
  fcs <- generate_fcs_data(n_points)
  data <- fcs$coords
  colors <- fcs$colors
  
  results <- bench::press(
    threads = threads_vector,
    {
      bench::mark(
        Kernel_Default = {
          png(tempfile(), width = 1000, height = 1000)
          rasterPoints(data, col = colors, colorder = "default", 
                       cex = 1, ncores = threads, force_new = TRUE)
          dev.off()
        },
        Kernel_Colorder = {
          png(tempfile(), width = 1000, height = 1000)
          rasterPoints(data, col = colors, colorder = c("gray70", "red3"), 
                       cex = 1, ncores = threads, force_new = TRUE)
          dev.off()
        },
        Kernel_Density = {
          png(tempfile(), width = 1000, height = 1000)
          rasterPoints(data, colorder = "density", 
                       ncores = threads, force_new = TRUE)
          dev.off()
        },
        iterations = 5,
        check = FALSE
      )
    }
  )
  
  df_res <- as.data.frame(results)
  df_res$threads <- results$threads
  df_res$expression <- as.character(results$expression)
  df_res$median_sec <- as.numeric(results$median)
  
  base_times <- aggregate(median_sec ~ expression, data = df_res[df_res$threads == 1, ], FUN = mean)
  
  df_res$speedup <- unlist(lapply(1:nrow(df_res), function(i) {
    expr <- df_res$expression[i]
    t1 <- base_times$median_sec[base_times$expression == expr]
    return(t1 / df_res$median_sec[i])
  }))
  
  return(df_res)
}

# ------------------------------------------------------------------------------
# 5. BENCHMARK EXECUTION
# ------------------------------------------------------------------------------
get_system_specs()

# Execute comparative rendering benchmark
bench_rendering <- run_rendering_benchmark(sample_sizes = c(1e5, 1e6, 5e6, 1e7))

# Execute parallel scalability benchmark (OpenMP)
bench_openmp <- run_openmp_scaling_benchmark(n_points = 1e7, threads_vector = c(1, 2, 4, 8, 16))

# Print console results
print("=== RENDERING AND MEMORY RESULTS ===")
print(bench_rendering, n = 24)

print("=== OPENMP SCALABILITY RESULTS ===")
clean_results <- bench_openmp
clean_results[, 3] <- as.numeric(clean_results[, 3])
clean_results[, 4] <- as.numeric(clean_results[, 4])
clean_results[, 6] <- as.numeric(clean_results[, 6])
clean_results[, 10] <- as.numeric(clean_results[, 10])
clean_results[, 15] <- as.numeric(clean_results[, 15])
print(clean_results[, -(12:14)])

# ------------------------------------------------------------------------------
# 6. MANUSCRIPT FIGURE GENERATION AND SAVING
# ------------------------------------------------------------------------------
generate_article_figures <- function(openmp_df) {
  message("Generating high-resolution figures for publication...")
  
  # Figure 1: Demonstration of rare population preservation using colorder
  fcs <- generate_fcs_data(1e6)
  png("Figure1_colorder_comparison.png", width = 1600, height = 800, res = 150)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  plot(fcs$coords, type = "n", xlab = "FSC-A", ylab = "SSC-A", 
       main = "A) Standard Rasterization (Random Overlap)")
  rasterPoints(fcs$coords, col = fcs$colors, colorder = "default", cex = 2, force_new = FALSE)
  
  plot(fcs$coords, type = "n", xlab = "FSC-A", ylab = "SSC-A", 
       main = "B) Priority-based Rendering (colorder = c('gray70', 'red3'))")
  rasterPoints(fcs$coords, col = fcs$colors, colorder = c("gray70", "red3"), cex = 2, force_new = FALSE)
  dev.off()
  
  # Figure 2: OpenMP Scalability Plot (Speedup)
  if (!missing(openmp_df)) {
    p_scaling <- ggplot(openmp_df, aes(x = threads, y = speedup, color = expression, group = expression)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
      scale_x_continuous(breaks = unique(openmp_df$threads)) +
      labs(
        title = "OpenMP Multi-threading Scalability (10 Million Events)",
        subtitle = "Dashed line represents theoretical ideal speedup",
        x = "Number of OpenMP Threads / Cores",
        y = "Speedup Factor (relative to 1 thread)",
        color = "C++ Kernel"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
    
    ggsave("Figure2_openmp_scaling.png", plot = p_scaling, width = 8, height = 6, dpi = 300)
  }
  
  message("Figures successfully generated!")
}

generate_article_figures(bench_openmp)
