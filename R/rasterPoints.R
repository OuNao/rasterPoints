#' @title rasterPoints.
#'
#' @description Plot scatter points using rasterImage. Open a graphic device if grDevices::dev.cur() == 1 (null device).
#' @param x Matrix of data to plot. Must have 2 columns.
#' @param col Character vector with lenght = nrow(x) or lenght = 1 (same color for all points).
#' @param cex Integer vector. See par("cex").
#' @param interpolate Logical. Passed to rasterImage.
#' @param ncores Integer. Number of cores used (OpenMP). Defaul = 0 (max avaiable cores).
#' @param colorder Character vector with the color priority order. Colors with higher colorder index has priority over the others. The "default" behavior is set the point color only the first time and ignore any new attempt to set the point color.
#' @param smooth Logical. Apply Gaussian neighborhood kernel smoothing in density mode?
#' @param smooth_radius Kernel neighborhood radius in pixels. Default is 4L.
#' @param smooth_sigma Gaussian standard deviation. Default is 2.0.
#' @param margin_pct Interior bounding box margin ratio (0.05). Prevents axis saturation artifacts.
#' @param force_new Logical. Force create a new plot.
#' @param ... Graphic parameters passed to rasterImage/plot (see par()).
#' @return Nothing.
#' @keywords raster, scatter, plot, points.
#' @examples 
#' data<-matrix(rnorm(1e6), ncol = 2)
#' system.time({
#'   plot(data, type="n")
#'   points(data, pch = ".", type = "p", cex=1, col = "blue")
#' })
#' system.time({
#'   plot(data, type="n")
#'   rasterPoints(data, cex=1, col = "blue")
#' })
#' @export
rasterPoints<-function(x, col="black", cex=1, 
                       interpolate = F, ncores = 0, 
                       colorder = "default", smooth = TRUE, 
                       smooth_radius = 4L, smooth_sigma = 2.0,
                       margin_pct = 0.05, force_new = FALSE, 
                       colramp = NULL, ...) {
  if (!is.matrix(x) || ncol(x)!=2 || nrow(x)<1) stop("x must be a matrix with 2 columns and >0 rows!", call. = F)
  if (any(!is.finite(x))) stop("The coordinate matrix 'x' cannot contain non-finite values (NA, NaN, Inf, -Inf)!", call. = F)
  if (!is.vector(col) || !is.character(col) || (length(col)!=1 && length(col) != nrow(x))) stop("col mus be a character vector of lenght 1 or nrow(x)", call. = F)
  if (any(is.na(col))) stop("col cannot contain NA values!", call. = F)
  if (!is.vector(cex) || !is.numeric(cex) || any(cex!=as.integer(cex)) || length(cex)==0) stop("cex must be a integer vector of lenght >0.", call. = F)
  if (!is.vector(interpolate) || !is.logical(interpolate) || length(interpolate)!=1) stop("interpolate must be a logical (TRUE or FALSE)", call. = F)
  if (!is.vector(ncores) || !is.numeric(ncores) || ncores!=as.integer(ncores) || length(ncores)!=1) stop("ncores must be a integer (integer vector of lenght = 1).", call. = F)
  if (length(col)==1) {
    col_rgba <- grDevices::col2rgb(col, alpha = TRUE)
    col <- rep(grDevices::rgb(col_rgba[1,]/255, col_rgba[2,]/255, col_rgba[3,]/255, col_rgba[4,]/255), nrow(x))
  }
  if (length(cex)==1) cex=rep(cex, nrow(x))
  if (force_new) graphics::plot(x, type = "n", ...)
  tryCatch(graphics::par(new=TRUE),error=function(e) e, warning=function(w) graphics::plot(x, type = "n", ...))
  col_factor <- factor(col)
  col_idx <- as.integer(col_factor)
  unique_colors <- levels(col_factor)
  if (colorder[1] != "default" && colorder[1] != "density"){
    if (any(is.na(colorder))) stop("colorder cannot contain NA values!", call. = F)
    if (!all(colorder %in% unique_colors)) warning(" Some colors in 'colorder' do not exist in the 'col' vector!")
  }
  if (colorder[1] == "default") {
    colorder_vec <- rep(0L, length(unique_colors))
  } else if (colorder[1] == "density") {
    if (is.null(colramp) || !is.function(colramp)) {
      colramp<-grDevices::colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"), alpha = TRUE)
    }
    unique_colors <- colramp(256)
  } else {
    colorder_vec <- rep(0L, length(unique_colors))
    m <- match(colorder, unique_colors)
    valid_m <- !is.na(m)
    colorder_vec[m[valid_m]] <- seq_along(colorder)[valid_m]
  }
  usr <- graphics::par('usr')
  psize<-grDevices::dev.size('px')
  pict<-c(graphics::grconvertX(usr[1:2],"user", "ndc")*psize[1], graphics::grconvertY(usr[3:4], "user", "ndc") * psize[2])
  pict<-as.integer(pict)
  usr<-c(graphics::grconvertX(pict[1:2]/psize[1], "ndc","user"), graphics::grconvertY(pict[3:4]/psize[2], "ndc", "user"))
  size<-c(pict[2]-pict[1]+1, pict[4]-pict[3]+1)
  if (any(psize <= 0) || any(size <= 0)) stop("Invalid plot device size. Ensure graphics device is open and has a valid positive size.", call. = F)
  if (colorder == "density") {
    if (smooth) {
      if (smooth_sigma <= 0) smooth_sigma <- 2.0
      if (smooth_radius <= 0) smooth_radius <- 4L
    }
    img_idx <- data2raster_density(
      x = x, 
      usr = usr, 
      width = as.integer(size[1]), 
      height = as.integer(size[2]), 
      n_bins = 256L,
      smooth = as.logical(smooth),
      smooth_radius = as.integer(smooth_radius),
      smooth_sigma = as.numeric(smooth_sigma),
      margin_pct = as.numeric(margin_pct),
      cex = cex[1],
      ncores = as.integer(ncores)
    )
  } else {
    img_idx <- data2raster(x, col_idx, colorder_vec, usr, size[1], size[2], cex, ncores)
  }
  mapping_colors <- c("#00000000", unique_colors)
  img <- matrix(mapping_colors[img_idx + 1L], nrow = nrow(img_idx))
  
  rast<-grDevices::as.raster(img)
  graphics::rasterImage(rast, 
                        xleft=usr[1],
                        xright=usr[2],
                        ybottom=usr[3],
                        ytop=usr[4], interpolate = interpolate, ...)
}