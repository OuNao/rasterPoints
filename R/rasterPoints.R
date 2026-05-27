#' @title rasterPoints.
#'
#' @description Plot scatter points using rasterImage. Open a graphic device if grDevices::dev.cur() == 1 (null device).
#' @param x Matrix of data to plot. Must have 2 columns.
#' @param col Character vector with lenght = nrow(x) or lenght = 1 (same color for all points).
#' @param cex Integer vector. See par("cex").
#' @param interpolate Logical. Passed to rasterImage.
#' @param ncores Integer. Number of cores used (OpenMP). Defaul = 0 (max avaiable cores).
#' @param colorder Character vector with the color priority order. Colors with higher colorder index has priority over the others. The "default" behavior is set the point color only the first time and ignore any new attempt to set the point color.
#' @param force_new Logical. Force create a new plot.
#' @param ... Graphic parameters passed to rasterImage/plot (see par()).
#' @return Nothing.
#' @keywords raster, scatter, plot, points.
#' @examples 
#' data<-matrix(rnorm(10^6), ncol = 2)
#' system.time({
#'   plot(data, type="n")
#'   points(data, pch = ".", type = "p", cex=1, col = "blue")
#' })
#' system.time({
#'   plot(data, type="n")
#'   rasterPoints(data, cex=1, col = "blue")
#' })
#' @export
rasterPoints<-function(x, col="black", cex=1, interpolate = F, ncores = 0, colorder = "default", force_new = FALSE, colramp = NULL, ...) {
  if (!is.matrix(x) || ncol(x)!=2 || nrow(x)<1) stop("x must be a matrix with 2 columns and >0 rows!", call. = F)
  if (!is.vector(col) || !is.character(col) || (length(col)!=1 && length(col) != nrow(x))) stop("col mus be a character vector of lenght 1 or nrow(x)", call. = F)
  if (!is.vector(cex) || !is.numeric(cex) || any(cex!=as.integer(cex)) || length(cex)==0) stop("cex must be a integer vector of lenght >0.", call. = F)
  if (!is.vector(interpolate) || !is.logical(interpolate) || length(interpolate)!=1) stop("interpolate must be a logical (TRUE or FALSE)", call. = F)
  if (!is.vector(ncores) || !is.numeric(ncores) || ncores!=as.integer(ncores) || length(ncores)!=1) stop("ncores must be a integer (integer vector of lenght = 1).", call. = F)
  if (length(col)==1) col=rep(do.call(grDevices::rgb, as.list(grDevices::col2rgb(col)/255)), nrow(x))
  if (length(cex)==1) cex=rep(cex, nrow(x))
  if (force_new) graphics::plot(x, type = "n", ...)
  tryCatch(graphics::par(new=TRUE),error=function(e) e, warning=function(w) graphics::plot(x, type = "n", ...))
  col_factor <- factor(col)
  col_idx <- as.integer(col_factor)
  unique_colors <- levels(col_factor)
  if (colorder[1] == "default") {
    colorder_vec <- rep(0L, length(unique_colors))
  } else if (colorder[1] == "density") {
    if (is.null(colramp) || !is.function(colramp)) {
      colramp<-grDevices::colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
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
  if (colorder[1] == "density") {
    img_idx <- data2raster_density(x, usr, size[1], size[2], 256, ncores)
  } else {
    img_idx <- data2raster(x, col_idx, colorder_vec, usr, size[1], size[2], cex, ncores)
  }
  img_idx[img_idx == 0] <- NA
  img <- matrix(unique_colors[img_idx], nrow = nrow(img_idx))
  
  rast<-grDevices::as.raster(img)
  graphics::rasterImage(rast, 
                        xleft=usr[1],
                        xright=usr[2],
                        ybottom=usr[3],
                        ytop=usr[4], interpolate = interpolate, ...)
}