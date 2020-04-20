#' @title rasterPoints.
#'
#' @description Plot scatter points using rasterImage.
#' @param x Matrix of data to plot. Must have 2 columns.
#' @param col Character vector with lenght = nrow(x) or lenght = 1 (same color for all points).
#' @param cex Integer. See par("cex").
#' @param interpolate Logical. Passed to rasterImage.
#' @return Nothing.
#' @keywords raster, scatter, plot, points.
#' @export
rasterPoints<-function(x, col, cex=1, interpolate = F) {
  if (!is.matrix(x) || ncol(x)!=2 || nrow(x)<1) stop("x must be a matrix with 2 columns and >0 rows!", call. = F)
  if (!is.vector(col) || !is.character(col) || (length(col)!=1 && length(col) != nrow(x))) stop("col mus be a character vector of lenght 1 or nrow(x)", call. = F)
  if (!is.vector(cex) || !is.numeric(cex) || cex!=as.integer(cex) || length(cex)!=1) stop("cex must be a integer (integer vector of lenght = 1).", call. = F)
  if (!is.vector(interpolate) || !is.logical(interpolate) || length(interpolate)!=1) stop("interpolate must be a logical (TRUE or FALSE)", call. = F)
  if (length(col)==1) col=do.call(rgb, as.list(col2rgb(col)/255))
  usr <- graphics::par('usr')
  psize<-grDevices::dev.size('px')
  pict<-c(graphics::grconvertX(usr[1:2],"user", "ndc")*psize[1], graphics::grconvertY(usr[3:4], "user", "ndc") * psize[2])
  pict<-as.integer(pict)
  usr<-c(graphics::grconvertX(pict[1:2]/psize[1], "ndc","user"), graphics::grconvertY(pict[3:4]/psize[2], "ndc", "user"))
  size<-c(pict[2]-pict[1], pict[4]-pict[3])
  img<-data2raster(x, col, usr, size[1], size[2], cex)
  rast<-as.raster(img)
  graphics::rasterImage(rast, 
                        xleft=usr[1],
                        xright=usr[2],
                        ybottom=usr[3],
                        ytop=usr[4], interpolate = interpolate)
}