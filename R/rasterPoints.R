#' @title A function to merge flow cytometry data files (.FCS).
#'
#' @description Plot points using rasterImage.
#' @param x Data to plot. 2 column matrix.
#' @param col Vector of common columns number
#' @param cex Vector of variable columns number
#' @param interpolate Quantile normalize the common parameters before imputation
#' @return Merged flowFrame
#' @keywords FCS
#' @export
rasterPoints<-function(x, col, cex=1, interpolate = F) {
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