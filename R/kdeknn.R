#' @title kdewithknn
#' @name kdek
#' @description plot kde function with knn estimation
#' @param x data
#' @param k number of knn estimate
#' @param xrange range of x
#' @param yrange range of y
#' @param phi angles defining the viewing direction. phi give the colatitude.
#' @param theta angles defining the viewing direction. theta gives the azimuthal direction.
#' @param col the color(s) of the surface facets. Transparent colours are ignored. This is recycled to the (nx-1)(ny-1) facets.
#' @param border the color of the line drawn around the surface facets.
#' @return picture
#' @import FNN
#' @export
plotknnde=function(x,k,xrange,yrange,phi,theta,col,border){
  p=ncol(x)
  n=nrow(x)
  est_pt=expand.grid(xrange,yrange)
  distance=knnx.dist(x,est_pt,k)
  est_de=matrix(k/(2*n*distance[,k]),nrow = length(xrange))
  return(persp(xrange,yrange,est_de,phi,theta,col,border))
}