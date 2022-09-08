#' Color palette `RdBl`
#'
#' This function generates a color palette consisting of red (low), white (mid), and
#' blue (high). 
#'
#' @param n Number of gradations between red and white and between white and blue.
#' @return Sequence of colors with red and blue at the ends and white at the middle.
#' The length is `2*n+1`.
#' @export
RdBl.colors = function(n=128){
    c(hsv(1, seq(1, 0,length.out=n+1), 1),
      hsv(2/3, seq(1/n, 1, length.out=n), 1))
}

#' Color palette `BlWh`
#'
#' This function generates a color palette consisting of blue (low) and white (high).
#'
#' @param n Number of gradations between blue and white.
#' @return Sequence of colors with blue and white at the ends.
#' The length is `n`.
#' @export
BlWh.colors = function(n=256){
    hsv(2/3, seq(1, 0, length.out=n), 1)
}

#' Color palette `WhBl`
#'
#' This function generates a color palette consisting of white (low) and blue (high).
#'
#' @param n Number of gradations between white and blue.
#' @return Sequence of colors with blue and white at the ends.
#' The length is `n`.
#' @export
WhBl.colors = function(n=256){
    hsv(2/3, seq(0, 1, length.out=n), 1)
}

#' Image plot of matrix data
#'
#' Wrapper function of `heatmap` having a similar interface with `plt.imshow` in Python.
#'
#' @param x Matrix data input
#' @param vmin Minimum value of the colorscale. If not provided, the default is `-vmax`.
#' @param vmax Maximum value of the colorscale. If not provided, the default is the maximum absolute value of `x`.
#' @param col Color palette at which the image is plotted.
#' @return Imageplot of the input in the matrix position.
#' @export
imshow = function(x, vmin=NULL, vmax=NULL, col=RdBl.colors(256), ...){
    if(is.null(vmax)){
        vmax = max(abs(x))
    }
    if(is.null(vmin)){
        vmin = -vmax
    }
    
    heatmap(x, Rowv=NA, Colv=NA, revC=TRUE, scale="none",
            col = col, zlim = c(vmin,vmax), ...)
}