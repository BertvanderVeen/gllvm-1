#' @title Plot covariate coefficients and confidence intervals
#' @description Plots covariate coefficients and their confidence intervals.
#'
#' @param object an object of class 'gllvm'.
#' @param y.label logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param which.lvs vector indicating which LVs will be plotted. Default is \code{NULL} when optima and tolerances of all LVs are plotted.
#' @param cex.ylab the magnification to be used for axis annotation relative to the current setting of cex.
#' @param mfrow same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param xlim.list list of vectors with length of two to define the intervals for an x axis in each covariate plot. Defaults to NULL when the interval is defined by the range of point estimates and confidence intervals
#' @param ... additional graphical arguments.
#'
#' @author Bert van der Veen
#'@aliases lvplot lvplot.gllvm.quadratic
#'@export
lvplot.gllvm.quadratic <- function(object, y.label = TRUE, which.lvs = NULL, cex.ylab = 0.5, mfrow = NULL, mar = c(4, 6, 2, 1), 
                                     xlim.list = rep(list(c(-10,10)),length(which.lvs)), ...) {
  
  if (any(class(object) != "gllvm.quadratic")) 
    stop("Class of the object isn't 'gllvm'.")

    if (is.null(which.lvs)) 
    which.lvs <- c(1:NCOL(object$lvs))
    optima <- as.matrix(summary(object)$Optima[, which.lvs])
    cnames <- paste("LV",which.lvs)
    labely <- rownames(optima)
    m <- length(labely)
    Xc <- optima
    if (is.null(mfrow) && length(which.lvs) > 1) 
      mfrow <- c(1, length(which.lvs))
    
    if (!is.null(mfrow)) 
      par(mfrow = mfrow, mar = mar)
    if (is.null(mfrow)) 
      par(mar = mar)
    for (i in which.lvs) {
      Xc <- optima[, i]
      sdoptima <- object$sd$optima[, i]
      lower <- Xc - 1.96 * sdoptima
      upper <- Xc + 1.96 * sdoptima
      Xc <- sort(Xc)
      lower <- lower[names(Xc)]
      upper <- upper[names(Xc)]
      
      col.seq <- rep("black", m)
      #grey out optima whose CI is 2*optima, we're too uncertain about them
      col.seq[lower < 2*Xc & upper > 2*Xc] <- "grey"
      
      #don't want to greyout optima. 
      # col.seq[lower < 2*Xc & upper > 2*Xc] <- "grey"
      
      At.y <- seq(1, m)
      if (!is.null(xlim.list[[i]])) {
        plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "o", cex.lab = 1.3, 
             ...)
      } else {
        plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "o", 
             cex.lab = 1.3, ...)
      }
      for(j in 1:ncol(object$y)){
      if(Xc[j]>(-10)&Xc[j]<10)
      segments(x0 = nlower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
      }
      
      #tolerances
      tolerances <- 1/sqrt(-2*object$params$theta[,-c(1:object$num.lv)][,i])
      sdtolerances <- object$sd$tolerances[, i]
      lower <- tolerances - 1.96 * sdtolerances
      upper <- tolerances + 1.96 * sdtolerances
      tolerances <- sort(tolerances)
      lower <- lower[names(tolerances)]
      upper <- upper[names(tolerances)]
      
      col.seq <- rep("black", m)
      #grey out tolerances as if they cross 0 it's unclear if we have a quadartic response.
      col.seq[lower < 0 & upper > 0] <- "grey"
      points(x=tolerances, y = At.y - 0.2, pch="t", col = col.seq)
      segments(x0 = lower, y0 = At.y - 0.2, x1 = upper, y1 = At.y -0.2, col = col.seq)
      
            abline(v = 0, lty = "dashed")
      if (y.label) 
        axis(2, at = At.y, labels = names(Xc), las = 1, cex.axis = cex.ylab)
    }
}

#'@export
lvplot <- function(object, ...) {
  UseMethod(generic = "lvplot")
}

