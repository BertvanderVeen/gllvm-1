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
#'@details
#'Optima of which the SE is larger than the optima are greyed out. Tolerances of which the CI of the quadratic coefficient crosses zero are greyed out. Statistical uncertainties of optima that are larger than  15 or smaller than -15 are not plotted by default.
#'
#' @author Bert van der Veen
#'@aliases lvplot lvplot.gllvm.quadratic
#'@export
lvplot.gllvm.quadratic <- function(object, y.label = TRUE, which.lvs = NULL, cex.ylab = 0.5, mfrow = NULL, mar = c(4, 6, 2, 1), 
                                     xlim.list = rep(list(c(-15,15)),length(which.lvs)), ...) {
  
  if (any(class(object) != "gllvm.quadratic")) 
    stop("Class of the object isn't 'gllvm'.")

    if (is.null(which.lvs)) 
    which.lvs <- c(1:NCOL(object$lvs))

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
    optima <- as.matrix(summary(object)$Optima[, which.lvs])
    for (i in which.lvs) {
      Xc <- optima[, i]
      sdoptima <- object$sd$optima[, i]
      lower <- Xc - 1.96 * sdoptima
      upper <- Xc + 1.96 * sdoptima
      Xc <- sort(Xc)
      sdoptima <- sdoptima[names(Xc)]
      lower <- lower[names(Xc)]
      upper <- upper[names(Xc)]
      
      col.seq <- rep("black", m)
      #grey out species with SD larger than optima
      col.seq[!sdoptima<abs(Xc)] <- "grey"
      
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
      if(Xc[j]>(-15)&Xc[j]<15)
      segments(x0 = lower[j], y0 = At.y[j], x1 = upper[j], y1 = At.y[j], col = col.seq[j])
      }
      
      #tolerances
      tolerances <- 1/sqrt(-2*object$params$theta[,-c(1:object$num.lv)][,i])
      sdtolerances <- object$sd$tolerances[, i]
      lower <- tolerances - 1.96 * sdtolerances
      upper <- tolerances + 1.96 * sdtolerances
      tolerances <- tolerances[names(Xc)]
      lower <- lower[names(tolerances)]
      upper <- upper[names(tolerances)]
      
      col.seq <- rep("black", m)
      #grey out tolerances as if they cross 0 it's unclear if we have a quadartic response.
      
      CIquad <- confint(object)[-c(1:object$num.lv)*ncol(object$y),][1:(object$num.lv*ncol(object$y)),][(1+(i-1)*p):(i*p),]
      #grey out species that are not sure to have a quadratic response
      col.seq[CIquad[,1]<0& CIquad[,2]>0] <- "grey"
      points(x=tolerances, y = At.y - 0.2, pch="t", col = col.seq)
      for(j in 1:ncol(object$y)){
        if(tolerances[j]>(-10)&tolerances[j]<10)
          segments(x0 = lower[j], y0 = At.y[j]-0.2, x1 = upper[j], y1 = At.y[j]-0.2, col = col.seq[j])
      }
      
            abline(v = 0, lty = "dashed")
      if (y.label) 
        axis(2, at = At.y, labels = names(Xc), las = 1, cex.axis = cex.ylab)
    }
}

#'@export
lvplot <- function(object, ...) {
  UseMethod(generic = "lvplot")
}

