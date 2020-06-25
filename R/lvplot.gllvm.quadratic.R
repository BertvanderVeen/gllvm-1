#' @title Plot covariate coefficients and confidence intervals
#' @description Plots covariate coefficients and their confidence intervals.
#'
#' @param object an object of class 'gllvm'.
#' @param plot.optima logical, defaults to \code{TRUE}. If \code{FALSE}, only plots species tolerances.
#' @param y.label logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param which.lvs vector indicating which LVs will be plotted. Default is \code{NULL} when optima and tolerances of all LVs are plotted.
#' @param cex.ylab the magnification to be used for axis annotation relative to the current setting of cex.
#' @param mfrow same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param xlim.list list of vectors with length of two to define the intervals for an x axis in each covariate plot. Defaults to NULL when the interval is defined by the range of point estimates and confidence intervals
#' @param level the confidence level. Scalar between 0 and 1.
#' @param ... additional graphical arguments.
#'
#' @details
#' Species optima and tolerances per latent variable are plotted with confidence intervals. Ablines at -2,2 and -5,5 are added for assistance, but do not have any meaning. Tolerances of which the CI of the quadratic coefficient crosses zero are greyed out. Tolerances are shown opposite of optima for readability, though are always positive. Optima of which the SE is larger than the optima are greyed out. Statistical uncertainties of optima that are larger than  15 or smaller than -15 are not plotted by default. Plots are ordered by variation explained (only from the quadratic term if plot.optima is FALSE).
#'
#' @author Bert van der Veen
#' @aliases lvplot lvplot.gllvm.quadratic
#' @export
lvplot.gllvm.quadratic <- function(object, plot.optima = TRUE, y.label = TRUE, which.lvs = NULL, cex.ylab = 0.5, mfrow = NULL, mar = c(4, 6, 2, 1),
                                   xlim.list = ifelse(plot.optima==T,rep(list(c(-5, 5)), length(which.lvs)),rep(list(c(0, 2)), length(which.lvs))), level = 0.95, ...) {
  if (any(class(object) != "gllvm.quadratic")) {
    stop("Class of the object isn't 'gllvm'.\n")
  }

  if (is.null(object$sd)) {
    warning("No standard errors present in model.\n")
  }
  if (is.null(which.lvs)) {
    which.lvs <- c(1:NCOL(object$lvs))
  }else{
    which.lvs <- sort(which.lvs)#for get resid cov below
  }

  cnames <- paste("LV", which.lvs)
  optima <- summary(object)$Optima[, which.lvs, drop = F]
  labely <- rownames(optima)
  p <- ncol(object$y)

  if (is.null(mfrow) && length(which.lvs) > 1) {
    mfrow <- c(1, length(which.lvs))
  }

  if (!is.null(mfrow)) {
    par(mfrow = mfrow, mar = mar)
  }
  if (is.null(mfrow)) {
    par(mar = mar)
  }
  LVidx <- 1:length(which.lvs)
  if(plot.optima==T){
    LVidx <- LVidx[order(getResidualCov(object)$trace.q[which.lvs]+getResidualCov(object)$trace.q2[which.lvs],decreasing=T)]
  }else{
    LVidx <- LVidx[order(getResidualCov(object)$trace.q2[which.lvs],decreasing=T)]
  }
  for (i in LVidx) {
    Xc <- optima[, i]
    sdoptima <- object$sd$optima[, which.lvs[i]]
    lower <- Xc + qnorm(level) * sdoptima
    upper <- Xc + qnorm(1 - level) * sdoptima
    Xc <- sort(Xc)
    sdoptima <- sdoptima[names(Xc)]
    lower <- lower[names(Xc)]
    upper <- upper[names(Xc)]

    col.seq <- rep("black", p)
    # grey out species with SD larger than optima
    col.seq[!sdoptima < abs(Xc)] <- "grey"

    # don't want to greyout optima.
    # col.seq[lower < 2*Xc & upper > 2*Xc] <- "grey"
    if (length(xlim.list) != length(which.lvs) & plot.optima == TRUE) {
      if (length(xlim.list) < which.lvs[i]) {
        xlim.list <- append(xlim.list, list(c(min(lower), max(upper))))
      }
    }

    At.y <- seq(1, p)

    if (plot.optima == TRUE) {
      plot(
        x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "o", cex.lab = 1.3,
        ...
      )
    }

    if (plot.optima == T) {
      for (j in 1:ncol(object$y)) {
        if (Xc[j] > (-15) & Xc[j] < 15) {
          segments(x0 = lower[j], y0 = At.y[j], x1 = upper[j], y1 = At.y[j], col = col.seq[j])
        }
      }
    }

    # tolerances
    tolerances <- 1 / sqrt(-2 * object$params$theta[, -c(1:object$num.lv)][, which.lvs[i]])
    sdtolerances <- object$sd$tolerances[, which.lvs[i]]
    lower <- tolerances + qnorm(level) * sdtolerances
    upper <- tolerances + qnorm(1 - level) * sdtolerances
    if (plot.optima == TRUE) tolerances <- tolerances[names(Xc)]
    if (plot.optima == FALSE) tolerances <- sort(tolerances)
    lower <- lower[names(tolerances)]
    upper <- upper[names(tolerances)]
    if (plot.optima == TRUE) sgn.opt <- -sign(Xc)
    if (plot.optima == FALSE) sgn.opt <- 1
    lower <- lower * sgn.opt
    upper <- upper * sgn.opt
    tolerances <- tolerances * sgn.opt

    col.seq <- rep("black", p)
    # grey out tolerances as if they cross 0 it's unclear if we have a quadratic response.

    CIquad <- confint(object)[-c(1:(object$num.lv * ncol(object$y))), ][((which.lvs[i] - 1) * p + 1):(which.lvs[i] * p), ]
    row.names(CIquad) <- colnames(object$y)
    CIquad <- CIquad[names(tolerances), ]
    # grey out species that are not sure to have a quadratic response
    col.seq[CIquad[, 1] < 0 & CIquad[, 2] > 0] <- "grey"

    if (length(xlim.list) != length(which.lvs) & plot.optima == FALSE) {
      if (length(xlim.list) < which.lvs[i]) {
        xlim.list <- append(xlim.list, list(c(min(lower), max(upper))))
      }
    }

    if (plot.optima == FALSE) {
      plot(
        x = tolerances, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "t", cex.lab = 1.3,
        ...
      )
    }

    if (plot.optima == TRUE) points(x = tolerances, y = At.y, pch = "t", col = col.seq)
    for (j in 1:ncol(object$y)) {
      if (tolerances[j] > (-10) & tolerances[j] < 10) {
        segments(x0 = lower[j], y0 = At.y[j], x1 = upper[j], y1 = At.y[j], col = col.seq[j])
      }
    }
    if (xlim.list[[i]][1] < (-5)) {
      abline(v = -5, lty = "dashed", col = "grey")
    }
    if (xlim.list[[i]][2] > 5) {
      abline(v = 5, lty = "dashed", col = "grey")
    }
    if (xlim.list[[i]][1] < (-2)) {
      abline(v = -2, lty = "dashed", col = "grey")
    }
    if (xlim.list[[i]][2] > 2) {
      abline(v = 2, lty = "dashed", col = "grey")
    }
    if (xlim.list[[i]][1] >= -2 & xlim.list[[i]][[1]] != 0 | xlim.list[[i]][2] <= 2 & xlim.list[[i]][[1]] != 0) {
      abline(v = 0, lty = "dashed", col = "grey")
    }

    if (y.label) {
      if (plot.optima == TRUE) axis(2, at = At.y, labels = names(Xc), las = 1, cex.axis = cex.ylab)
      if (plot.optima == FALSE) axis(2, at = At.y, labels = names(tolerances), las = 1, cex.axis = cex.ylab)
    }
  }
}

#' @export
lvplot <- function(object, ...) {
  UseMethod(generic = "lvplot")
}
