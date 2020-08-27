#' @title Summarizing gllvm model fits
#' @description A summary of the fitted 'gllvm' object, including function call, distribution family and model parameters.
#'
#' @param object   an object of class 'gllvm'
#' @param ... not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @details
#' Also gives species optimum, tolerances and maxima per latent variable. Maxima are on the link scale.
#'
#' @examples
#' ## Load a dataset from the mvabund package
#' data(antTraits)
#' y <- as.matrix(antTraits$abund)
#' # Fit gllvm model
#' fit <- gllvm(y = y, family = poisson())
#' summary(fit)
#' @export

summary.gllvm.quadratic <- function(object, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- dim(object$X)[2]
  nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  family <- object$family

  M <- cbind(object$params$beta0, object$params$theta)
  opt <- -object$params$theta[, 1:object$num.lv, drop = F] / (2 * object$params$theta[, -c(1:object$num.lv), drop = F])
  colnames(opt) <- paste("LV", 1:object$num.lv, sep = "")
  if (is.null(row.names(opt))) row.names(opt) <- names(object$params$beta0)
  tol <- 1 / sqrt(-2 * object$params$theta[, -c(1:object$num.lv), drop = F])

  max <- object$params$beta0 + opt * object$params$theta[, 1:object$num.lv, drop = F] + opt^2 * object$params$theta[, -c(1:object$num.lv), drop = F]

  row.names(max) <- row.names(tol) <- colnames(object$y)
  colnames(max) <- colnames(tol) <- colnames(opt)

  sumry <- list()
  sumry$"log-likelihood" <- object$logL
  crit <- inf.criteria(object)
  sumry$df <- crit$k
  sumry$AIC <- crit$AIC
  sumry$AICc <- crit$AICc
  sumry$BIC <- crit$BIC

  newnams <- c("Intercept", c(paste("theta.LV", 1:num.lv, sep = ""), paste("theta.LV^2", 1:num.lv, sep = "")))
  colnames(M) <- newnams
  rownames(M) <- colnames(object$y)
  sumry$Call <- object$call
  sumry$family <- object$family
  sumry$Coefficients <- M
  sumry$Optima <- opt
  sumry$Tolerances <- tol
  sumry$Maxima <- max
  if(object$common.tolerances==FALSE){
  sumry$Gradient.length <- 4*sqrt(0.5)*1/apply(tol,2,median)
  }else{
  sumry$Gradient.length <- 4*sqrt(0.5)*1/tol[1,]
  }
  names(sumry$Gradient.length) <- paste("LV",1:object$num.lv,sep="")

  if (!is.null(object$TR)) {
    if (!is.null(object$X)) {
      sumry$"Covariate coefficients" <- object$params$B
    }
  } else {
    if (!is.null(object$X)) {
      sumry$"Environmental coefficients" <- object$params$Xcoef
    }
  }
  if (!is.null(object$params$row.params)) {
    sumry$"Row intercepts" <- object$params$row.params
  }

  if (object$row.eff == "random") {
    object$params$sigma2 <- object$params$sigma^2
    names(object$params$sigma2) <- "sigma^2"
    sumry$"Variance of random row intercepts" <- object$params$sigma2
  }

  if (object$family == "negative.binomial") {
    sumry$"Dispersion parameters" <- object$params$phi
  }
  class(sumry) <- "summary.gllvm.quadratic"
  return(sumry)
}
