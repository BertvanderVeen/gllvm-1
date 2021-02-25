#' @title Plot species optima for latent variables from a gllvm
#' @description Plots latent variables and their corresponding species optima.
#'
#' @param object   an object of class 'gllvm'.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot (sorted by varianve explained). The default is all. Can supply a vector, in which case those species will be picked from the order of variation explained on the latent variable (e.g. the species that explain fourth and sixth most variation).
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
#' @param s.colors colors for sites
#' @param s.labels Can be FALSE, for no labels, "row.names" or "rug".
#' @param cex.spp size of species labels in biplot
#' @param scale For 2D plots, to scale optima, tolerances, and site scores by the ratio variance explained. OLD: Alternatively can be "tolerances" to scale optima by tolerances and site scores by average tolerances per latent variable.
#' @param opt.region Only for 2D plots, efaults to FALSE. If "confidence", plots confidence intervals for species optima as ellipses. If "distribution" plots predicted species distributions.
#' @param type Can be used to predict on the response or link scale. Default is response (except for ordinal, for which the only option is "link").
#' @param intercept Can be used to include species-intercepts in the plot. Default is TRUE
#' @param legend when \code{TRUE} adds legend in the topleft corner of the plot, instead of species names in the plot
#' @param site.region logical, if \code{TRUE} prediction regions for the predicted latent variables are plotted, defaults to \code{FALSE}.
#' @param level level for prediction regions. Can be a vector of size two, then the first index is the level for the site prediction regions and the second for optima.
#' @param lty.ellips line type for prediction ellipses. Can be a vector of size two. See graphical parameter lty.
#' @param lwd.ellips line width for prediction ellipses. Can be a vector of size two. See graphical parameter lwd.
#' @param col.ellips colors for site prediction ellipses.
#' @param level level for prediction regions. Can be a vector of size two. Defaults to 0.95.
#' @param alpha.col used to control the transparency of confidence interval ribbons in 1D plot, and prediction regions in 2D plot.
#' @param length.out number of points to evaluate the function at for 1D plots. Defaults to 100.
#' @param xlim
#' @param ylim
#' @param ylab
#' @param ... additional graphical arguments.
#'
#' @details
#' Function constructs a plot of species optima. in 1D a dashed line indicates the optimum of species. In 2D species distributions can be added as ellipses.
#'
#' @author Bert van der Veen
#'
#' @examples
#' #' ## Load a dataset from the mvabund package
#' data(antTraits)
#' y <- as.matrix(antTraits$abund)
#' fit <- gllvm.quadratic(y, family = poisson())
#' # Ordination plot:
#' optiplot(fit)
#' # Optiplot with 10 species
#' optiplot(fit, ind.spp = 10)
#' @aliases optiplot optiplot.gllvm.quadratic
#' @export
#' @export optiplot.gllvm.quadratic
optiplot.gllvm.quadratic <- function(object, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = NULL,
                                     s.colors = 1, s.labels = "rug", cex.spp = 0.7, opt.region = ifelse(length(which.lvs) == 1, "confidence", "distribution"), type = "response", intercept = TRUE, legend = FALSE, scale = FALSE, site.region = FALSE, lty.ellips = c("solid", "dashed"), lwd.ellips = 1, col.ellips = "gray", alpha.col = 0.4, level = 0.95, ylim = NULL, xlim = NULL, ylab = NULL, length.out = 100, ...) {
  if (class(object) != "gllvm.quadratic") {
    stop("Class of the object isn't 'gllvm.quadratic'. linear GLLVM not implemented yet.\n")
  }
  if (!is.list(object$sd)) {
    warning("No standard errors present in model, setting `opt.region = FALSE`.\n")
    opt.region <- FALSE
  }
  if (!opt.region %in% c("distribution", "confidence", FALSE)) {
    stop("Wrong input for `opt.region`.\n")
  }
  # if(object$family=="binomial"&opt.region=="confidence"&type=="response"){
  #   stop("Plots on the response scale with confidence intervals not yet supported for the binomial distribution. \n")
  # }
  if (!is.null(object$X)) {
    if (any(round(apply(object$X, 2, mean), 0) != 0) & length(which.lvs)==1) {
      warning("Your plot will look best with standardized covariates.\n")
    }
  }
  # if(object$row.eff!=FALSE&opt.region=="confidence"&length(which.lvs)==1){
  #   warning("Confidence interval bands not yet implement with row-effects.\n")
  #   opt.region <- FALSE
  # }
  n <- NROW(object$y)
  p <- NCOL(object$y)

  if (!is.null(ind.spp)) {
    if (length(ind.spp) == 1) {
      if (ind.spp > p) {
        stop("There are not that many species in the model. \n")
      } else {
        ind.spp <- 1:ind.spp
      }
    }
  } else {
    ind.spp <- 1:p
  }

  if (object$num.lv == 0) {
    stop("No latent variables to plot.")
  }
  if (is.null(which.lvs)) {
    if (object$num.lv > 1) {
      which.lvs <- 1:2
    } else {
      which.lvs <- 1
    }
  }
  if (length(which.lvs) > object$num.lv) stop("More latent variables select than included in the model.")
  if (is.null(rownames(object$params$theta))) {
    rownames(object$params$theta) <- names(object$params$beta0)
  }
  if (is.null(rownames(object$params$theta))) {
    rownames(object$params$theta) <- paste("v", 1:p)
  }

  if (length(which.lvs) == 1) {
    # if(!is.null(object$X)){
    #   intercept <- FALSE
    #   warning("Setting intercept to FALSE due to included covariates.")
    # }
    resid.cov <- object$params$theta[, which.lvs, drop = F]^2 + 2 * object$params$theta[, -c(1:object$num.lv), drop = F][, which.lvs, drop = F]^2
    largest.lnorms <- order(rowSums(resid.cov), decreasing = TRUE)[ind.spp]
    cols <- (grDevices::rainbow(ncol(object$y) + 1)[2:(ncol(object$y) + 1)])[largest.lnorms]
    if (object$num.lv == 1) {
      quadr.coef <- object$params$theta[largest.lnorms, -1, drop = F]
    } else {
      quadr.coef <- object$params$theta[, -c(1:object$num.lv), drop = F][largest.lnorms, which.lvs, drop = F]
    }

    quadr.coef[which(round(quadr.coef, 3) == 0)] <- 0

    if (type == "response") {
      if (object$family == "gaussian") {
        linkinv <- function(x) x
      } else if (object$family == "binomial") {
        linkinv <- pnorm
      } else if (object$family %in% c("poisson", "negative.binomial", "gamma")) {
        linkinv <- exp
      } else if (object$family == "ordinal") {
        warning("Ordinal model plot not yet implemented on response scale")
        type <- "link"
      }
    }

    # mu <- predict(object, newLV = newLV, LVonly = T, which.lvs = which.lvs,type = type,intercept=intercept)[,largest.lnorms,drop=F]
    # have to rebuiltin the backtransformation of the linear predictor now I included curve
    mu <- predict(object, LVonly = T, which.lvs = which.lvs, type =F, intercept = F)[, largest.lnorms, drop = F]
    
    opt <- object$params$theta[largest.lnorms,which.lvs,drop=F]/(2*abs(object$params$theta[,-c(1:object$num.lv),drop=F][largest.lnorms,which.lvs,drop=F]))
    c <- object$params$beta0[largest.lnorms] + rowSums(opt*object$params$theta[largest.lnorms,which.lvs,drop=F]) 
    tol <- 1/sqrt(2*abs(object$params$theta[,-c(1:object$num.lv),drop=F][largest.lnorms,which.lvs,drop=F]))
    if(intercept==TRUE){
    for(j in 1:ncol(mu)){
      mu[,j] <- mu[,j] + c[j] - sum(opt[j,]^2/(2*tol[j,]^2))
      }
      }
    mu <- linkinv(mu)

    if (is.null(ylim)) {
      ylim <- range(mu)
    }
    if(is.null(ylab)){
    ylab <- "Predicted"
    }
    if (legend == F) {
      # pdf(NULL)
      # plot(NA, xlim = c(range(object$lvs)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", ...)
      # maxstr<-max(strwidth(colnames(mu)))*1.25
      # invisible(dev.off())
      # previpously range below was range(object$lvs)+maxstr, might still want to add this again
      if (is.null(xlim)) {
        xlim <- c(c(min(object$lvs[, which.lvs]) - .05, max(object$lvs[, which.lvs]) + .05))
      }
      plot(NA, xlim = xlim, ylim = ylim, ylab = ylab, xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", main = main, ...)
    } else {
      if (is.null(xlim)) {
        xlim <- c(range(object$lvs[, which.lvs]))
      }
      plot(NA, xlim = xlim, ylim = ylim, ylab = ylab, xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", main = main, ...)
    }

    if (legend == T) {
      legend(x = "topleft", text.col = cols, colnames(mu), cex = cex.spp)
    }

    if (opt.region == "confidence") {
      V <- -object$Hess$cov.mat.mod
      colnames(V) <- colnames(object$Hess$Hess.full[object$Hess$incl, object$Hess$incl])
      row.names(V) <- row.names(object$Hess$Hess.full[object$Hess$incl, object$Hess$incl])
      V <- V[colnames(V) != "r0", colnames(V) != "r0"]
      # remove everything before the intercept
      # V<-solve(object$Hess$Hess.full[object$Hess$incl,object$Hess$incl][idx,idx])
      # stilll need to reduce V to dim(mu) here
      if (intercept == F) {
        for (q in 1:object$num.lv) {
          for (j in 1:ncol(object$y)) {
            if (q > 1 & j < q) {
              # add zeros where necessary
              V <- cbind(cbind(V[, 1:c(p * q + j - 1)], 0), V[, (p * q + j):ncol(V)])
              V <- rbind(rbind(V[1:c(p * q + j - 1), ], 0), V[(p * q + j):nrow(V), ])
            }
          }
        }
      } else {
        for (q in 1:object$num.lv) {
          for (j in 1:ncol(object$y)) {
            if (q > 1 & j < q) {
              # add zeros where necessary
              V <- cbind(cbind(V[, 1:c(p + p * (q - 1) + j - 1)], 0), V[, c(p + p * (q - 1) + j):ncol(V)])
              V <- rbind(rbind(V[1:c(p + p * (q - 1) + j - 1), ], 0), V[c(p + p * (q - 1) + j):nrow(V), ])
            }
          }
        }
      }
      colnames(V)[colnames(V) == ""] <- "lambda"
      row.names(V)[row.names(V) == ""] <- "lambda"
    }

    if (!is.null(xlim)) {
      minLV <- xlim[1]
      maxLV <- xlim[2]
    } else {
      minLV <- min(object$lvs[, which.lvs])
      maxLV <- max(object$lvs[, which.lvs])
    }
    if (intercept == F) {
      if (type == "link") {
        func <- function(x, u, u2) x * u + x^2 * u2
      } else {
        func <- function(x, u, u2) linkinv(x * u + x^2 * u2)
      }
    } else {
      if (type == "link") {
        func <- function(x, beta, u, u2) beta + x * u + x^2 * u2
      } else {
        func <- function(x, beta, u, u2) linkinv(beta + x * u + x^2 * u2)
      }
    }
    for (j in 1:ncol(mu)) {
      curvePlot <- list(x = NA, y = NA)
      curvePlot$x <- seq(minLV, maxLV, length.out = length.out)
      if (intercept == F) {
        curvePlot$y <- func(curvePlot$x, u = object$params$theta[largest.lnorms, , drop = F][j, which.lvs], u2 = object$params$theta[largest.lnorms, , drop = F][j, -(1:object$num.lv), drop = F][, which.lvs])
      } else {
        curvePlot$y <- func(curvePlot$x, beta = object$params$beta0[largest.lnorms][j], u = object$params$theta[largest.lnorms, , drop = F][j, which.lvs], u2 = object$params$theta[largest.lnorms, , drop = F][j, -(1:object$num.lv), drop = F][, which.lvs])
      }

      lines(x = curvePlot$x, y = curvePlot$y, col = cols[j])

      if (opt.region == "confidence") {
        if (intercept == F) {
          if (object$ridge[[2]] > 0 & object$common.tolerances == F) {
            X <- cbind(curvePlot$x, curvePlot$x^2, curvePlot$x^2)
          } else {
            X <- cbind(curvePlot$x, curvePlot$x^2)
          }
        } else {
          if (object$ridge[[2]] > 0 & object$common.tolerances == F) {
            X <- cbind(1, curvePlot$x, curvePlot$x^2, curvePlot$x^2)
          } else {
            X <- cbind(1, curvePlot$x, curvePlot$x^2)
          }
        }
        if (intercept == F) {
          if (object$ridge[[2]] > 0) {
            idx <- colnames(V) == "lambda" | colnames(V) == "lambda2" | colnames(V) == "lambda3"
            V.theta <- V[idx, idx]
            if (object$common.tolerances == T) {
              idx <- c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + ncol(V.theta) - object$num.lv):ncol(V.theta))[which.lvs])
            } else {
              idx <- c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + p * object$num.lv + (object$num.lv * (j - 1))):(p * object$num.lv + (object$num.lv * (j - 1)) + object$num.lv))[which.lvs], ((1 + ncol(V.theta) - object$num.lv):ncol(V.theta))[which.lvs])
            }
          } else {
            idx <- colnames(V) == "lambda" | colnames(V) == "lambda2"
            V.theta <- V[idx, idx]
            if (object$common.tolerances == T) {
              idx <- c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + ncol(V.theta) - object$num.lv):ncol(V.theta))[which.lvs])
            } else {
              idx <- c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + p * object$num.lv + (object$num.lv * (j - 1))):(p * object$num.lv + (object$num.lv * (j - 1)) + object$num.lv))[which.lvs])
            }
          }
          V.theta2 <- V.theta[idx, idx]
        } else {
          if (object$ridge[[2]] > 0) {
            idx <- colnames(V) == "b" | colnames(V) == "lambda" | colnames(V) == "lambda2"
            V.theta <- V[idx, idx]
            if (object$common.tolerances == T) {
              idx <- c(j, p + c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], c(p * object$num.lv + 1:object$num.lv)[which.lvs]))
            } else { # old index for quadratic was c(1+object$num.lv*p+(j-1*2)+j:(j+object$num.lv-1))[which.lvs]
              idx <- c(j, p + c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + p * object$num.lv + (object$num.lv * (j - 1))):(p * object$num.lv + (object$num.lv * (j - 1)) + object$num.lv))[which.lvs], ((1 + ncol(V.theta) - object$num.lv):ncol(V.theta))[which.lvs]))
            }

            V.theta2 <- V.theta[idx, idx]
          } else {
            idx <- colnames(V) == "b" | colnames(V) == "lambda" | colnames(V) == "lambda2"
            V.theta <- V[idx, idx]
            if (object$common.tolerances == T) {
              idx <- c(j, p + c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], c(p * object$num.lv + 1:object$num.lv)[which.lvs]))
            } else {
              idx <- c(j, p + c(c((c(1:object$num.lv) - 1) * p + j)[which.lvs], ((1 + p * object$num.lv + (object$num.lv * (j - 1))):(p * object$num.lv + (object$num.lv * (j - 1)) + object$num.lv))[which.lvs]))
            }

            V.theta2 <- V.theta[idx, idx]
          }
        }

        # need to do this per species and X needs to be the x evaluated at the grid
        se <- sqrt(abs(rowSums(X * (X %*% V.theta2))))

        curveLowCI <- curvePlot$y + se * qnorm(level)
        curveUpCI <- curvePlot$y + se * qnorm(1 - level)

        cols2 <- scales::alpha(cols[j], ifelse(length(alpha.col) > 2, alpha.col[2], alpha.col[1]))
        polygon(c(curvePlot$x, rev(curvePlot$x)), c(curveLowCI, rev(curveUpCI)), border = NA, col = cols2)
      }
      if (legend == F) {
        opt <- summary(object)$Optima[largest.lnorms, which.lvs]
        if (intercept == TRUE) {
          maximum <- summary(object)$Maxima[largest.lnorms, which.lvs][j]
          if (type == "response") maximum <- linkinv(maximum)
        } else {
          maximum <- (summary(object)$Optima * object$params$theta[, 1:object$num.lv, drop = F] + summary(object)$Optima^2 * object$params$theta[, -c(1:object$num.lv), drop = F])[largest.lnorms, which.lvs][j]
          if (type == "response") {
            maximum <- linkinv(maximum)
          }
        }

        if (type == "response") {
          if (opt[j] < max(object$lvs[, which.lvs]) & opt[j] > min(object$lvs[, which.lvs])) {
            text(x = opt[j], y = maximum, labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 3) # should adjust "adj" rather than adding 0.1 to the maximum.
            segments(x0 = opt[j], x1 = opt[j], y0 = 0, y1 = maximum, lty = "dashed", col = cols[j])
          } else if (opt[j] < min(object$lvs[, which.lvs])) {
            text(x = min(object$lvs[, which.lvs]), y = apply(mu, 2, max)[j], labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 4)
          } else if (opt[j] > max(object$lvs[, which.lvs])) {
            text(x = max(object$lvs[, which.lvs]), y = apply(mu, 2, max)[j], labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 2)
          }
        } else { # not going correct yet for poisson
          if (opt[j] < max(object$lvs[, which.lvs]) & opt[j] > min(object$lvs[, which.lvs])) {
            text(x = opt[j], y = maximum, labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 3) # should adjust "adj" rather than adding 0.1 to the maximum.
            segments(x0 = opt[j], x1 = opt[j], y0 = min(mu), y1 = maximum, lty = "dashed", col = cols[j])
          } else if (opt[j] < min(object$lvs[, which.lvs])) {
            text(x = min(object$lvs[, which.lvs]), y = apply(mu, 2, max)[j], labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 4)
          } else if (opt[j] > max(object$lvs[, which.lvs])) {
            text(x = max(object$lvs[, which.lvs]), y = apply(mu, 2, max)[j], labels = colnames(mu)[j], col = cols[j], cex = cex.spp, pos = 2)
          }
        }
      }
    }
    if (s.labels == "numeric") {
      text(x = object$lvs[, which.lvs], y = -1, labels = 1:nrow(object$y), col = s.colors)
    } else if (s.labels == "rug") {
      rug(object$lvs[, which.lvs], col = s.colors)
    }
  } else if (length(which.lvs) > 1 & object$num.lv > 1) {
    resid.cov <- object$params$theta[, which.lvs, drop = F]^2 + 2 * object$params$theta[, -c(1:object$num.lv), drop = F][, which.lvs, drop = F]^2
    largest.lnorms <- order(rowSums(resid.cov), decreasing = TRUE)[ind.spp]
    optima <- -object$params$theta[, 1:object$num.lv, drop = F][largest.lnorms, which.lvs, drop = F] / (2 * object$params$theta[largest.lnorms, -c(1:object$num.lv),
      drop = F
    ][, which.lvs, drop = F])
    cols <- (grDevices::rainbow(ncol(object$y) + 1)[2:(ncol(object$y) + 1)])[largest.lnorms]
    if (is.null(colnames(object$y))) row.names(optima) <- largest.lnorms
    if (!is.null(colnames(object$y))) row.names(optima) <- colnames(object$y)[largest.lnorms]

    if (class(object) == "gllvm.quadratic") quadr.coef <- object$params$theta[, -c(1:object$num.lv), drop = F][largest.lnorms, which.lvs, drop = F]

    # quadr.coef[which(round(quadr.coef, 2) == 0)] <- 0
    excl <- sapply(1:nrow(optima), function(j) any(optima[j, ] > 10 | optima[j, ] < -10))

    if (opt.region == "distribution") optSD <- 1 / sqrt(-2 * quadr.coef)
    if (opt.region == "confidence") optSD <- object$sd$optima[, which.lvs] # should also include covariances of the parameters..i need to do this properly still
    # need to sort optSD by largestlnorms
    # and take ind.spp
    if (opt.region == "confidence") {
      optSD <- optSD[largest.lnorms, ]
    }
    if (length(which(excl)) != 0) {
      if (length(which(excl)) == ncol(object$y)) {
        stop("Optima are too far removed from the latent varaibles to visualize")
      } else {
        message("There are species that lack optima on one (or multiple) of the chosen latent variables, which will thus not be visualized.",
          sep = " "
        )
        optima <- optima[!excl, , drop = F]
        if (opt.region != F) optSD <- optSD[!excl, , drop = F]

        cols <- cols[!excl]
      }
    }
    lvs <- object$lvs[, which.lvs]


    if (scale == TRUE) {
      optima <- sweep(optima, 2, ((getResidualCov(object)$trace.q + getResidualCov(object)$trace.q2) / sum(getResidualCov(object)$trace.q))[which.lvs], "*")
      lvs <- sweep(lvs, 2, ((getResidualCov(object)$trace.q + getResidualCov(object)$trace.q2) / sum(getResidualCov(object)$trace.q))[which.lvs], "*")
      if (opt.region != F) optSD <- sweep(optSD, 2, ((getResidualCov(object)$trace.q + getResidualCov(object)$trace.q2) / getResidualCov(object)$trace)[which.lvs], "*")
    }

    if (opt.region == F) {
      plot(rbind(optima, lvs),
        xlab = paste("LV", which.lvs[1]),
        ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...
      )
    } else if (opt.region %in% c("confidence", "distribution")) {
      lower <- optima + qnorm(level) * optSD # need to adapt this, not the right size at the moment
      upper <- optima + qnorm(1 - level) * optSD
      if (any(!apply(cbind(lower, upper), 1, function(x) all(x > -100 & x < 100)))) {
        flag <- T
      } else {
        flag <- F
      }
      optima <- optima[apply(cbind(lower, upper), 1, function(x) all(x > -100 & x < 100)), ]
      cols <- cols[apply(cbind(lower, upper), 1, function(x) all(x > -100 & x < 100))]
      optSD <- optSD[apply(cbind(lower, upper), 1, function(x) all(x > -100 & x < 100)), ]
      lower <- optima + qnorm(level) * optSD # need to adapt this, not the right size at the moment
      upper <- optima + qnorm(1 - level) * optSD
      if (opt.region == "confidence") {
        if (is.null(ylim)) {
          ylim <- range(c(rbind(upper[apply(abs(optima) < optSD, 1, all), ], lower[apply(abs(optima) < optSD, 1, all), ])[, 2], lvs[, 2]))
          if (any(ylim < (-20) | ylim > 20)) {
            ylim <- ifelse(ylim < (-20), -10, ylim)
            ylim <- ifelse(ylim > 20, 10, ylim)
          }
        }
        if (is.null(xlim)) {
          # zoom in on optima that have decently sized SD. If the SD are too large, we're not really interested.
          xlim <- range(c(rbind(upper[apply(abs(optima) < optSD, 1, all), ], lower[apply(abs(optima) < optSD, 1, all), ])[, 1], lvs[, 1]))

          if (any(xlim < (-20) | xlim > 20)) {
            xlim <- ifelse(xlim < (-20), -10, xlim)
            xlim <- ifelse(xlim > 20, 10, xlim)
          }
        }

        plot(NA,
          xlim = xlim, ylim = ylim, xlab = paste("LV", which.lvs[1]),
          ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...
        )
      }
      if (opt.region == "distribution") {
        if (is.null(xlim)) xlim <- range(c(rbind(upper, lower)[, 1], lvs[, 1]))
        if (is.null(ylim)) ylim <- range(c(rbind(upper, lower)[, 2], lvs[, 2]))
        plot(NA,
          xlim = xlim, ylim = ylim, xlab = paste("LV", which.lvs[1]),
          ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...
        )
      }
    }

    if (site.region) {
      if (length(col.ellips) != n) {
        col.ellips <- rep(col.ellips, n)
      }

      sdb <- sdA(object)
      object$A <- sdb + object$A
      ellipse <- gllvm:::ellipse
      for (i in 1:n) {
        # covm <- diag(diag(object$A[i,which.lvs,which.lvs]));
        covm <- object$A[i, which.lvs, which.lvs]
        ellipse(object$lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(ifelse(length(level) > 1, level[1], level), df = object$num.lv)), col = scales::alpha(col.ellips[i], ifelse(length(alpha.col) > 1, alpha.col[1], alpha.col)), lwd = ifelse(length(lwd.ellips) > 1, lwd.ellips[1], lwd.ellips), lty = ifelse(length(lty.ellips) > 1, lty.ellips[1], lty.ellips))
      }
    }


    abline(v = 0, h = 0, lty = "dotted")

    if (is.null(row.names(object$y))) {
      row.names(object$y) <- 1:nrow(object$y)
    }

    if (s.labels != FALSE) {
      text(lvs, labels = row.names(object$y), col = s.colors)
    }
    text(optima, labels = row.names(optima), col = cols, cex = cex.spp)
    if (opt.region == "distribution") {
      ellipse <- gllvm:::ellipse
      for (j in 2:nrow(optima)) {
        ellipse(optima[j, ], covM = diag(optSD[j, ]), rad = sqrt(qchisq(ifelse(length(level) > 1, level[2], level), df = object$num.lv)), col = scales::alpha(cols[j], ifelse(length(alpha.col) > 1, alpha.col[2], alpha.col)), lty = ifelse(length(lty.ellips) > 1, lty.ellips[2], lty.ellips), lwd = ifelse(length(lwd.ellips) > 1, lwd.ellips[2], lwd.ellips)) # these ignore scaling for now
        # need to draw a line for the species optima that are fixed at zero
        # car::ellipse(c(optima[j, 1], optima[j, 2]), s, env.range[j, ], center.pch = NULL, col=scales::alpha(cols[j], 0.7), lty = "dashed", lwd=1)
      }
    } else if (opt.region == "confidence") {
      ellipse <- gllvm:::ellipse
      for (j in 1:nrow(optima)) {
        if (!all(optima[j, ] == 0)) {
          if (any(optima[j, ] == 0)) {
            if (optima[j, 1] == 0) {
              segments(x0 = 0, x1 = 0, y0 = lower[1, 2], y1 = upper[1, 2], col = scales::alpha(cols[j], ifelse(length(alpha.col) > 1, alpha.col[2], alpha.col)), lty = ifelse(length(lty.ellips) > 1, lty.ellips[2], lty.ellips), lwd = ifelse(length(lwd.ellips) > 1, lwd.ellips[2], lwd.ellips))
            } else {
              segments(y0 = 0, y1 = 0, x0 = lower[1, 1], x1 = upper[1, 1], col = scales::alpha(cols[j], ifelse(length(alpha.col) > 1, alpha.col[2], alpha.col)), lty = ifelse(length(lty.ellips) > 1, lty.ellips[2], lty.ellips), lwd = ifelse(length(lwd.ellips) > 1, lwd.ellips[2], lwd.ellips))
            }
          } else {
            ellipse(optima[j, ], covM = diag(optSD[j, ]), rad = sqrt(qchisq(ifelse(length(level) > 1, level[2], level), df = object$num.lv)), col = scales::alpha(cols[j], ifelse(length(alpha.col) > 1, alpha.col[2], alpha.col)), lty = ifelse(length(lty.ellips) > 1, lty.ellips[2], lty.ellips), lwd = ifelse(length(lwd.ellips) > 1, lwd.ellips[2], lwd.ellips)) # these ignore scaling for now
          }
        }
      }
    }
  } else {
    stop("Not enough LVs for an optiplot")
  }
}


#' @export optiplot
optiplot <- function(object, ...) {
  UseMethod(generic = "optiplot")
}
