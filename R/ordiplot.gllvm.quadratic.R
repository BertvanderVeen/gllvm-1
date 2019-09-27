#' @title Plot latent variables from gllvm model
#' @description Plots latent variables and their corresponding coefficients (biplot).
#'
#' @param object   an object of class 'gllvm'.
#' @param biplot   \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot. The default is none, or all if \code{biplot = TRUE}.
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
#' @param jitter   if \code{TRUE}, jittering is applied on points.
#' @param jitter.amount   numeric, positive value indicating an amount of jittering for each point, defaults to 0.2 (jitter range).
#' @param s.colors colors for sites
#' @param symbols logical, if \code{TRUE} sites are plotted using symbols, if \code{FALSE} (default) site numbers are used
#' @param cex.spp size of species labels in biplot
#' @param scale For 2D plots, either "species" or "sites" to scale optima or site scores by the ratio variance explained. Alternatively can be "tolerances" to scale optima by tolerances and site scores by average tolerances per latent variable.
#' @param bell logical, if TRUE plots bell-shapes (1D) or biplot with (scaled) predicted optima and 95% of the predicted environmental ranges (2D)
#' @param env.ranges logical, if TRUE plots predicted species distributions in 2D
#' @param type when predicting bell-shapes, can be used to predict on the response or link scale. Default is response.
#' @param intercept when predicting bell-shapes, can be used to include species-intercepts in the plot. Default is TRUE
#' @param legend when \code{TRUE} adds legend in the topleft corner of the plot, instead of species names in the plot
#' @param ...\tadditional graphical arguments.
#'
#' @details
#' Function constructs a scatter plot of two latent variables, i.e. an ordination plot. If only one latent
#' variable is in the fitted model, latent variables are plotted against their row indices \code{bell=F}.
#' or if \code{bell=T} against their marginal predictions on the link scale.
#'
#' Coefficients related to latent variables are plotted in the same figure with the latent
#' variables if \code{biplot = TRUE, bell = F}. They are labeled using the column names of y. The number
#' of latent variable coefficients to be plotted can be controlled by ind.spp. An argument alpha
#' is used to control the relative scaling of the latent variables and their coefficients.
#' If \code{alpha = 0.5}, the latent variables and their coefficients are on the same scale.
#'
#' If \code{bell=T} and two latent variables are selected (which.lvs), species optima are plotted with predicted environmental ranges.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui
#'
#' @examples
#' #'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'fit <- gllvm(y, family = poisson())
#'# Ordination plot:
#'ordiplot(fit)
#'# Biplot with 10 species
#'ordiplot(fit, biplot = TRUE, ind.spp = 10)
#'
#'@aliases ordiplot ordiplot.gllvm.quadratic
#'@export
#'@export ordiplot.gllvm.quadratic
ordiplot.gllvm.quadratic <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = NULL, jitter = FALSE, 
    jitter.amount = 0.2, s.colors = 1, symbols = FALSE, cex.spp = 0.7, bell = TRUE,env.ranges=FALSE, type = "response", intercept = TRUE, legend=FALSE,scale="species",...) {
    if (any(class(object) != "gllvm.quadratic")) 
        stop("Class of the object isn't 'gllvm.quadratic'.")
    a <- jitter.amount
    n <- NROW(object$y)
    p <- NCOL(object$y)
    if (!is.null(ind.spp)) {
        ind.spp <- min(c(p, ind.spp))
    } else {
        ind.spp <- p
    }
    if (object$num.lv == 0) 
        stop("No latent variables to plot.")
    if (is.null(which.lvs)) {
        if (object$num.lv > 1) {
            which.lvs <- 1:2
        } else {
            which.lvs <- 1
        }
    }
    if (is.null(rownames(object$params$theta))) 
        rownames(object$params$theta) = paste("V", 1:p)
    
    if (bell == F) {
        if (length(which.lvs) == 1) {
            plot(1:n, object$lvs[, which.lvs], ylab = "LV1", xlab = "Row index")
        }
        
        if (length(which.lvs) > 1) {
            testcov <- predict(object, LVonly = T)
            do.svd <- svd(testcov, length(which.lvs), length(which.lvs))
            choose.lvs <- do.svd$u * matrix(do.svd$d[1:length(which.lvs)]^alpha, nrow = n, ncol = length(which.lvs), byrow = TRUE)
            choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:length(which.lvs)]^(1 - alpha), nrow = p, ncol = length(which.lvs), byrow = TRUE)
            
            
            if (!biplot) {
                sdd <- diag(sqrt(diag(cov(object$lvs[, which.lvs]))), nrow = length(which.lvs))
                choose.lvs <- scale(choose.lvs) %*% sdd
                plot(choose.lvs[, which.lvs], xlab = paste("Latent variable ", which.lvs[1]), ylab = paste("Latent variable ", which.lvs[2]), 
                  main = main, type = "n", ...)
                if (!jitter) 
                  if (symbols) {
                    points(choose.lvs[, which.lvs], col = s.colors, ...)
                  } else {
                    text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
                  }
                if (jitter) 
                  if (symbols) {
                    points(choose.lvs[, which.lvs][, 1] + runif(n, -a, a), choose.lvs[, which.lvs][, 2] + runif(n, -a, a), col = s.colors, 
                      ...)
                  } else {
                    text((choose.lvs[, which.lvs][, 1] + runif(n, -a, a)), (choose.lvs[, which.lvs][, 2] + runif(n, -a, a)), label = 1:n, 
                      cex = 1.2, col = s.colors)
                  }
            }
            
            if (biplot) {
                largest.lnorms <- order(apply(object$params$theta[, -c(which.lvs, which.lvs * 2)]^2, 1, sum), decreasing = TRUE)[1:ind.spp]
                
                plot(rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs]), xlab = paste("Latent variable ", which.lvs[1]), 
                  ylab = paste("Latent variable ", which.lvs[2]), main = main, type = "n", ...)
                
                if (!jitter) {
                  if (symbols) {
                    points(choose.lvs[, which.lvs], col = s.colors, ...)
                  } else {
                    text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
                  }
                  text(matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)), label = rownames(object$params$theta)[largest.lnorms], 
                    col = 4, cex = cex.spp)
                }
                if (jitter) {
                  if (symbols) {
                    points(choose.lvs[, which.lvs[1]] + runif(n, -a, a), (choose.lvs[, which.lvs[2]] + runif(n, -a, a)), col = s.colors, 
                      ...)
                  } else {
                    text((choose.lvs[, which.lvs[1]] + runif(n, -a, a)), (choose.lvs[, which.lvs[2]] + runif(n, -a, a)), label = 1:n, 
                      cex = 1.2, col = s.colors)
                  }
                  text((matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)) + runif(2 * length(largest.lnorms), 
                    -a, a)), label = rownames(object$params$theta)[largest.lnorms], col = 4, cex = cex.spp)
                }
            }
            
        }
    } else {
        if (length(which.lvs) == 1) {
            if (object$num.lv == 1) {
                quadr.coef <- object$params$theta[1:ind.spp, -1,drop=F]
            } else {
                quadr.coef <- object$params$theta[,-c(1:object$num.lv),drop=F][1:ind.spp, which.lvs,drop=F]
            }
            
            quadr.coef[which(round(quadr.coef, 3) == 0)] <- 0
            
            newLV <- matrix(NA, nrow = 1000, ncol = length(which.lvs))
            newLV[, 1] <- seq(from = min(object$lvs[, which.lvs]), max(object$lvs[, which.lvs]), length.out = 1000)
            
            mu <- predict(object, newLV = newLV, LVonly = T, which.lvs = which.lvs,type = type,intercept=intercept)[,1:ind.spp,drop=F]
            
            if(legend==F){
              pdf(NULL)
              plot(NA, xlim = c(min(newLV), max(newLV)), ylim = range(mu), ylab = "(marginal) Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
              maxstr<-max(strwidth(colnames(mu)))*1.25  
              invisible(dev.off())
                                                                                                             
              plot(NA, xlim = c(min(newLV), max(newLV)+maxstr), ylim = range(mu), ylab = "(marginal) Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
            }else{
              plot(NA, xlim = c(min(newLV), max(newLV)), ylim = range(mu), ylab = "(marginal) Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
            }
            cols <- (grDevices::rainbow(ncol(mu) + 1)[2:(ncol(mu) + 1)])
            if(legend==T){
              legend(x="topleft",text.col=cols,colnames(mu))
            }
            
            for (j in 1:ncol(mu)) {
                lines(x = sort(newLV[, 1]), y = mu[order(newLV[, 1]), j], col = cols[j])
                col <- col2rgb(cols[j], alpha = TRUE)
                col[4] <- 127
                col <- rgb(col[1], col[2], col[3], col[4], maxColorValue = 255)
                if(legend==F){
                  text(x = max(newLV[, 1]), y = tail(mu[order(newLV[, 1]), j])[1], labels = colnames(mu)[j], col = cols[j], adj = 0)  
                }
            }
            abline(v = 0, h = 1, col = "black", lty = "dashed")
            text(x = object$lvs[, which.lvs], y = range(mu)[1], labels = 1:nrow(object$y), col = "grey")
        } else if (length(which.lvs) > 1 & object$num.lv > 1) {
            optima <- -object$params$theta[, 1:object$num.lv, drop = F][1:ind.spp, which.lvs, drop = F]/(2 * object$params$theta[1:ind.spp, -c(1:object$num.lv), 
                drop = F][, which.lvs, drop = F])
            quadr.coef <- object$params$theta[, -c(1:object$num.lv), drop = F][1:ind.spp, which.lvs, drop = F]
            quadr.coef[which(round(quadr.coef, 3) == 0)] <- 0
            excl <- which(sapply(1:nrow(quadr.coef), function(j) any(quadr.coef[j, ] == 0)))
            
            if (any(excl)) {
                message(paste("The species", paste(row.names(quadr.coef[excl, ]), collapse = ", "), "lack a quadratic response on the chosen LV(s) and won't be plotted.", 
                  sep = " "))
                optima <- optima[-excl, , drop = F]
                quadr.coef <- quadr.coef[-excl, , drop = F]
                
            }
            lvs <- object$lvs
            tolerances <- 1/sqrt(-2 * quadr.coef)
            if(scale=="species"){
              optima <- optima/(getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q))
              env.upper <-  env.upper / (getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q))
              env.lower <-  env.lower / (getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q))
            }else if(scale=="sites"){
              lvs<-lvs/(getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q))
            }else if (scale == "tolerances") {
              optima <- optima/tolerances
              lvs <- lvs/apply(tolerances, 2, mean)
              tolerances <- tolerances/tolerances
            }

            
            env.lower <- optima - 1.96 * tolerances
            env.upper <- optima + 1.96 * tolerances
            
            if(env.ranges==F){
              plot(rbind(optima, lvs), xlab = paste("Latent variable ", which.lvs[1]), 
                   ylab = paste("Latent variable ", which.lvs[2]), main = main, type = "n", ...)
            }else{
              plot(rbind(rbind(env.lower, env.upper), rbind(env.lower, env.upper)), xlab = paste("Latent variable ", which.lvs[1]), 
                   ylab = paste("Latent variable ", which.lvs[2]), main = main, type = "n", ...)
            }
            
            
            if(is.null(row.names(object$y))){
              row.names(object$y)<-1:nrow(object$y)
            }
            col <- grDevices::rainbow(nrow(tolerances))
            
            
            text(lvs, labels = row.names(object$y),col="gray")
            text(optima, labels = row.names(optima), col = col)
            env.range <- env.upper - env.lower
            if(env.ranges==T){
              for (j in 1:nrow(optima)) {
                s = diag(2)
                car::ellipse(c(optima[j, 1], optima[j, 2]), s, env.range[j, ], center.pch = NULL, col = col[j], lty = "dashed")
              }
            }

            
            
        } else {
            stop("Not enough LVs for a biplot")
        }
    }
}


#'@export ordiplot
ordiplot <- function(object, ...) {
    UseMethod(generic = "ordiplot")
}

