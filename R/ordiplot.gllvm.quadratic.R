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
#' @param bell logical, if TRUE plots bell-shapes (1D) or biplot with (scaled) predicted optima and 95 percent of the predicted environmental ranges (2D)
#' @param env.ranges logical, if TRUE plots predicted species distributions in 2D
#' @param type when predicting bell-shapes, can be used to predict on the response or link scale. Default is response (except for ordinal, for which the only option is "link").
#' @param intercept when predicting bell-shapes, can be used to include species-intercepts in the plot. Default is TRUE
#' @param legend when \code{TRUE} adds legend in the topleft corner of the plot, instead of species names in the plot
#' @param predict.region logical, if \code{TRUE} prediction regions for the predicted latent variables are plotted, defaults to \code{FALSE}.
#' @param leve level for prediction regions.
#' @param ...\tadditional graphical arguments.
#'
#' @details
#' Function constructs a scatter plot of two latent variables, i.e. an ordination plot. If only one latent
#' variable is in the fitted model, latent variables are plotted against their row indices \code{bell=F}.
#' or if \code{bell=T} against their predictions.
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
                                     jitter.amount = 0.2, s.colors = 1, symbols = FALSE, cex.spp = 0.7, bell = TRUE,env.ranges=FALSE, type = "response", intercept = TRUE, legend=FALSE,scale="species", predict.region = FALSE, level = 0.95, ...) {
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
      testcov <- predict(object, LVonly = T) #LV effect only
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
        resid.cov <- object$params$theta[,which.lvs,drop=F]^2 + 2*object$params$theta[,-c(1:object$num.lv),drop=F][,which.lvs,drop=F]^2
        largest.lnorms <- order(rowSums(resid.cov), decreasing = TRUE)[1:ind.spp]
        
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
        if(predict.region){
          sdb<-sdA(object)
          object$A<-sdb+object$A
          r=0
          #if(object$row.eff=="random") r=1
          
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- diag(object$A[i,which.lvs+r]);
            } else {
              #covm <- diag(diag(object$A[i,which.lvs,which.lvs]));
              covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)))
          }
        }
      }
      
    }
  } else {
    if (length(which.lvs) == 1) {
      resid.cov <- object$params$theta[,which.lvs,drop=F]^2 + 2*object$params$theta[,-c(1:object$num.lv),drop=F][,which.lvs,drop=F]^2
      largest.lnorms <- order(rowSums(resid.cov), decreasing = TRUE)[1:ind.spp]
      cols <- (grDevices::rainbow(ncol(object$y) + 1)[2:(ncol(object$y) + 1)])[largest.lnorms]
      if (object$num.lv == 1) {
        quadr.coef <- object$params$theta[largest.lnorms, -1,drop=F]
      } else {
        quadr.coef <- object$params$theta[,-c(1:object$num.lv),drop=F][largest.lnorms, which.lvs,drop=F]
      }
      
      quadr.coef[which(round(quadr.coef, 3) == 0)] <- 0
      
      newLV <- matrix(NA, nrow = 1000, ncol = length(which.lvs))
      newLV[, 1] <- seq(from = min(object$lvs[, which.lvs]), max(object$lvs[, which.lvs]), length.out = 1000)
      if(object$family=="ordinal")type="link"
      mu <- predict(object, newLV = newLV, LVonly = T, which.lvs = which.lvs,type = type,intercept=intercept)[,largest.lnorms,drop=F]
      
      if(legend==F){
        pdf(NULL)
        plot(NA, xlim = c(min(newLV), max(newLV)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
        maxstr<-max(strwidth(colnames(mu)))*1.25  
        invisible(dev.off())
        
        plot(NA, xlim = c(min(newLV), max(newLV)+maxstr), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
      }else{
        plot(NA, xlim = c(min(newLV), max(newLV)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i")
      }
      
      if(legend==T){
        legend(x="topleft",text.col=cols,colnames(mu))
      }
      
      for (j in 1:ncol(mu)) {
        lines(x = sort(newLV[, 1]), y = mu[order(newLV[, 1]), j], col = cols[j])
        if(legend==F){
          text(x = max(newLV[, 1]), y = tail(mu[order(newLV[, 1]), j])[1], labels = colnames(mu)[j], col = cols[j], adj = 0)  
        }
      }
      text(x = object$lvs[, which.lvs], y = range(mu)[1], labels = 1:nrow(object$y), col = "grey")
    } else if (length(which.lvs) > 1 & object$num.lv > 1) {
      resid.cov <- object$params$theta[,which.lvs,drop=F]^2 + 2*object$params$theta[,-c(1:object$num.lv),drop=F][,which.lvs,drop=F]^2
      largest.lnorms <- order(rowSums(resid.cov), decreasing = TRUE)[1:ind.spp]
      optima <- -object$params$theta[, 1:object$num.lv, drop = F][largest.lnorms, which.lvs, drop = F]/(2 * object$params$theta[largest.lnorms, -c(1:object$num.lv), 
                                                                                                                                drop = F][, which.lvs, drop = F])
      cols <- (grDevices::rainbow(ncol(object$y) + 1)[2:(ncol(object$y) + 1)])[largest.lnorms]
      if(is.null(colnames(object$y)))row.names(optima)<-largest.lnorms
      if(!is.null(colnames(object$y)))row.names(optima)<-colnames(object$y)[largest.lnorms]
      
      
      quadr.coef <- object$params$theta[, -c(1:object$num.lv), drop = F][largest.lnorms, which.lvs, drop = F]
      # quadr.coef[which(round(quadr.coef, 2) == 0)] <- 0
      excl <- sapply(1:nrow(optima), function(j) any(optima[j, ] > 10 | optima[j, ] < -10))
      
      if (length(which(excl))!=0) {
        if(length(which(excl))==ncol(object$y)){
          stop("Optima are too far removed from the latent varaibles to visualize")
        }else{
          message(paste("Columns", paste(row.names(optima[excl, ]), collapse = ", "), "have optima too far removed from the latent variables and will not be visualized (if bell=T)", 
                        sep = " "))
          optima <- optima[!excl, , drop = F]
          quadr.coef <- quadr.coef[!excl, , drop = F] 
          cols <- cols[!excl]
        }
        
      }
      lvs <- object$lvs
      tolerances <- 1/sqrt(-2 * quadr.coef)
      
      if(scale=="species"){
        optima <- sweep(optima,2,(getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q)),"*")
        tolerances <-  sweep(tolerances,2,(getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q)),"*")
      }else if(scale=="sites"){
        lvs<-sweep(lvs,2, (getResidualCov(object)$trace.q/sum(getResidualCov(object)$trace.q)),"*")
      }else if (scale == "tolerances") {
        optima <- optima/tolerances
        lvs <- lvs/apply(tolerances, 2, mean)
        tolerances <- tolerances/tolerances
      }
      if(predict.region){
        sdb<-sdA(object)
        object$A<-sdb+object$A
        r=0
        #if(object$row.eff=="random") r=1
        
        for (i in 1:n) {
          if(object$Lambda.struc == "diagonal"){
            covm <- diag(object$A[i,which.lvs+r]);
          } else {
            #covm <- diag(diag(object$A[i,which.lvs,which.lvs]));
            covm <- object$A[i,which.lvs+r,which.lvs+r];
          }
          ellipse(object$lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)))#these ignore scaling for now
        }
      }
      env.lower <- optima - 1.96 * tolerances
      env.upper <- optima + 1.96 * tolerances
      if(env.ranges==F){
        plot(rbind(optima, lvs), xlab = paste("Latent variable ", which.lvs[1]), 
             ylab = paste("Latent variable ", which.lvs[2]), main = main, type = "n", ...)
      }else{
        env.range <- env.upper - env.lower
        xlim<-range(c(rbind(optima+env.range,optima-env.range)[,which.lvs[1]],lvs[which.lvs[[1]]]))
        ylim<-range(c(rbind(optima+env.range,optima-env.range)[,which.lvs[2]],lvs[which.lvs[[2]]]))
        plot(NA, xlim=xlim,ylim=ylim,xlab = paste("Latent variable ", which.lvs[1]), 
             ylab = paste("Latent variable ", which.lvs[2]), main = main, type = "n", ...)
      }
      abline(v=0,h=0,lty="dotted")
      
      if(is.null(row.names(object$y))){
        row.names(object$y)<-1:nrow(object$y)
      }

      text(lvs, labels = row.names(object$y),col="gray")
      text(optima, labels = row.names(optima), col = cols)
      if(env.ranges==T){
        for (j in 1:nrow(optima)) {
          s = diag(2)
          car::ellipse(c(optima[j, 1], optima[j, 2]), s, env.range[j, ], center.pch = NULL, col = cols[j], lty = "dashed")
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

