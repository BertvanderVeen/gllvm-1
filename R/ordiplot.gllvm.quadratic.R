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
      plot(1:n, object$lvs[, which.lvs], ylab = paste("LV", which.lvs), xlab = "Row index")
    }
    
    if (length(which.lvs) > 1) {
      testcov <- predict(object, LVonly = T) #LV effect only
      do.svd <- svd(testcov, length(which.lvs), length(which.lvs))
      choose.lvs <- do.svd$u * matrix(do.svd$d[1:length(which.lvs)]^alpha, nrow = n, ncol = length(which.lvs), byrow = TRUE)
      choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:length(which.lvs)]^(1 - alpha), nrow = p, ncol = length(which.lvs), byrow = TRUE)
      
      
      if (!biplot) {
        sdd <- diag(sqrt(diag(cov(object$lvs[, which.lvs]))), nrow = length(which.lvs))
        choose.lvs <- scale(choose.lvs) %*% sdd
        plot(choose.lvs[, which.lvs], xlab = paste("LV", which.lvs[1]), ylab = paste("LV", which.lvs[2]), 
             main = main, type = "n", ...)
        if (!jitter) 
          if (symbols) {
            points(choose.lvs[, which.lvs], col = s.colors)
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
        
        plot(rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs]), xlab = paste("LV", which.lvs[1]), 
             ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...)
        
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
          if(is.null(object$sd)){
            cat("Cannot plot prediction regions if no standard errors were calculated.")
          }else{
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
      
      if(object$family=="ordinal")type="link"
      #mu <- predict(object, newLV = newLV, LVonly = T, which.lvs = which.lvs,type = type,intercept=intercept)[,largest.lnorms,drop=F]
      #have to rebuiltin the backtransformation of the linear predictor now I included curve
      mu<-predict(object, LVonly = T, which.lvs = which.lvs,type = type,intercept=intercept)[,largest.lnorms,drop=F]
      
      if(legend==F){
        #pdf(NULL)
        #plot(NA, xlim = c(range(object$lvs)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", ...)
        #maxstr<-max(strwidth(colnames(mu)))*1.25  
        #invisible(dev.off())
        #previpously range below was range(object$lvs)+maxstr, might still want to add this again
        plot(NA, xlim = c(c(min(object$lvs)-.1,max(object$lvs)+.1)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", ...)
      }else{
        plot(NA, xlim = c(range(object$lvs)), ylim = range(mu), ylab = "Predicted ", xlab = paste("LV", which.lvs, sep = " "), xaxs = "i", ...)
      }
      
      if(legend==T){
        legend(x="topleft",text.col=cols,colnames(mu), cex=cex.spp)
      }
      if(type=="response"){
      if(object$family=="binomial"){
        linkinv<-pnorm
      }else if(object$family%in%c("poisson", "negative.binomial")){
        linkinv <- exp
      }else{
        stop("Ordinal model plot not yet implemented on response scale")
      }
      }
      for (j in 1:ncol(mu)) {
        if(intercept==F){
          if(type=="link"){func <-function(x,u,u2)x*u+x^2*u2}else{func <-function(x,u,u2)linkinv(x*u+x^2*u2)}
          curve(func(x,u=object$params$theta[largest.lnorms,,drop=F][j,which.lvs],u2=object$params$theta[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j],add=T)  
        }else{
          if(type=="link"){func <-function(x,beta,u,u2)beta+x*u+x^2*u2}else{func <-function(x,beta,u,u2)linkinv(beta+x*u+x^2*u2)}
          curve(func(x,beta=object$params$beta0[largest.lnorms][j],u=object$params$theta[largest.lnorms,,drop=F][j,which.lvs],u2=object$params$theta[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j],add=T)  
        }
        #lines(x = sort(newLV[, 1]), y = mu[order(newLV[, 1]), j], col = cols[j])
        if(legend==F){
          opt <- summary(object)$Optima[largest.lnorms,which.lvs]
          if(opt[j]<max(object$lvs)&opt[j]>min(object$lvs)){
            text(x = opt[j], y = apply(mu,2,max)[j], labels = colnames(mu)[j], col = cols[j], cex=cex.spp, adj=c(0.5,-0.5))    
          }else if(opt[j]<min(object$lvs)){
            text(x = min(object$lvs), y = apply(mu,2,max)[j], labels = colnames(mu)[j], col = cols[j], cex=cex.spp, adj=c(0,-0.5))    
          }else if(opt[j]>max(object$lvs))
            text(x = max(object$lvs), y = apply(mu,2,max)[j], labels = colnames(mu)[j], col = cols[j], cex=cex.spp, adj=c(1,-0.5))    
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
          message("Species ", paste(paste(row.names(optima[excl, ,drop=F]), collapse = ", "), " optima too far removed from the latent variables and will not be visualized (if bell=T)", 
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
      
      if(env.ranges==F){
        plot(rbind(optima, lvs), xlab = paste("LV", which.lvs[1]), 
             ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...)
      }
      if(env.ranges==T){
        env.lower <- optima - 1.96 * tolerances#need to adapt this, not the right size at the moment
        env.upper <- optima + 1.96 * tolerances
        env.range <- env.upper - env.lower
        if(any(!apply(env.range,1,function(x)all(x>-100&x<100)))){
          flag<-T
        }else{
          flag<-F
        }
        optima<-optima[apply(env.range,1,function(x)all(x>-100&x<100)),]
        cols<-cols[apply(env.range,1,function(x)all(x>-100&x<100))]
        env.range<-env.range[apply(env.range,1,function(x)all(x>-100&x<100)),]
        xlim<-range(c(rbind(optima+env.range,optima-env.range)[,which.lvs[1]],lvs[which.lvs[[1]]]))
        ylim<-range(c(rbind(optima+env.range,optima-env.range)[,which.lvs[2]],lvs[which.lvs[[2]]]))
        plot(NA, xlim=xlim,ylim=ylim,xlab = paste("LV", which.lvs[1]), 
             ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...)
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
          ellipse(lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)), col="gray", lty="solid")#these ignore scaling for now
        }
      }
      
      
      abline(v=0,h=0,lty="dotted")
      
      if(is.null(row.names(object$y))){
        row.names(object$y)<-1:nrow(object$y)
      }
      

      text(lvs, labels = row.names(object$y),col="gray")
      text(optima, labels = row.names(optima), col = cols, cex=cex.spp)
      if(env.ranges==T){
        for (j in 1:nrow(optima)) {
          s = diag(2)
          ellipse(optima[j,which.lvs], covM = diag(env.range[j,]), rad = sqrt(qchisq(level, df=object$num.lv)), col=scales::alpha(cols[j], 0.7), lty="dashed")#these ignore scaling for now
          #car::ellipse(c(optima[j, 1], optima[j, 2]), s, env.range[j, ], center.pch = NULL, col=scales::alpha(cols[j], 0.7), lty = "dashed", lwd=1)
        }
        if(flag)message("Some tolerances are too large to plot.")
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

