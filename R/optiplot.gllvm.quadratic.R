  #' @title Plot species optima for latent variables from a gllvm
  #' @description Plots latent variables and their corresponding species optima.
  #'
  #' @param object   an object of class 'gllvm'.
  #' @param ind.spp  the number of response variables (usually, species) to include on the biplot (sorted by varianve explained). The default is all.
  #' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
  #' @param main  main title.
  #' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
  #' @param s.colors colors for sites
  #' @param cex.spp size of species labels in biplot
  #' @param scale For 2D plots, either "FALSE",species" or "sites" to scale optima or site scores by the ratio variance explained. Alternatively can be "tolerances" to scale optima by tolerances and site scores by average tolerances per latent variable.
  #' @param opt.ranges Only for 2D plots, efaults to FALSE. If "statistical", plots statistical uncertainties for species optima. If "environmental" plots preicted environmental ranges
  #' @param type Can be used to predict on the response or link scale. Default is response (except for ordinal, for which the only option is "link").
  #' @param intercept Can be used to include species-intercepts in the plot. Default is TRUE
  #' @param legend when \code{TRUE} adds legend in the topleft corner of the plot, instead of species names in the plot
  #' @param predict.region logical, if \code{TRUE} prediction regions for the predicted latent variables are plotted, defaults to \code{FALSE}.
  #' @param level level for prediction regions of sites and optima.
  #' @param alpha.col used to control the transparency of confidence interval ribbons in 1D plot
  #' @param ... additional graphical arguments.
  #'
  #' @details
  #' Function constructs a plot of species optima.
  #'
  #'Plots species optima in one or two dimensions, potentially with predicted environmental ranges.
  #'
  #' @author Bert van der Veen
  #'
  #' @examples
  #' #'## Load a dataset from the mvabund package
  #'data(antTraits)
  #'y <- as.matrix(antTraits$abund)
  #'fit <- gllvm.quadratic(y, family = poisson())
  #'# Ordination plot:
  #'optiplot(fit)
  #'# Optiplot with 10 species
  #'optiplot(fit, ind.spp = 10)
  #'
  #'@aliases optiplot optiplot.gllvm.quadratic
  #'@export
  #'@export optiplot.gllvm.quadratic
  optiplot.gllvm.quadratic <- function(object,  ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = NULL, 
                                       s.colors = 1, cex.spp = 0.7, opt.ranges=FALSE, type = "response", intercept = TRUE, legend=FALSE,scale=FALSE, predict.region = FALSE, level = 0.95, alpha.col = 0.4, ...) {
    if(class(object)!="gllvm.quadratic")
      stop("Class of the object isn't 'gllvm.quadratic'. linear GLLVM not implemented yet.")
    
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
    if(length(which.lvs)>object$num.lv)stop("More latent variables select than included in the model.")
    if (is.null(rownames(object$params$theta))) 
      rownames(object$params$theta) = paste("V", 1:p)
    
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
        if(opt.ranges=="statistical"){
          thetaCIlower <- object$params$theta - 1.96*object$sd$theta
          thetaCIupper <- object$params$theta + 1.96*object$sd$theta
          if(intercept==F){
            curveLowCI<-curve(func(x,u=thetaCIlower[largest.lnorms,,drop=F][j,which.lvs],u2=thetaCIlower[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j], add=T, lty="dashed")  
            curveUpCI<-curve(func(x,u=thetaCIupper[largest.lnorms,,drop=F][j,which.lvs],u2=thetaCIupper[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j], add=T, lty="dashed")  
          }else{
            curveLowCI<-curve(func(x,beta=object$params$beta0[largest.lnorms][j],u=thetaCIlower[largest.lnorms,,drop=F][j,which.lvs],u2=thetaCIlower[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j], lty="dashed", add=T)  
            curveUpCI<-curve(func(x,beta=object$params$beta0[largest.lnorms][j],u=thetaCIupper[largest.lnorms,,drop=F][j,which.lvs],u2=thetaCIupper[largest.lnorms,,drop=F][j,-(1:object$num.lv),drop=F][,which.lvs]),col=cols[j], lty="dashed", add=T)  
          }
          
          cols2 <- scales::alpha(cols[j], alpha.col)
          polygon(c(curveLowCI$x,rev(curveLowCI$x)),c(curveLowCI$y,rev(curveUpCI$y)),border=NA,col=cols2)
        }
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
      
      if(class(object)=="gllvm.quadratic")quadr.coef <- object$params$theta[, -c(1:object$num.lv), drop = F][largest.lnorms, which.lvs, drop = F]
      
      # quadr.coef[which(round(quadr.coef, 2) == 0)] <- 0
      excl <- sapply(1:nrow(optima), function(j) any(optima[j, ] > 10 | optima[j, ] < -10))
      
      if(opt.ranges=="environmental")optSD <- 1/sqrt(-2 * quadr.coef)
      if(opt.ranges=="statistical")optSD <- object$sd$optima[,which.lvs]
      #need to sort optSD by largestlnorms
      #and take ind.spp
      if(opt.ranges!=F){
        optSD <- optSD[largest.lnorms,]
      }
      if (length(which(excl))!=0) {
        if(length(which(excl))==ncol(object$y)){
          stop("Optima are too far removed from the latent varaibles to visualize")
        }else{
          message("Species ", paste(paste(row.names(optima[excl, ,drop=F]), collapse = ", "), " optima too far removed from the latent variables and will not be visualized (if bell=T)", 
                                    sep = " "))
          optima <- optima[!excl, , drop = F]
          if(opt.ranges!=F)optSD <- optSD[!excl, , drop = F] 
          
          cols <- cols[!excl]
        }
        
      }
      lvs <- object$lvs
      
      
      if(scale=="species"){
        optima <- sweep(optima,2,((getResidualCov(object)$trace.q+getResidualCov(object)$trace.q2)/sum(getResidualCov(object)$trace.q)),"*")
        if(opt.ranges!=F)optSD <-  sweep(optSD,2,((getResidualCov(object)$trace.q+getResidualCov(object)$trace.q2)/getResidualCov(object)$trace),"*")
      }else if(scale=="sites"){
        lvs<-sweep(lvs,2, ((getResidualCov(object)$trace.q+getResidualCov(object)$trace.q2)/sum(getResidualCov(object)$trace.q)),"*")
      }
      #else if (scale == "tolerances"&opt.ranges!=FALSE) {
      #  optima <- optima/tolerances
      #  lvs <- lvs/apply(tolerances, 2, mean)
      #  tolerances <- tolerances/tolerances
      #}
      
      if(opt.ranges==F){
        plot(rbind(optima, lvs), xlab = paste("LV", which.lvs[1]), 
             ylab = paste("LV", which.lvs[2]), main = main, type = "n", ...)
      }
      
      if(opt.ranges%in%c("statistical","environmental")){
        lower <- optima - 1.96 * optSD#need to adapt this, not the right size at the moment
        upper <- optima + 1.96 * optSD
        if(any(!apply(cbind(lower,upper),1,function(x)all(x>-100&x<100)))){
          flag<-T
        }else{
          flag<-F
        }
        optima<-optima[apply(cbind(lower,upper),1,function(x)all(x>-100&x<100)),]
        cols<-cols[apply(cbind(lower,upper),1,function(x)all(x>-100&x<100))]
        optSD<-optSD[apply(cbind(lower,upper),1,function(x)all(x>-100&x<100)),]
        xlim<-range(c(rbind(upper,lower)[,which.lvs[1]],lvs[which.lvs[[1]]]))
        ylim<-range(c(rbind(upper,lower)[,which.lvs[2]],lvs[which.lvs[[2]]]))
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
      if(opt.ranges!=F){
        for (j in 1:nrow(optima)) {
          ellipse(optima[j,which.lvs], covM = diag(optSD[j,]), rad = sqrt(qchisq(level, df=object$num.lv)), col=scales::alpha(cols[j], alpha.col), lty="dashed")#these ignore scaling for now
          #need to draw a line for the species optima that are fixed at zero
          #car::ellipse(c(optima[j, 1], optima[j, 2]), s, env.range[j, ], center.pch = NULL, col=scales::alpha(cols[j], 0.7), lty = "dashed", lwd=1)
        }
        if(flag)message("Some ranges are too large to plot.")
      }
    } else {
      stop("Not enough LVs for an optiplot")
    }
    
  }
  
  
  #'@export optiplot
  optiplot <- function(object, ...) {
    UseMethod(generic = "optiplot")
  }