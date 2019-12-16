#' @title Niche overlap
#' @description  Calculates the overlap between niches of species
#'
#' @param object   An object of class 'gllvm.quadratic'.
#' @param which.lvs The latent variables niche overlap should be calculated for. 
#' @param spp Species index, if specified calculates the probability of species occuring jointly.
#' @param option either 1 or 2. If 1, calculates the matrix of pairwise occurence probabilities. If 2 calculates the niche overlap coefficient as specified in Mutshinda and O' Hara (2011). Essentially, the probability of joint occurrence for species i and j, relative to the probability of the niche belonging to species i.
#' @details
#' The niche overlap is calculated per latent variable. If multiple latent variables are selected, calculates the average niche overlap.
#' 
#' @author Bert van der Veen
#' @references
#' Inman, H.F., Bradley Jr, E.L. (1989). The overlapping coefficient as a measure of agreement between probability distributions and point estimation of the overlap of two normal densities. Communications in Statistics - Theory and methods, 18:3851-3874.
#' Mutshinda C.M., O'Hara, R.B. (2011) Integrating the niche and neutral perspectives on community structure and dynamics. Oecologia 166: 241-251.
#' 
#' @examples
#'# Load a dataset from the mvabund package
#'data(spider)
#'# Fit quadratic gllvm
#'fit <- gllvm.quadratic(y = y, family = poisson(),num.lv=2)
#'# calculate the niche overlap for the first latent variable
#' overlap <- niche.overlap(fit, which.lvs = 1)
#' 
#'# plot results
#'library(corrplot)
#'par(mfrow=c(1,2))
#'corrplot(overlap, cl.lim = c(0,1))
#'ordiplot(fit,which.lvs=1)
#'
#'# plot residual correlations of species with >20% niche overlap in two dimensions
#'corrplot(gllvm.quadratic::getResidualCor(test),p.mat=1-niche.overlap(test,which.lvs=1:2),sig.level=.8)
#'
#'# calculate probability of multiple species occuring together in two dimensions
#'niche.overlap(test,which.lvs=1:2,spp=9:12)
#'
#'# calculate probability of a species occurring in the combined latent space
#'niche.overlap(test,which.lvs=1:2,spp=1)
#'@aliases niche.overlap niche.overlap.gllvm.quadratic
#'@method niche.overlap gllvm.quadratic
#'@export

niche.overlap.gllvm.quadratic <- function(object,which.lvs=1, spp=NULL, option=1){
  if(class(object)!="gllvm.quadratic")stop("Can only calculate niche overlap for an object of class 'gllvm.quadratic'")
  p <- ncol(object$y)
  if(!is.null(spp)){
   nicheO <- integrate(niche.integrand2, -Inf, Inf, spp=spp, which.lvs=which.lvs, object=object)$value
  }else{
  
  nicheO <- matrix(0,ncol=p, nrow=p)
  if(option==1){
    for(j in p:1){
      for(j2 in 1:(j-1)){
        if(j!=j2&j>1){
          overlap<-numeric(length(which.lvs))
          for(q in which.lvs){
            opt1<-summary(object)$Optima[j,q]
            opt2<-summary(object)$Optima[j2,q]
            tol2<-summary(object)$Tol[j2,q]
            tol1<-summary(object)$Tol[j,q]
            overlap[q] <- integrate(niche.integrand,-Inf, Inf, opt1=opt1, opt2=opt2, tol1=tol1, tol2=tol2)$value
          }
          overlap[is.nan(overlap)|is.infinite(overlap)]<-0
          
          nicheO[j,j2]<-sum(overlap)/length(which.lvs)
        }
      }
    }  
    nicheO[upper.tri(nicheO)]<-t(nicheO)[upper.tri(t(nicheO))]
  }else{
    for(j in 1:p){
      for(j2 in 1:p){
        if(j!=j2){
          overlap<-numeric(length(which.lvs))
          for(q in which.lvs){
            opt1<-summary(object)$Optima[j,q]
            opt2<-summary(object)$Optima[j2,q]
            tol2<-summary(object)$Tol[j2,q]
            tol1<-summary(object)$Tol[j,q]
            overlap[q] <- integrate(niche.integrand,-Inf, Inf, opt1=opt1, opt2=opt2, tol1=tol1, tol2=tol2)$value/integrate(niche.integrand,-Inf, Inf, opt1=opt1, tol1=tol1, opt2=opt1, tol2=tol1)$value
          }
          overlap[is.nan(overlap)|is.infinite(overlap)]<-0
          
          nicheO[j,j2]<-sum(overlap)/length(which.lvs)
        }
      }
    }  
  }
  diag(nicheO)<-1
  colnames(nicheO) <- row.names(nicheO) <- colnames(object$y)
  }
  return(nicheO)
}

#corrplot(niche.overlap2(test,1:2),cl.lim=c(0,1))
#' @export niche.overlap
niche.overlap <- function(object, ...) {
  UseMethod(generic = "niche.overlap")
}

#Pairwise probability of occuring jointly
niche.integrand <- function(x,opt1,opt2,tol1,tol2){
  f1 <- dnorm(x, mean=opt1, sd=tol1)
  f2 <- dnorm(x, mean=opt2, sd=tol2)
  result <- f1*f2
  return(result)
}

#probability of spp occuring jointly across LVs
niche.integrand2 <- function(x,object, which.lvs, spp){
  result<-1
  spp <- unique(spp)
  for(j in spp){
    for(q in which.lvs){
    result<-result*dnorm(x, mean=summary(object)$Optima[j,q], sd=summary(object)$Tol[j,q])
    }
  }
  return(result)
}

#FOREACH CAN USE OPENMP. ALSO ON WINDOWS??

