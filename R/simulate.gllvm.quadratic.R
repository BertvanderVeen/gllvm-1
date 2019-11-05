#' @title Simulate data from gllvm fit
#' @description Generate new data using the fitted values of the parameters
#'
#' @param object an object of class 'gllvm'.
#' @param nsim an optional positive integer specifying the number of simulated datasets. Defaults to 1.
#' @param seed an optional integer to set seed number, passed to set.seed. Defaults to a random seed number.
#' @param conditional if \code{FALSE} the simulation is performed marginally over the latent variables.
#' @param ... not used.
#'
#' @details
#' simulate function for gllvm objects.  Note this simulates marginally over LVs.
#' an option is to add a conditional argument, which would fix the LVs.
#' David Warton--change this
#' 
#' @return A matrix containing generated data.
#' @author David Warton, Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'  \donttest{
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- scale(antTraits$env[, 1:3])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, family = poisson())
#'# Simulate data
#'newdata <- simulate(fit)
#'}
#'@aliases simulate
#'@method simulate gllvm.quadratic
#'@export

simulate.gllvm.quadratic = function(object, nsim = 1, conditional = FALSE, seed = NULL, ...) {
    # code chunk from simulate.lm to sort out the seed thing:
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv) else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    
    nRows = dim(object$lvs)[1]
    nCols = dim(object$params$theta)[1]
    if(conditional == FALSE){
    # generate new latent variables
    lvsNew = matrix(rnorm(nsim * nRows * object$num.lv), ncol = object$num.lv)
    }else{
      lvsNew <- object$lvs
    }
    if (is.null(object$X)) {
        prs = predict.gllvm.quadratic(object, newLV = lvsNew, type = "response")
    } else prs = predict.gllvm.quadratic(object, newX = object$X[rep(1:nRows, nsim), ], newLV = lvsNew, type = "response")
    
    # generate new data
    nTot = nsim * nRows * nCols  # total number of values to generate
    if (object$family == "negative.binomial") 
        phis = matrix(rep(object$params$phi, each = nsim * nRows), ncol = nCols)
    if (object$family == "gaussian") 
        phis = matrix(rep(object$params$phi, each = nsim * nRows), ncol = nCols)
    if(object$family == "ordinal"){
      sims = matrix(0, nrow = nsim * nRows, ncol = nCols)
        for(j in 1:nCols){
        k <- unique(object$y[,j])
        for(i in 1:(nsim * nRows)){
            sims[i,j] <- sample(k,1,prob=prs[,i,j][!is.na(prs[,i,j])])
          }
        }
      dimnames(prs)[[3]] <- colnames(object$y)
      dimnames(prs)[[2]] <- 1:(nsim * nRows)
      prs <- prs[1,,]
      
    }
    newDat = switch(object$family, binomial = rbinom(nTot, size = 1, prob = prs), poisson = rpois(nTot, prs), negative.binomial = rnbinom(nTot, 
        size = phis*prs, mu = prs), ordinal = sims, stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
    # reformat as data frame with the appropriate labels
    newDat = as.data.frame(matrix(newDat, ncol = nCols))
    dimnames(newDat) = dimnames(prs)
    return(newDat)
}

#' @export simulate
simulate <- function(object, ...) {
    UseMethod(generic = "simulate")
}
