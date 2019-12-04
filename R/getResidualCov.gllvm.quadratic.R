#' @title Extract residual covariance matrix from gllvm object
#' @description  Calculates the residual covariance matrix for gllvm model.
#'
#' @param object  an object of class 'gllvm'.
#' @param adjust  The type of adjustment used for negative binomial and binomial distribution when computing residual correlation matrix. Options are 0 (no adjustment) and 1 (the default adjustment), see details.
#'
#' @return Function returns following components:
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of the residual covariance matrix}
#'  \item{trace.q }{trace of the residual covariance matrix per latent variable}
#'
#' @details 
#' Residual covariance matrix, storing information on species co-occurrence that is not explained by the environmental variables (if included), is calculated using the matrix of latent variables loadings, that is,  \Sigma_{j,k} =\eqn{\Theta_j\Theta_k' + 2diag(D_j)diag(\D_k)'} + .
#' 
#' When the responses are modelled using the negative binomial distribution, the residual variances for each species must be adjusted for overdispersion. This is done by the term \eqn{\psi^{(1)}(1/\phi_j)} (\code{adjust = 1}), where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' 
#' The residual covariance matrix with \code{adjust = 1} can be obtained by using Poisson-Gamma parametrization
#' \deqn{Y_{ij} \sim Poisson(\mu_{ij} \lambda_j),}
#' where \eqn{\lambda_j \sim Gamma(\phi_j, \phi_j)} and \eqn{\mu_{ij}} is as above. The mean and the variance are of similar form as above and we have that
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = \theta_j'\theta_j + 2D_jD_j' + \psi^{(1)}(\phi_j),}
#' where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' In the case of binomial distribution, the adjustment terms (\code{adjust = 1}) are 1 for probit link.
#' These are obtained by treating binomial model as latent variable model. Assume
#' \deqn{Y^*_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} \sim N(0, 1)} for probit model, and \eqn{e_{ij} ~ logistic(0, 1)} for logit model.
#' Then binary response is defined as \eqn{Y_{ij} = 1}, if \eqn{Y^*_{ij} > 0} and 0 otherwise.
#' Now we have that \eqn{\mu_{ij} = P(Y_{ij} = 1) = P(Y^*_{ij} > 0) = P(\eta_{ij} > -e_{ij}) = P(e_{ij} <= \eta_{ij})} which leads to the probit model.
#' On linear predictor scale we then have that
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For the probit model, j,k = 1\ldots p \Sigma_{j,k} \eqn{\Theta_j\Theta_k' + 2diag(D_j)diag(\D_k)'+ I_m}.
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual covariance:
#'rescov <- getResidualCov(fit)
#'rescov$cov
#'# Trace of the covariance matrix
#'rescov$tr
#'# Trace per latent variable
#'rescov$trace.q
#'
#'@aliases getResidualCov getResidualCov.gllvm.quadratic
#'@method getResidualCov gllvm.quadratic
#'@export
#'@export getResidualCov.gllvm.quadratic
getResidualCov.gllvm.quadratic = function(object, adjust = 1) {
    ResCov <- object$params$theta[, 1:object$num.lv] %*% t(object$params$theta[, 1:object$num.lv]) + 2 * object$params$theta[, -c(1:object$num.lv)] %*% 
        t(object$params$theta[, -c(1:object$num.lv)])
    ResCov.q <- sapply(1:object$num.lv, function(q) object$params$theta[, q] %*% t(object$params$theta[, q]) + 2 * object$params$theta[, 
        q + object$num.lv] %*% t(object$params$theta[, q + object$num.lv]), simplify = F)
    if (adjust > 0 && object$family %in% c("negative.binomial", "binomial")) {
        if (object$family == "negative.binomial") {
            ResCov <- ResCov + diag(trigamma(object$params$phi))  #adjusted for gamma parameterization with phi rather than inv.phi
            ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + diag(trigamma(object$params$phi))/object$num.lv, simplify = F)
        }
        
        if (object$family == "binomial") {
            if (object$link == "probit") 
                ResCov <- ResCov + diag(ncol(object$y))
            ResCov.q <- sapply(1:object$num.lv, function(q) ResCov.q[[q]] + diag(ncol(object$y))/object$num.lv, simplify = F)
        }
    }
    ResCov.q <- sapply(1:object$num.lv, function(q) sum(diag(ResCov.q[[q]])))
    names(ResCov.q) <- paste("LV", 1:object$num.lv, sep = "")
    colnames(ResCov) <- colnames(object$y)
    rownames(ResCov) <- colnames(object$y)
    out <- list(cov = ResCov, trace = sum(diag(ResCov)), trace.q = ResCov.q)
    return(out)
}

#'@export getResidualCov
getResidualCov <- function(object, adjust) {
    UseMethod(generic = "getResidualCov")
}
