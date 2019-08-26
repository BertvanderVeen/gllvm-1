#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' 
#'@aliases getLV getLV.gllvm.quadratic
#'@method getLV gllvm.quadratic
#'
#'@export getLV.gllvm.quadratic

getLV.gllvm.quadratic <- function(object)
{
  return(object$lvs)
}

#'@export getLV

getLV <- function(object)
{
  UseMethod(generic = "getLV")
}