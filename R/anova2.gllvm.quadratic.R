#' @title Analysis Of Deviance for quadratic and linear gllvms
#' @description  Computes an analysis of deviance table for two generalized linear latent variable model fits.
#'
#' @param object   an object of class 'gllvm'.
#' @param ...   one or more objects of class 'gllvm'
#'
#' @details
#' Computes likelihood-ratio test for two or more gllvm models.
#' Test results makes sense only for nested models.
#' Notice also that this test is not designed for testing models which have
#' degrees of freedom difference larger than 20. For such models
#' the P-value should be treated as very approximate.
#'
#' @author Jenni Niku
#'
#' @examples
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- antTraits$abund
#'X <- antTraits$env
#'TR <- antTraits$traits
#'# Fit gllvm model
#'fit1 <- gllvm(y, X, TR, formula = ~ Bare.ground + Shrub.cover, family = poisson())
#'fit2 <- gllvm(y, X, TR, formula = ~ Bare.ground + Shrub.cover +
#'              (Bare.ground + Shrub.cover) : Webers.length, family = poisson())
#'# Test if the model with fourth corner interaction terms is significantly
#'# better using likelihood-ratio test:
#'anova(fit1, fit2)
#'
#'@aliases anova2 anova2.gllvm.quadratic
#'@method anova2 gllvm.quadratic
#'@export

anova2.gllvm.quadratic <- function(object, ...) {
  objects <- list(object, ...)
  if (length(objects) < 2)
    stop("At least two objects are needed for tests.")
  if (!any((sapply(objects, class) %in% c("gllvm.quadratic"))))
    stop("The function 'anova.gllvm.quadratic' needs to include at least one object of class gllvm.quadratic.")
  if (any(!(sapply(objects, class) %in% c("gllvm.quadratic","gllvm"))))
    stop("Unknown class object in function.")
  
  tt <- sapply(objects, function(x)
    x$method)
  if (!(all(tt == "VA") == !all(tt == "LA")))
    stop("The objects are not comparable when they are fitted using different methods.")

  y <- object$y
  n <- NROW(y)
  p <- NCOL(y)
  diff <- sapply(objects, function(x) sum(x$y - y))
  if (any(!(diff == 0)))
    stop("The objects can not be compared")


  df.list <- sapply(objects, function(x)
    attr(logLik(x), "df"))
  objects_order <- objects[order(df.list)]
  formulas <- sapply(objects_order, function(x)
    formula(x$terms))

  #if(test=="LR"){
  df.list <- sapply(objects_order, function(x) attr(logLik(x), "df"))
  ll.list <- sapply(objects_order, logLik)

  D <- 2 * (ll.list[-1] - ll.list[1:(length(df.list) - 1)])
  df.chisq <- (df.list[-1] - df.list[1:(length(df.list) - 1)])
  Pval <- 1 - pchisq(D, df.chisq)
  paste("Model", 1:length(objects_order))
  result <- data.frame( Resid.Df = n * p - df.list, D = c(0, D), Df.diff = c(0, df.chisq), P.value = c("", signif(Pval)) )
  #}
  if (any(result$Df > 20))
    warning( "This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.\n")
  for (i in 1:length(objects_order)) {
    formchar <- as.character(formulas[[i]])
    if (length(formchar) == 3)
      formchar <- formchar[c(2, 1, 3)]
    cat("Model ", i, ": ", formchar, "\n")
  }
  return(result)
}

#' @export anova2
anova2 <- function(object, ...)
{
  UseMethod(generic = "anova2")
}
