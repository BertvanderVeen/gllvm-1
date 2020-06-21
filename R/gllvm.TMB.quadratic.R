########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku, Bert van der Veen
##########################################################################################
gllvm.TMB.quadratic <- function(y, X = NULL, formula = NULL, num.lv = 2, family = "poisson",
                                Lambda.struc = "unstructured", row.eff = FALSE, reltol = 1e-10, trace = FALSE, trace2 = FALSE,
                                seed = NULL, maxit = 2000, start.lvs = NULL, offset = NULL, sd.errors = TRUE,
                                n.init = 1, start.params = NULL,
                                optimizer = "optim", starting.val = "res", diag.iter = 1,
                                Lambda.start = c(0.1, 0.5), jitter.var = 0, par.scale = 1, fn.scale = 1, zeta.struc = "species", maxit.lingllvm = NULL, starting.val.lingllvm = "res", common.tolerances = FALSE, parallel = FALSE, start.struc = "species", gamma1 = 0, gamma2 = 0, theta4 = NULL, Lambda2.start = 0.01) {
  n <- dim(y)[1]
  p <- dim(y)[2]
  tr <- NULL
  num.lv <- num.lv
  y <- as.matrix(y)
  formula1 <- formula
  if (any(diag.iter > 1)) stop("'diag.iter' can only be 0 or 1.") # Lambda.start at 0.1 will cause issues as that's what lambda2 is sometimes started at. Then A-D is 0.
  if (!is.numeric(y)) {
    stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  }
  if (is.null(rownames(y))) {
    rownames(y) <- paste("Row", 1:n, sep = "")
  }
  if (is.null(colnames(y))) {
    colnames(y) <- paste("Col", 1:p, sep = "")
  }
  if (family == "binomial") {
    link <- "probit"
  }
  if (num.lv == 0) {
    stop("Need to at least include one LV in the model. If you want to fit a GLM, please use gllvm() instead.")
  }
  if (family == "ordinal") {
    y00 <- y
    if (min(y) == 0) {
      y <- y + 1
    }
    max.levels <- apply(y, 2, function(x) length(min(x):max(x)))
    if (any(max.levels == 1) & zeta.struc == "species" || all(max.levels == 2) & zeta.struc == "species") {
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
    }

    if (any(!apply(y, 2, function(x) all(diff(sort(unique(x))) == 1))) & zeta.struc == "species") {
      stop("Can't fit ordinal model if there are species with missing classes. Please reclassify per species or use zeta.struc = `common` ")
    }

    if (any(diff(sort(unique(c(y)))) != 1) & zeta.struc == "common") {
      stop("Can't fit ordinal model if there are missing classes. Please reclassify.")
    }
  }
  num.X <- 0

  if (!is.null(X)) {
    if (!is.null(formula)) {
      xb <- as.matrix(model.matrix(formula, data = data.frame(X)))
      X <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      colnames(X) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      Xd <- X1 <- X

      num.X <- dim(X)[2]
    } else {
      n1 <- colnames(X)
      formula <- paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula <- paste(formula, n1[i1], sep = "+")
        }
      }
      formula <- formula(formula)
      xb <- as.matrix(model.matrix(formula, data = data.frame(X)))
      X <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      num.X <- dim(X)[2]
      colnames(X) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      Xd <- X1 <- X

      nxd <- colnames(Xd)
      formulab <- paste("~", nxd[1], sep = "")

      for (i in 2:length(nxd)) {
        formulab <- paste(formulab, nxd[i], sep = "+")
      }
      formula1 <- formulab
    }
  }
  if (is.null(formula) && is.null(X)) {
    formula <- "~ 1"
  }
  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if (!is.null(seed)) {
    set.seed(seed)
  }

  out <- list(y = y, X = X, logL = Inf, X.design = X)
  old.logL <- Inf
  if (starting.val == "lingllvm") {
    if (length(diag.iter) > 1) {
      diag.iter2 <- diag.iter[2]
      diag.iter <- diag.iter[1]
    } else {
      diag.iter2 <- diag.iter
      diag.iter <- 0
    }
    if (length(jitter.var) > 1) {
      jitter.var2 <- jitter.var[1]

      jitter.var <- jitter.var[2]
    } else {
      jitter.var2 <- jitter.var
      jitter.var <- 0
    }
    if (length(n.init) == 1) {
      n.init2 <- n.init
      n.init <- 1
    } else {
      n.init2 <- n.init[1]
      n.init <- n.init[2]
    }
  } else {
    if (length(n.init) > 1) n.init <- n.init[2]
    if (length(jitter.var) > 1) {
      jitter.var2 <- jitter.var[1]

      jitter.var <- jitter.var[2]
    } else {
      jitter.var2 <- jitter.var
      jitter.var <- 0
    }
    if (length(diag.iter) > 1) diag.iter <- diag.iter[2]
  }



  if (is.null(theta4)) {
    theta4 <- rep(0, num.lv)
  } else if (length(theta4) == 1) {
    theta4 <- rep(theta4, num.lv)
  } else if (length(theta4) != num.lv) {
    stop("Wrong length theta4 supplied.")
  }

  if (starting.val == "lingllvm") {
    if (is.null(maxit.lingllvm)) {
      maxit.lingllvm <- maxit
    }
    if (trace2) cat("Running linear gllvm..")
    fit <- gllvm(y, formula = formula, X = X, num.lv = num.lv, family = family, row.eff = row.eff, n.init = n.init2, maxit = maxit.lingllvm, reltol = reltol, optimizer = optimizer, diag.iter = diag.iter2, jitter.var = jitter.var2, starting.val = starting.val.lingllvm, Lambda.start = Lambda.start, Lambda.struc = Lambda.struc, method = "VA", sd.errors = FALSE, offset = offset, zeta.struc = zeta.struc, seed = seed)
    start.params <- fit
    if (trace2) cat("Done! \n")
  }
  if (n.init[1] > 1) seed <- sample(1:10000, n.init)

  # helper function for parallel optimization
  makeMod <- function(i) {
    sigma <- 1
    if (starting.val != "lingllvm") {
      start.values.gllvm.TMB <- gllvm:::start.values.gllvm.TMB
      fit <- start.values.gllvm.TMB(y = y, X = X, TR = NULL, family = family, offset = offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[i], starting.val = starting.val, jitter.var = jitter.var2, row.eff = row.eff, zeta.struc = zeta.struc)
    }
    if (is.null(start.params)) {
      beta0 <- fit$params[, 1]
      betas <- NULL
      if (!is.null(X)) {
        betas <- c(fit$params[, 2:(num.X + 1)])
      }
      lambdas <- NULL

      lambdas <- as.matrix(fit$params[, (ncol(fit$params) - num.lv + 1):ncol(fit$params)])
      lambdas[upper.tri(lambdas)] <- 0
      fit$params <- cbind(fit$params, matrix(Lambda2.start, ncol = num.lv, nrow = p))
      row.params <- NULL

      if (row.eff != FALSE) {
        row.params <- fit$row.params
        if (row.eff == "random") {
          sigma <- sd(row.params) # 1;#
        }
      } # rep(0,n)
      lvs <- NULL
      lvs <- matrix(fit$index, ncol = num.lv)
      # if(!is.null(X)){
      #   #lambdas <- fit$params[,(ncol(fit$params) - num.lv*2 + 1):(ncol(fit$params)-num.lv)]
      #   if(common.tolerances==F)lambda2 <- fit$params[,(ncol(fit$params)-num.lv+1):ncol(fit$params)]
      # }else if(is.null(X)){
      #   #lambdas <- fit$params[,(ncol(fit$params) - num.lv*2 + 1):(ncol(fit$params)-num.lv)]
      #   if(common.tolerances==F)lambda2 <- fit$params[,-c(1:(num.lv+1))]
      # }

      # subtract a fraction from the 0 quadratic scores, otherwise the optimization can't get away from the 0s where necessary.
    } else {
      if (dim(start.params$y) == dim(y) &&
        is.null(X) == is.null(start.params$X) &&
        (row.eff == start.params$row.eff)) {
        if (start.params$family == "ordinal") {
          if (start.params$zeta.struc == "species") zeta <- start.params$params$zeta[, -1]
          if (start.params$zeta.struc == "common") zeta <- start.params$params$zeta[-1]
        }
        if (start.params$family %in% c("negative.binomial", "gaussian", "gamma")) {
          phi <- start.params$params$phi
        }
        beta0 <- start.params$params$beta0 ## column intercepts
        betas <- NULL
        if (!is.null(X)) {
          if (!(all(dim(X) == dim(start.params$X)))) stop("Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that predictors X are the same in both models.")
        }
        betas <- c(start.params$params$Xcoef) ## covariates coefficients
        lambdas <- NULL
        lambda2 <- NULL
        lambdas <- start.params$params$theta[, 1:num.lv]
        lambdas[upper.tri(lambdas)] <- 0
        if (class(start.params) == "gllvm.quadratic") {
          lambda2 <- start.params$params$theta[, -c(1:num.lv)]
        } else {
          if (starting.val != "lingllvm") {
            if (!is.null(X)) {
              lambda2 <- fit$params[, (ncol(fit$params) - num.lv + 1):ncol(fit$params)]
            } else if (is.null(X)) {
              lambda2 <- fit$params[, -c(1:(num.lv + 1))]
            } # this still needs to be implement for traits.
          } else {
            lambda2 <- matrix(0, nrow = p, ncol = num.lv)
            for (j in 1:p) {
              for (q in 1:num.lv) {
                lambda2[j, q] <- -.5 / (sum((start.params$lvs[, q] - (sum(y[, j] * start.params$lvs[, q]) / sum(y[, j])))^2 * y[, j]) / sum(y[, j]))
              }
            }
            if (any(is.infinite(lambda2))) {
              lambda2[is.infinite(lambda2)] <- -0.5
            }
            fit$params$theta <- cbind(fit$params$theta, lambda2)
          }
        }
        row.params <- NULL
        if (start.params$row.eff != FALSE) {
          row.params <- start.params$params$row.params
          if (row.eff == "fixed") {
            row.params[1] <- 0
          }
          if (row.eff == "random") {
            sigma <- start.params$params$sigma
          }
        } ## row parameters
        lvs <- NULL
        lvs <- matrix(start.params$lvs, ncol = num.lv)
      } else {
        stop("Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.")
      }
    }
    phis <- NULL
    if (family %in% c("gaussian", "gamma")) {
      phi <- fit$phi
    }

    if (family == "negative.binomial") {
      phis <- fit$phi
      if (starting.val != "lingllvm") {
        if (any(phis > 25)) phis[phis > 25] <- 25
        if (any(phis < 0.04)) phis[phis < 0.04] <- 0.04
        fit$phi <- phis
        phis <- 1 / phis
      }
    }
    if (is.null(start.params) | family != "ordinal") {
      if (family == "ordinal") {
        K <- max(y00) - min(y00)
        if (zeta.struc == "species") {
          zeta <- c(t(fit$zeta[, -1]))
          zeta <- zeta[!is.na(zeta)]
        } else {
          zeta <- fit$zeta[-1]
        }
      } else {
        zeta <- 0
      }
    }

    if (is.null(offset)) {
      offset <- matrix(0, nrow = n, ncol = p)
    }

    if (!is.null(row.params)) {
      r0 <- row.params
    } else {
      r0 <- rep(0, n)
    }
    a <- c(beta0)
    b <- NULL
    if (!is.null(X)) b <- matrix(betas, ncol(X), p, byrow = TRUE)
    # diag(lambdas) <- log(diag(lambdas)) !!!
    lambda <- lambdas[lower.tri(lambdas, diag = TRUE)]
    u <- lvs

    if (!is.null(phis)) {
      phi <- phis
    } else {
      phi <- rep(1, p)
      fit$phi <- phi
    }

    optr <- NULL
    timeo <- NULL
    se <- NULL

    if (is.null(start.params) || start.params$method != "VA") {
      if (Lambda.struc == "diagonal" || diag.iter > 0) {
        Au <- log(rep(Lambda.start[1], num.lv * n)) # 1/2, 1
      } else {
        Au <- c(log(rep(Lambda.start[1], num.lv * n)), rep(0, num.lv * (num.lv - 1) / 2 * n)) # 1/2, 1
      }
    } else {
      Au <- NULL
      for (d in 1:num.lv) {
        if (start.params$Lambda.struc == "unstructured" || length(dim(start.params$A)) == 3) {
          Au <- c(Au, log(start.params$A[, d, d]))
        } else {
          Au <- c(Au, log(start.params$A[, d]))
        }
      }
      if (Lambda.struc != "diagonal" && diag.iter == 0) {
        if (start.params$Lambda.struc == "unstructured") {
          for (d1 in 1:num.lv) {
            for (d2 in 1:num.lv) {
              if (d1 > d2) { # only get unique off-diagonal
                Au <- c(Au, start.params$A[, d1, d2])
              }
            }
          }
        } else {
          Au <- c(Au, rep(0, num.lv * (num.lv - 1) / 2 * n)) ## this is redundant
        }
      }
    }
    if (start.struc != "common") {
      if (starting.val == "lingllvm") lambda2 <- matrix(-Lambda2.start, ncol = num.lv, nrow = p)
      if (!is.null(start.params)) {
        if (class(start.params) == "gllvm") lambda2 <- matrix(Lambda2.start, ncol = num.lv, nrow = p)
      }
      if (common.tolerances == TRUE) start.struc <- "common"
    }

    if (start.struc == "common") {
      lambda3 <- 0
      lambda2 <- matrix(Lambda2.start, ncol = num.lv, nrow = 1)
    } else {
      lambda2 <- matrix(Lambda2.start, ncol = num.lv, nrow = p)
      lambda3 <- rep(0, num.lv)
    }

    u <- u + mvtnorm::rmvnorm(n, rep(0, num.lv), diag(num.lv) * jitter.var) # mostly makes sense when using lingllvm

    # lambda2 <- matrix(apply(lambda2,2,mean),ncol=num.lv)#perhaps think about doing something else here still? i.e. start at..tol from lingllvm?
    if (length(Lambda.start) < 2) {
      Ar <- rep(1, n)
    } else {
      Ar <- rep(Lambda.start[2], n)
    }
    if (row.eff == FALSE) {
      xr <- matrix(0, 1, p)
    } else {
      xr <- matrix(1, 1, p)
    }
    if (!is.null(X)) {
      Xd <- cbind(1, X)
    } else {
      Xd <- matrix(1, n)
    }
    if (starting.val == "lingllvm" & row.eff == "random") sigma <- sigma[1]
    if (family == "poisson") {
      familyn <- 0
    }
    if (family == "negative.binomial") {
      familyn <- 1
    }
    if (family == "binomial") {
      familyn <- 2
    }
    if (family == "ordinal") {
      familyn <- 7
    }
    if (family == "gaussian") {
      familyn <- 3
    }
    if (family == "gamma") {
      familyn <- 4
    }

    if (starting.val != "zero" & start.struc != "common" | class(start.params) == "lingllvm" & start.struc != "common") {
      mp <- list(r0 = factor(rep(NA, length(r0))), b = factor(rep(NA, length(rbind(a, b)))), B = factor(rep(NA, 1)), lambda = factor(rep(NA, length(lambda))), lambda3 = factor(rep(NA, num.lv)), u = factor(rep(NA, length(u))), lg_phi = factor(rep(NA, length(phi))), log_sigma = factor(rep(NA, length(sigma))), Au = factor(rep(NA, length(Au))), lg_Ar = factor(rep(NA, length(Ar))), zeta = factor(rep(NA, length(zeta))))

      if (row.eff == "random") {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 1, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          map = mp, parameters = list(r0 = matrix(r0), b = rbind(a, b), B = matrix(0), lambda = lambda, lambda2 = t(lambda2), lambda3 = lambda3, u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, lg_Ar = log(Ar), zeta = zeta),
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        )
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 0, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          map = mp, parameters = list(r0 = matrix(r0), b = rbind(a, b), B = matrix(0), lambda = lambda, lambda2 = t(lambda2), lambda3 = lambda3, u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), zeta = zeta),
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        ) ## GLLVM
      }

      if (optimizer == "nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, iter.max = maxit, eval.max = maxit, trace = trace2)), silent = !trace2))
      }
      if (optimizer == "optim") {
        if (!is.null(par.scale)) {
          if (par.scale == "coef") {
            parscale <- abs(objr$par) # this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
            parscale[parscale == 0] <- 1
          } else if (is.numeric(par.scale)) {
            parscale <- rep(par.scale, length(objr$par))
          }
        } else {
          parscale <- rep(1, length(objr$par))
        }
        if (is.null(fn.scale) | !is.numeric(fn.scale)) {
          fnscale <- 1
        } else {
          fnscale <- fn.scale
        }
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr, method = "BFGS", control = list(reltol = reltol, maxit = maxit, parscale = parscale, fnscale = fnscale, trace = trace2), hessian = FALSE), silent = !trace2))
      }
      lambda2 <- matrix(optr$par, byrow = T, ncol = num.lv, nrow = p)
    }


    if (row.eff == "random") {
      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 1, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
        parameters = list(r0 = matrix(r0), b = rbind(a, b), B = matrix(0), lambda = lambda, lambda2 = t(lambda2), lambda3 = lambda3, u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, lg_Ar = log(Ar), zeta = zeta),
        inner.control = list(mgcmax = 1e+200, maxit = maxit),
        DLL = "qgllvm"
      )
    } else {
      objr <- TMB::MakeADFun(
        data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 0, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
        parameters = list(r0 = matrix(r0), b = rbind(a, b), B = matrix(0), lambda = lambda, lambda2 = t(lambda2), lambda3 = lambda3, u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), zeta = zeta),
        inner.control = list(mgcmax = 1e+200, maxit = maxit),
        DLL = "qgllvm"
      ) ## GLLVM
    }

    if (optimizer == "nlminb") {
      timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, iter.max = maxit, eval.max = maxit, trace = trace2)), silent = !trace2))
    }
    if (optimizer == "optim") {
      if (!is.null(par.scale)) {
        if (par.scale == "coef") {
          parscale <- abs(objr$par) # this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
          parscale[parscale == 0] <- 1
        } else if (is.numeric(par.scale)) {
          parscale <- rep(par.scale, length(objr$par))
        }
      } else {
        parscale <- rep(1, length(objr$par))
      }
      if (is.null(fn.scale) | !is.numeric(fn.scale)) {
        fnscale <- 1
      } else {
        fnscale <- fn.scale
      }
      timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr, method = "BFGS", control = list(reltol = reltol, maxit = maxit, parscale = parscale, fnscale = fnscale, trace = trace2), hessian = FALSE), silent = !trace2))
    }
    if (inherits(optr, "try-error")) warning(optr[1])
    # if(!inherits(optr,"try-error")&common.tolerances!=TRUE){
    if (diag.iter > 0 && Lambda.struc == "unstructured" && num.lv > 1 && !inherits(optr, "try-error")) {
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam == "r0"])
      b1 <- matrix(param1[nam == "b"], num.X + 1, p)
      lambda1 <- param1[nam == "lambda"]
      u1 <- matrix(param1[nam == "u"], n, num.lv)
      lg_phi1 <- param1[nam == "lg_phi"]
      log_sigma1 <- param1[nam == "log_sigma"]
      # previously  c(pmax(param1[nam=="Au"],rep(log(0.001), num.lv*n)), rep(0.01,num.lv*(num.lv-1)/2*n))
      # this line adds the covariance parameters after diag iter, it didn't start though, this does, sometimes.
      Au1 <- c(pmax(param1[nam == "Au"], rep(log(1e-4), num.lv * n)), rep(1e-4, num.lv * (num.lv - 1) / 2 * n)) # c(rep(0,length(param1[names(param1)=="Au"])), rep(0.01,num.lv*(num.lv-1)/2*n))
      # this is causing issues I think...
      lg_Ar1 <- param1[nam == "lg_Ar"]

      zeta <- param1[nam == "zeta"]

      if (common.tolerances == F & start.struc != "common") {
        lambda2 <- t(matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = p))
      } else {
        lambda2 <- t(matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = 1))
      }

      lambda3 <- param1[nam == "lambda3"]
      # lambda3 <- ifelse(lambda3<0.1,0.5,lambda3)#added this line for the binomial

      if (row.eff == "random") {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 1, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = log_sigma1, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        )
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 0, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        ) # GLLVM#
      }

      if (optimizer == "nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, iter.max = maxit, eval.max = maxit, trace = trace2)), silent = !trace2))
      }
      if (optimizer == "optim") {
        if (!is.null(par.scale)) {
          if (par.scale == "coef") {
            parscale <- abs(objr$par) # this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
            parscale[parscale == 0] <- 1
          } else if (is.numeric(par.scale)) {
            parscale <- rep(par.scale, length(objr$par))
          }
        } else {
          parscale <- rep(1, length(objr$par))
        }
        if (is.null(fn.scale) | !is.numeric(fn.scale)) {
          fnscale <- 1
        } else {
          fnscale <- fn.scale
        }
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr, method = "BFGS", control = list(reltol = reltol, maxit = maxit, parscale = parscale, fnscale = fnscale, trace = trace2), hessian = FALSE), silent = !trace2))
      }
      if (optimizer == "optim") {
        if (inherits(optr, "try-error") || is.nan(optr$value) || is.na(optr$value) || is.infinite(optr$value)) {
          optr <- optr1
          objr <- objr1
          Lambda.struc <- "diagonal"
        }
      } else {
        if (inherits(optr, "try-error") || is.nan(optr$objective) || is.na(optr$objective) || is.infinite(optr$objective)) {
          optr <- optr1
          objr <- objr1
          Lambda.struc <- "diagonal"
        }
      }
    }
    if (inherits(optr, "try-error")) warning(optr[1])
    if (common.tolerances == F & start.struc == "common") {
      objr1 <- objr
      optr1 <- optr
      param1 <- optr$par
      nam <- names(param1)
      r1 <- matrix(param1[nam == "r0"])
      b1 <- matrix(param1[nam == "b"], num.X + 1, p)
      lambda1 <- param1[nam == "lambda"]
      u1 <- matrix(param1[nam == "u"], n, num.lv)
      lg_phi1 <- param1[nam == "lg_phi"]
      log_sigma1 <- param1[nam == "log_sigma"]
      # previously  c(pmax(param1[nam=="Au"],rep(log(0.001), num.lv*n)), rep(0.01,num.lv*(num.lv-1)/2*n))
      # this line adds the covariance parameters after diag iter, it didn't start though, this does, sometimes.
      Au1 <- param1[nam == "Au"]

      lg_Ar1 <- param1[nam == "lg_Ar"]

      zeta <- param1[nam == "zeta"]
      # if(gamma2==0){
      lambda2 <- t(matrix(param1[nam == "lambda2"], byrow = T, ncol = num.lv, nrow = p))
      # }else{
      #   lambda2<- t(matrix(param1[nam=="lambda2"],byrow=T,ncol=num.lv,nrow=p))
      # }
      # if(gamma2==0&starting.val!="zero"){
      #   lambda3 <- rep(0,num.lv)
      # }else{
      if (gamma2 != 0) {
        lambda3 <- rep(Lambda2.start, num.lv)
      } else {
        lambda3 <- rep(0, num.lv)
      }
      # }
      # lambda3 <- ifelse(lambda3<0.1,0.5,lambda3)#added this line for the binomial

      if (starting.val == "zero") {
        # starting values for unequal tolerances, weighted averaging at the moment
        # lambda2<-matrix(lambda3,ncol=p,nrow=num.lv)#this didnt work for the NB
        # based on weighted average species SD, see canoco
        # lambda2<-matrix(1e-3,nrow=num.lv,ncol=p)
        # for(j in 1:p){
        #   for(q in 1:num.lv){
        #     lambda2[q,j]<--.5/(sum((u[,q]-(sum(y[,j]*u[,q])/sum(y[,j])))^2*y[,j])/sum(y[,j]))
        #   }
        # }
        # if(any(is.infinite(lambda2))){
        #   lambda2[is.infinite(lambda2)]<--0.5
        # }
      } else {
        mp <- list(r0 = factor(rep(NA, length(r1))), b = factor(rep(NA, length(b1))), B = factor(rep(NA, 1)), lambda = factor(rep(NA, length(lambda1))), lambda3 = factor(rep(NA, num.lv)), u = factor(rep(NA, length(u1))), lg_phi = factor(rep(NA, length(lg_phi1))), log_sigma = factor(rep(NA, length(log_sigma1))), Au = factor(rep(NA, length(Au1))), lg_Ar = factor(rep(NA, length(lg_Ar1))), zeta = factor(rep(NA, length(zeta))))
        if (row.eff == "random") {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 1, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
            map = mp, parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = log_sigma1, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
            inner.control = list(mgcmax = 1e+200, maxit = maxit),
            DLL = "qgllvm"
          )
        } else {
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 0, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
            map = mp, parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
            inner.control = list(mgcmax = 1e+200, maxit = maxit),
            DLL = "qgllvm"
          ) # GLLVM#
        }

        if (optimizer == "nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, iter.max = maxit, eval.max = maxit, trace = trace2)), silent = !trace2))
        }
        if (optimizer == "optim") {
          if (!is.null(par.scale)) {
            if (par.scale == "coef") {
              parscale <- abs(objr$par) # this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
              parscale[parscale == 0] <- 1
            } else if (is.numeric(par.scale)) {
              parscale <- rep(par.scale, length(objr$par))
            }
          } else {
            parscale <- rep(1, length(objr$par))
          }
          if (is.null(fn.scale) | !is.numeric(fn.scale)) {
            fnscale <- 1
          } else {
            fnscale <- fn.scale
          }
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr, method = "BFGS", control = list(reltol = reltol, maxit = maxit, parscale = parscale, fnscale = fnscale, trace = trace2), hessian = FALSE), silent = !trace2))
        }
        lambda2 <- t(matrix(optr$par[names(optr$par) == "lambda2"], byrow = T, ncol = num.lv, nrow = p)) # doesnt always work for NB
      }


      if (row.eff == "random") {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 1, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = log_sigma1, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        )
      } else {
        objr <- TMB::MakeADFun(
          data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, model = 0, random = 0, zetastruc = ifelse(zeta.struc == "species", 1, 0), gamma = gamma1, gamma2 = gamma2, theta4 = theta4), silent = TRUE,
          parameters = list(r0 = r1, b = b1, B = matrix(0), lambda = lambda1, lambda2 = lambda2, lambda3 = lambda3, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = lg_Ar1, zeta = zeta), # log(phi)
          inner.control = list(mgcmax = 1e+200, maxit = maxit),
          DLL = "qgllvm"
        ) # GLLVM#
      }

      if (optimizer == "nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, iter.max = maxit, eval.max = maxit, trace = trace2)), silent = !trace2))
      }
      if (optimizer == "optim") {
        if (!is.null(par.scale)) {
          if (par.scale == "coef") {
            parscale <- abs(objr$par) # this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
            parscale[parscale == 0] <- 1
          } else if (is.numeric(par.scale)) {
            parscale <- rep(par.scale, length(objr$par))
          }
        } else {
          parscale <- rep(1, length(objr$par))
        }
        if (is.null(fn.scale) | !is.numeric(fn.scale)) {
          fnscale <- 1
        } else {
          fnscale <- fn.scale
        }
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr, method = "BFGS", control = list(reltol = reltol, maxit = maxit, parscale = parscale, fnscale = fnscale, trace = trace2), hessian = FALSE), silent = !trace2))
      }
      if (optimizer == "optim") {
        if (inherits(optr, "try-error") || is.nan(optr$value) || is.na(optr$value) || is.infinite(optr$value)) {
          optr <- optr1
          objr <- objr1
          common.tolerances <- T
        }
      } else {
        if (inherits(optr, "try-error") || is.nan(optr$objective) || is.na(optr$objective) || is.infinite(optr$objective)) {
          optr <- optr1
          objr <- objr1
          common.tolerances <- T
        }
      }
    }
    return(list(objr = objr, optr = optr, fit = fit, timeo = timeo))
  }
  #
  # #find best ridge parameter values?
  # data <- objr$env$data
  # pars <- objr$env$last.par.best
  # parslist <- objr$env$parList(objr$env$last.par.best)
  # LL <- objr$env$value.best
  #
  # for(i in 2:15){
  #   data$gamma2 <- matrix(i,ncol=p,nrow=num.lv)
  #     objr <- TMB::MakeADFun(
  #       data = data, silent=TRUE,
  #       parameters = parslist,
  #       inner.control=list(mgcmax = 1e+200,maxit = maxit),
  #       DLL = "qgllvm")
  #
  #     if(optimizer=="nlminb") {
  #       timeo <- system.time(optr <- try(nlminb(pars, objr$fn, objr$gr,control = list(rel.tol=reltol, iter.max=maxit, eval.max=maxit,trace=trace2)),silent = !trace2))
  #     }
  #     if(optimizer=="optim") {
  #       if(!is.null(par.scale)){
  #         if(par.scale=="coef"){
  #           parscale<-abs(objr$par)#this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
  #           parscale[parscale==0]<-1
  #         }else if(is.numeric(par.scale)){
  #           parscale<-rep(par.scale,length(objr$par))
  #         }
  #       }else{
  #         parscale <- rep(1,length(objr$par))
  #       }
  #       if(is.null(fn.scale)|!is.numeric(fn.scale)){
  #         fnscale<-1
  #       }else{
  #         fnscale<-fn.scale
  #       }
  #       timeo <- system.time(optr <- try(optim(pars, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit,parscale=parscale,fnscale=fnscale, trace=trace2),hessian = FALSE),silent = !trace2))
  #     }
  #
  #   if(LL>objr$env$value.best){
  #     gamma2 <- data$gamma2
  #     LL <- objr$env$value.best
  #   }
  #
  # }
  #

  if (n.init > 1 & parallel == TRUE) {
    # clusterEvalQ(cl,library.dynam("qgllvm","gllvm.quadratic",lib.loc="C:/Users/beve/Documents/R/win-library/3.6/"))
    # start.values.gllvm.TMB.quadratic<-getFromNamespace("start.values.gllvm.TMB.quadratic","gllvm.quadratic")
    try(results <- foreach(i = 1:n.init, .export = ls(), .multicombine = T, .inorder = F, .packages = "gllvm") %dopar% {
      madeMod <- makeMod(i)
      # found the issue, was exporting packages. Now I need to find out how to set a seed inside a foreach..might have to set outside the function inside the forach loop due to openmp
      return(madeMod)
    }, silent = T) # to muffle export error. Need to keep ls() in to make sure all objects are exported to workers as some objects are in another environment as they're passed from the main gllvm function.
    # on.exit(stopCluster(cl))
  } else if (n.init == 1) {
    results <- makeMod(1)
  } else if (n.init > 1 & parallel == FALSE) {
    results <- vector("list", n.init)
    for (i in 1:n.init) {
      results[[i]] <- makeMod(i)
    }
  }

  if (n.init > 1) {
    try(
      {
        bestLL <- lapply(results, function(x) x$objr$fn(x$optr$par))
        objr <- results[[which.min(unlist(bestLL))]]$objr
        optr <- results[[which.min(unlist(bestLL))]]$optr
        fit <- results[[which.min(unlist(bestLL))]]$fit
        timeo <- results[[which.min(unlist(bestLL))]]$timeo
      },
      silent = T
    )
    LL <- unlist(bestLL)
  } else {
    objr <- results$objr
    optr <- results$optr
    LL <- objr$fn(optr$par)
    fit <- results$fit
    timeo <- results$timeo
  }

  if (inherits(optr, "try-error")) warning(optr[1])
  param <- optr$par
  if (family %in% c("negative.binomial", "gaussian", "gamma")) {
    phis <- exp(param[names(param) == "lg_phi"])
  }
  if (family == "ordinal") {
    K <- max(y00) - min(y00)
    zetas <- param[names(param) == "zeta"]
    if (zeta.struc == "species") {
      zetanew <- matrix(NA, nrow = p, ncol = K)
      idx <- 0
      for (j in 1:ncol(y)) {
        k <- max(y[, j]) - 2
        if (k > 0) {
          for (l in 1:k) {
            zetanew[j, l + 1] <- zetas[idx + l]
          }
        }
        idx <- idx + k
      }
      zetanew[, 1] <- 0
      row.names(zetanew) <- colnames(y00)
      colnames(zetanew) <- paste(min(y):(max(y00) - 1), "|", (min(y00) + 1):max(y00), sep = "")
    } else {
      zetanew <- c(0, zetas)
      names(zetanew) <- paste(min(y00):(max(y00) - 1), "|", (min(y00) + 1):max(y00), sep = "")
    }

    zetas <- zetanew
    out$y <- y00
  }

  bi <- names(param) == "b"
  li <- names(param) == "lambda"
  l2i <- names(param) == "lambda2"
  if (common.tolerances == F & gamma2 > 0) l3i <- names(param) == "lambda3"
  ui <- names(param) == "u"
  if (row.eff != FALSE) {
    ri <- names(param) == "r0"
    if (row.eff == "fixed") row.params <- param[ri] # c(0,param[ri])
    if (row.eff == "random") {
      sigma <- exp(param["log_sigma"])
      row.params <- param[ri]
    }
  }
  betaM <- matrix(param[bi], p, num.X + 1, byrow = TRUE)
  beta0 <- betaM[, 1]
  if (!is.null(X)) betas <- betaM[, -1]
  lvs <- matrix(param[ui], n, num.lv)
  theta <- matrix(0, p, num.lv)
  if (p > 1) {
    theta[lower.tri(theta, diag = TRUE)] <- param[li]
    if (common.tolerances == F & gamma2 > 0) theta3 <- abs(param[l3i]) + theta4
    if (common.tolerances == F & gamma2 > 0) {
      theta2 <- -(abs(matrix(param[l2i], nrow = p, ncol = num.lv, byrow = T)) + matrix(theta3, ncol = num.lv, nrow = p, byrow = T))
    } else {
      theta2 <- -matrix(abs(param[l2i]) + theta4, ncol = num.lv, nrow = p, byrow = T)
    }

    theta <- cbind(theta, theta2)
  } else {
    theta <- param[li]
    theta <- c(theta, param[l2i])
  }
  # diag(theta) <- exp(diag(theta)) !!!


  if (family %in% c("negative.binomial", "gaussian", "gamma")) {
    phis <- exp(param[names(param) == "lg_phi"])
  }

  out$start <- fit
  try(out$convergence <- optr$convergence, silent = T)
  out$logL <- objr$env$value.best[1]
  out$lvs <- lvs
  out$params$theta <- theta
  rownames(out$lvs) <- rownames(out$y)

  if (num.lv > 1) {
    colnames(out$lvs) <- paste("LV", 1:num.lv, sep = "")
    colnames(out$params$theta) <- c(paste("LV", 1:num.lv, sep = ""), paste("LV", 1:num.lv, "^2", sep = ""))
    rownames(out$params$theta) <- colnames(out$y)
  }
  if (common.tolerances == F & gamma2 != 0) out$params$theta2 <- -theta3
  names(beta0) <- colnames(out$y)
  out$params$beta0 <- beta0
  if (!is.null(X)) {
    betas <- matrix(betas, ncol = ncol(X))
    out$params$Xcoef <- betas
    rownames(out$params$Xcoef) <- colnames(out$y)
    colnames(out$params$Xcoef) <- colnames(X)
  }
  if (family == "ordinal") {
    out$params$zeta <- zetas
  }
  if (family %in% c("negative.binomial")) {
    out$params$inv.phi <- phis
    names(out$params$inv.phi) <- colnames(out$y)
    out$params$phi <- 1 / phis
    names(out$params$phi) <- colnames(out$y)
  }
  if (family %in% c("gaussian", "gamma")) {
    out$params$phi <- phis
    names(out$params$phi) <- colnames(out$y)
  }
  if (row.eff != FALSE) {
    if (row.eff == "random") {
      out$params$sigma <- sigma
      names(out$params$sigma) <- "sigma"
    }
    out$params$row.params <- row.params
    names(out$params$row.params) <- rownames(out$y)
  }
  if (family == "binomial") out$link <- link

  out$row.eff <- row.eff
  out$time <- timeo
  pars <- optr$par
  param <- objr$env$last.par.best
  Au <- param[names(param) == "Au"]
  A <- array(0, dim = c(n, num.lv, num.lv))
  for (d in 1:num.lv) {
    for (i in 1:n) {
      A[i, d, d] <- exp(Au[(d - 1) * n + i])
    }
  }
  if (length(Au) > num.lv * n) {
    k <- 0
    for (c1 in 1:num.lv) {
      r <- c1 + 1
      while (r <= num.lv) {
        for (i in 1:n) {
          A[i, r, c1] <- Au[num.lv * n + k * n + i]
          A[i, c1, r] <- A[i, r, c1]
        }
        k <- k
        r <- r + 1
      }
    }
  }
  out$A <- A

  if (row.eff == "random") {
    Ar <- exp(param[names(param) == "lg_Ar"])
    names(Ar) <- paste("Ar", 1:n, sep = "")
    out$Ar <- Ar^2
  }


  tr <- try(
    {
      if (sd.errors && !is.infinite(out$logL)) {
        if (trace) cat("Calculating standard errors for parameters...\n")
        # set parallel workers with openmp for the calculation of standard errors, if n.init>1
        # if(getDoParWorkers()>1&n.init>1&parallel==T){
        #   n.cores.old<-openmp()
        #   openmp(getDoParWorkers())
        # } #later for parallel computation of SE
        sdr <- optimHess(pars, objr$fn, objr$gr, control = list(reltol = reltol, maxit = maxit)) # maxit=maxit
        # if(getDoParWorkers()>1&n.init>1&parallel==T)openmp(n.cores.old)#reset back to old configuration
        m <- dim(sdr)[1]
        incl <- rep(TRUE, m)
        incld <- rep(FALSE, m)
        incl[names(objr$par) == "B"] <- FALSE
        incl[names(objr$par) == "lg_Ar"] <- FALSE
        incl[names(objr$par) == "Au"] <- FALSE
        incl[names(objr$par) == "r0"] <- FALSE
        incl[names(objr$par) == "log_sigma"] <- FALSE

        if (family != "ordinal") incl[names(objr$par) == "zeta"] <- FALSE

        if (common.tolerances == TRUE | gamma2 == 0) {
          incl[names(objr$par) == "lambda3"] <- FALSE
        }

        if (row.eff == "random") {
          incld[names(objr$par) == "lg_Ar"] <- TRUE
          incld[names(objr$par) == "r0"] <- FALSE
          incl[names(objr$par) == "log_sigma"] <- TRUE
        }
        if (row.eff == FALSE) {
          incl[names(objr$par) == "r0"] <- FALSE
          incl[names(objr$par) == "log_sigma"] <- FALSE
        }
        if (row.eff == "fixed") incl[names(objr$par) == "r0"] <- TRUE
        incl[1] <- FALSE
        incl[names(objr$par) == "log_sigma"] <- FALSE

        incl[names(objr$par) == "u"] <- FALSE
        incld[names(objr$par) == "u"] <- TRUE
        incld[names(objr$par) == "Au"] <- TRUE

        if (family == "binomial" || family == "ordinal" || family == "poisson") incl[names(objr$par) == "lg_phi"] <- FALSE

        if (family == "ordinal") incl[names(objr$par) == "zeta"] <- TRUE

        A.mat <- -sdr[incl, incl] # a x a
        D.mat <- -sdr[incld, incld] # d x d
        B.mat <- -sdr[incl, incld] # a x d
        cov.mat.mod <- try(MASS::ginv(A.mat - B.mat %*% solve(D.mat) %*% t(B.mat)), silent = TRUE)

        se <- sqrt(diag(abs(cov.mat.mod)))

        incla <- rep(FALSE, length(incl))
        incla[names(objr$par) == "u"] <- TRUE
        out$Hess <- list(Hess.full = sdr, incla = incla, incl = incl, incld = incld, cov.mat.mod = cov.mat.mod)

        if (row.eff == "fixed") {
          se.row.params <- c(0, se[1:(n - 1)])
          names(se.row.params) <- rownames(out$y)
          se <- se[-(1:(n - 1))]
        }
        sebetaM <- matrix(se[1:((num.X + 1) * p)], p, num.X + 1, byrow = TRUE)
        se <- se[-(1:((num.X + 1) * p))]
        se.lambdas <- matrix(0, p, num.lv)
        se.lambdas[lower.tri(se.lambdas, diag = TRUE)] <- se[1:(p * num.lv - sum(0:(num.lv - 1)))]
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep = "")
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas
        se <- se[-(1:(p * num.lv - sum(0:(num.lv - 1))))]
        # diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) !!!
        if (common.tolerances == F) {
          se.lambdas2 <- matrix(se[1:(p * num.lv)], p, num.lv, byrow = T)
        }
        if (common.tolerances == F & gamma2 > 0) {
          se.lambdas3 <- matrix(se[1:num.lv], p, num.lv, byrow = T)
          se.lambdas2 <- se.lambdas2 + se.lambdas3
          colnames(se.lambdas3) <- paste("LV", 1:num.lv, "^2", sep = "")
          se <- se[-(1:(num.lv))]
          out$sd$theta2 <- se.lambdas3[1, ]
        }
        if (common.tolerances == T) {
          se.lambdas2 <- matrix(se[1:num.lv], p, num.lv, byrow = T)
        }
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
        rownames(se.lambdas2) <- colnames(out$y)
        se <- se[-(1:(p * num.lv))]
        out$sd$theta <- cbind(out$sd$theta, se.lambdas2)


        # Calculate SE optimate
        out$sd$optima <- matrix(NA, nrow = p, ncol = num.lv)
        out$sd$tolerances <- matrix(NA, nrow = p, ncol = num.lv)

        if (gamma2 > 0 & common.tolerances == F) {
          idx <- c(
            which(colnames(sdr[incl, incl]) == "lambda"),
            which(colnames(sdr[incl, incl]) == "lambda2"),
            which(colnames(sdr[incl, incl]) == "lambda3")
          )
        } else {
          idx <- c(
            which(colnames(sdr[incl, incl]) == "lambda"),
            which(colnames(sdr[incl, incl]) == "lambda2")
          )
        }

        V <- -cov.mat.mod
        colnames(V) <- colnames(sdr[incl, incl])
        row.names(V) <- row.names(sdr[incl, incl])

        for (q in 1:num.lv) {
          for (j in 1:ncol(y)) {
            if (q > 1 & j < q) {
              # add zeros where necessary
              V <- cbind(cbind(V[, 1:c(p * q + j - 1)], 0), V[, (p * q + j):ncol(V)])
              V <- rbind(rbind(V[1:c(p * q + j - 1), ], 0), V[(p * q + j):nrow(V), ])
            }
          }
        }

        colnames(V)[colnames(V) == ""] <- "lambda"
        row.names(V)[row.names(V) == ""] <- "lambda"

        for (j in 1:p) {
          if (gamma2 > 0) {
            idx <- colnames(V) == "lambda" | colnames(V) == "lambda2" | colnames(V) == "lambda3"
            V.theta <- V[idx, idx]
            if (common.tolerances == T) {
              idx <- c((c(1:num.lv) - 1) * p + j, p * num.lv + 1:num.lv)
            } else {
              idx <- c((c(1:num.lv) - 1) * p + j, ((1 + p * num.lv + (num.lv * (j - 1))):(p * num.lv + (num.lv * (j - 1)) + num.lv)), ((1 + ncol(V.theta) - num.lv):ncol(V.theta)))
            }
          } else {
            idx <- colnames(V) == "lambda" | colnames(V) == "lambda2"
            V.theta <- V[idx, idx]
            if (common.tolerances == T) {
              idx <- c((c(1:num.lv) - 1) * p + j, p * num.lv + 1:num.lv)
            } else {
              idx <- c((c(1:num.lv) - 1) * p + j, ((1 + p * num.lv + (num.lv * (j - 1))):(p * num.lv + (num.lv * (j - 1)) + num.lv))) # The last is because the order from the tolerances is different from the lambdas1, its per species not per lv.
            }
          }

          V.theta2 <- V.theta[idx, idx]

          for (i in 1:num.lv) {
            du <- c((2 * out$params$theta[j, num.lv + i, drop = F])^-1, 2 * (out$params$theta[j, i, drop = F] / (2 * out$params$theta[j, num.lv + i, drop = F])^2))
            out$sd$optima[j, i] <- sqrt(abs(t(du) %*% V.theta2[c(i, num.lv + i), c(i, num.lv + i)] %*% du))
            # sd tolerances also
            dt <- 1 / (2 * out$params$theta[, -c(1:num.lv)][j, i] * (sqrt(-2 * out$params$theta[, -c(1:num.lv)][j, i]))) # need to be calculated with covariance of gamma3 if gamma2>0..that also requires subtracting theta3 from theta2
            out$sd$tolerances[j, i] <- sqrt(abs(V.theta2[-c(1:num.lv), -c(1:num.lv)][i, i] * dt^2))
          }
        }
        row.names(out$sd$optima) <- row.names(out$sd$tolerances) <- colnames(y)


        out$sd$theta <- out$sd$theta

        out$sd$beta0 <- sebetaM[, 1]
        names(out$sd$beta0) <- colnames(out$y)
        if (!is.null(X)) {
          out$sd$Xcoef <- matrix(sebetaM[, -1], nrow = nrow(sebetaM))
          rownames(out$sd$Xcoef) <- colnames(y)
          colnames(out$sd$Xcoef) <- colnames(X)
        }
        if (row.eff == "fixed") {
          out$sd$row.params <- se.row.params
        }

        if (family %in% c("negative.binomial")) {
          se.lphis <- se[1:p]
          out$sd$inv.phi <- se.lphis * out$params$inv.phi
          out$sd$phi <- se.lphis * out$params$phi
          names(out$sd$phi) <- colnames(y)
          se <- se[-(1:p)]
        }
        if (family %in% c("gamma", "gaussian")) {
          se.lphis <- se[1:p]
          out$sd$phi <- se.lphis * out$params$phi
          out$sd$phi <- se.lphis * out$params$phi
          names(out$sd$phi) <- colnames(y)
          se <- se[-(1:p)]
        }
        if (row.eff == "random") {
          out$sd$sigma <- se * out$params$sigma
          names(out$sd$sigma) <- "sigma"
        }

        if (family %in% c("ordinal")) {
          se.zetanew <- se.zetas <- se
          if (zeta.struc == "species") {
            se.zetanew <- matrix(NA, nrow = p, ncol = K)
            idx <- 0
            for (j in 1:ncol(y)) {
              k <- max(y[, j]) - 2
              if (k > 0) {
                for (l in 1:k) {
                  se.zetanew[j, l + 1] <- se.zetas[idx + l]
                }
              }
              idx <- idx + k
            }
            se.zetanew[, 1] <- 0
            out$sd$zeta <- se.zetanew
            row.names(out$sd$zeta) <- colnames(y00)
            colnames(out$sd$zeta) <- paste(min(y00):(max(y00) - 1), "|", (min(y00) + 1):max(y00), sep = "")
          } else {
            se.zetanew <- c(0, se.zetanew)
            out$sd$zeta <- se.zetanew
            names(out$sd$zeta) <- paste(min(y00):(max(y00) - 1), "|", (min(y00) + 1):max(y00), sep = "")
          }
        }
      }
    },
    silent = T
  )

  if (inherits(tr, "try-error")) {
    cat("Standard errors for parameters could not be calculated, due to singular fit.\n")
  }

  if (is.null(formula1)) {
    out$formula <- formula
  } else {
    out$formula <- formula1
  }

  out$TMBfn <- objr
  out$TMBfn$par <- optr$par
  out$logL <- -out$logL
  out$LL <- -LL
  out$start.struc <- start.struc
  out$common.tolerances <- common.tolerances
  out$ridge <- list(gamma1, gamma2)
  if (row.eff == "random") out$logL <- out$logL + n * 0.5
  if (family == "gaussian") {
    out$logL <- out$logL - n * p * log(pi) / 2
  }

  return(out)
}
