########################################################################################## GLLVM fourth corner model, with estimation done via Laplace and Variational approximation using TMB-package Original author:
########################################################################################## Jenni Niku, Bert van der Veen
gllvm.TMB.trait.quadratic <- function(y, X = NULL, TR = NULL, formula = NULL, num.lv = 2, family = "poisson", Lambda.struc = "unstructured", 
    row.eff = FALSE, reltol = 1e-10, seed = NULL, maxit = 2000, start.lvs = NULL, offset = NULL, sd.errors = TRUE, trace = FALSE, trace2 = FALSE, 
    n.init = 10, start.params = NULL, start0 = FALSE, optimizer = "optim", starting.val = "lingllvm", randomX = NULL, diag.iter = 1, Lambda.start = c(0.1,0.5), jitter.var = 0, yXT = NULL, ridge = FALSE, ridge.quadratic = FALSE, start.method="FA", parscale=1,fnscale=1, zeta.struc = "species", starting.val.lingllvm = "res") {
    if (is.null(X) && !is.null(TR)) 
        stop("Unable to fit a model that includes only trait covariates")
    
    term <- NULL
    n <- dim(y)[1]
    p <- dim(y)[2]
    y <- as.data.frame(y)
    formula1 <- formula
    if (family == "binomial") {
        link <- "probit"
    }
    if (num.lv == 0) {
        stop("Can't fit the species packing model without latent variables")
    }
    if(family == "ordinal") {
      y00<-y
      if(min(y)==0){ y=y+1}
      max.levels <- apply(y,2,function(x) length(min(x):max(x)))
      if(any(max.levels == 1) || all(max.levels == 2))
        stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
                          
      if(any(!apply(y,2,function(x)all(diff(sort(unique(x)))==1)))&zeta.struc=="species")
        stop("Can't fit ordinal model if there are species with missing classes. Please reclassify per species or use zeta.struc = `common` ")
      
      if(any(diff(sort(unique(c(y))))!=1)&zeta.struc=="common")
        stop("Can't fit ordinal model if there are missing classes. Please reclassify.")
    }
    if (NCOL(X) < 1) 
        stop("No covariates in the model, fit the model using gllvm(y,family=", family, "...)")
    
    # change categorical variables to dummy variables
    num.X <- 0
    X.new <- NULL
    if (!is.null(X)) {
        num.X <- dim(X)[2]
        for (i in 1:num.X) {
            if (!is.factor(X[, i])) {
                if (length(unique(X[, i])) > 2) {
                  Xi <- scale(X[, i])
                } else {
                  Xi <- X[, i]
                }
                X[, i] <- Xi
                X.new <- cbind(X.new, Xi)
                if (!is.null(colnames(X)[i])) 
                  colnames(X.new)[dim(X.new)[2]] <- colnames(X)[i]
            } else {
                dum <- model.matrix(~X[, i])
                dum <- dum[, !(colnames(dum) %in% c("(Intercept)"))]
                colnames(dum) <- paste(colnames(X)[i], levels(X[, i])[-1], sep = "")
                X.new <- cbind(X.new, dum)
            }
        }
        X.new <- data.frame(X.new)
    }
    if (!is.null(randomX)) {
        method <- "LA"
        xb <- as.matrix(model.matrix(randomX, data = X.new))
        xb <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
        Br <- matrix(0, ncol(xb), p)
        sigmaB <- diag(ncol(xb))
    }
    
    num.T <- 0
    T.new <- NULL
    if (!is.null(TR)) {
        num.T <- dim(TR)[2]
        T.new <- matrix(0, p, 0)
        if (num.T > 0) {
            for (i in 1:num.T) {
                if (!is.factor(TR[, i]) && length(unique(TR[, i])) > 2) {
                  TR[, i] <- scale(TR[, i])
                  T.new <- cbind(T.new, scale(TR[, i]))
                  colnames(T.new)[dim(T.new)[2]] <- colnames(TR)[i]
                } else {
                  dum <- model.matrix(~TR[, i] - 1)
                  colnames(dum) <- paste(colnames(TR)[i], levels(TR[, i]), sep = "")
                  T.new <- cbind(T.new, dum)
                }
            }
            T.new <- data.matrix(T.new)
        }
    }
    
    
    if (is.null(formula)) {
        n1 <- colnames(X)
        n2 <- colnames(TR)
        form1 <- paste("", n1[1], sep = "")
        if (length(n1) > 1) {
            for (i1 in 2:length(n1)) {
                form1 <- paste(form1, n1[i1], sep = "+")
            }
        }
        formula <- paste("y~", form1, sep = "")
        formula <- paste(formula, form1, sep = " + (")
        
        formula <- paste(formula, ") : (", sep = "")
        formula <- paste(formula, n2[1], sep = "")
        if (length(n2) > 1) {
            for (i2 in 2:length(n2)) {
                formula <- paste(formula, n2[i2], sep = "+")
            }
        }
        formula1 <- paste(formula, ")", sep = "")
        formula <- formula(formula1)
    }
    
    if (!is.null(X) || !is.null(TR)) {
        yX <- reshape(data.frame(cbind(y, X)), direction = "long", varying = colnames(y), v.names = "y")
        TR2 <- data.frame(time = 1:p, TR)
        if (is.null(yXT)) {
            yXT <- merge(yX, TR2, by = "time")
        }
        data <- yXT
        
        m1 <- model.frame(formula, data = data)
        term <- terms(m1)
        
        Xd <- as.matrix(model.matrix(formula, data = data))
        nXd <- colnames(Xd)
        Xd <- as.matrix(Xd[, !(nXd %in% c("(Intercept)"))])
        colnames(Xd) <- nXd[!(nXd %in% c("(Intercept)"))]
        if (!is.null(X.new)) 
            fx <- apply(matrix(sapply(colnames(X.new), function(x) {
                grepl(x, colnames(Xd))
            }), ncol(Xd), ncol(X.new)), 2, any)
        ft <- NULL
        if (NCOL(T.new) > 0) {
            ft <- apply(matrix(sapply(colnames(T.new), function(x) {
                grepl(x, colnames(Xd))
            }), ncol(Xd), ncol(T.new)), 2, any)
        }
        X1 <- as.matrix(X.new[, fx])
        TR1 <- as.matrix(T.new[, ft])
        colnames(X1) <- colnames(X.new)[fx]
        colnames(TR1) <- colnames(T.new)[ft]
        nxd <- colnames(Xd)
        formulab <- paste("~", nxd[1], sep = "")
        for (i in 2:length(nxd)) formulab <- paste(formulab, nxd[i], sep = "+")
        formula1 <- formulab
    }
    
    
    if (!(family %in% c("poisson", "negative.binomial", "binomial", "ordinal"))) 
        stop("Selected family not permitted...sorry!")
    if (!(Lambda.struc %in% c("unstructured", "diagonal"))) 
        stop("Lambda matrix (covariance of vartiational distribution for latent variable) not permitted...sorry!")
    if (num.lv == 1) 
        Lambda.struc <- "diagonal"  ## Prevents it going to 'unstructured' loops and causing chaos
    trial.size <- 1
    
    y <- as.matrix(y)
    if (!is.numeric(y)) 
        stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
    
    if (is.null(rownames(y))) 
        rownames(y) <- paste("Row", 1:n, sep = "")
    if (is.null(colnames(y))) 
        colnames(y) <- paste("Col", 1:p, sep = "")
    if (!is.null(X)) {
        if (is.null(colnames(X))) 
            colnames(X) <- paste("x", 1:ncol(X), sep = "")
    }
    
    out <- list(y = y, X = X1, TR = TR1, num.lv = num.lv, row.eff = row.eff, logL = Inf, family = family, offset = offset, randomX = randomX, 
        X.design = Xd, terms = term)
    old.logL <- Inf
    if (is.null(formula) && is.null(X) && is.null(TR)) {
        formula = "~ 1"
    }
    
    n.i <- 1
    if (n.init > 1) 
        seed <- sample(1:10000, n.init)
    if(starting.val=="lingllvm"){
      n.init2<-n.init
      n.init<-1
      #check if I've covered all options
      fit <- gllvm(y, formula = formula, X = X, TR = TR, num.lv = num.lv, family = family, row.eff = row.eff, n.init = n.init2, maxit = maxit, reltol=reltol, optimizer = optimizer, start.fit = start.params, diag.iter = diag.iter, jitter.var = jitter.var, starting.val = starting.val.lingllvm, Lambda.start = Lambda.start, seed = seed, Lambda.struc = Lambda.struc)
    }
    
    while (n.i <= n.init) {
        
        num.X <- dim(X)[2]
        num.T <- dim(TR)[2]
        
        if (n.init > 1 && trace) 
          if(n.init > 1 && trace)
            if(n.i==2|old.logL>out$logL){
              cat("Initial run ", n.i, "LL",out$logL , "\n")
            }else{
              cat("Initial run ", n.i, "\n")
            }
        old.logL <- out$logL
        if(starting.val!="lingllvm"){
        res <- start.values.gllvm.TMB.quadratic(y = y, X = X1, TR = TR1, family = family, offset = offset, trial.size = trial.size, 
            num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i], starting.val = starting.val, formula = formula, jitter.var = jitter.var, 
            yXT = yXT, row.eff = row.eff, start.method=start.method, zeta.struc = zeta.struc)
        }else{
          start.params <- fit
        }
        if (is.null(start.params)) {
            beta0 <- res$params[, 1]
            # common env params or different env response for each spp
            B <- NULL
            if (!is.null(TR) && !is.null(X)) {
                B <- c(res$B)[1:ncol(Xd)]
                if (any(is.na(B))) 
                  B[is.na(B)] <- 0
            }
            row.params <- NULL
            if (row.eff != FALSE) {
                row.params <- res$row.params
                if (row.eff == "random") {
                  sigma <- sd(row.params)
                }
            }
            vameans <- theta <- lambda <- lambda2 <- NULL
            
            
            vameans <- res$index
            theta <- matrix(res$params[, (ncol(res$params) - num.lv * 2 + 1):ncol(res$params)], ncol = num.lv * 2)  #fts$coef$theta#
            theta[, 1:num.lv][upper.tri(theta[, 1:num.lv])] <- 0
            if (Lambda.struc == "unstructured") {
                lambda <- array(NA, dim = c(n, num.lv, num.lv))
                for (i in 1:n) {
                  lambda[i, , ] <- diag(rep(1, num.lv))
                }
            }
            if (Lambda.struc == "diagonal") {
                lambda <- matrix(1, n, num.lv)
            }
            zero.cons <- which(theta == 0)
            
        } else {
            if (dim(start.params$y) == dim(y) && is.null(X) == is.null(start.params$X) && is.null(T) == is.null(start.params$TR) && 
                row.eff == start.params$row.eff) {
              if(start.params$family=="ordinal"){
                  zeta <- start.params$zeta 
              }
                beta0 <- start.params$params$beta0
                # common env params or different env response for each spp
                B <- NULL
                if (!is.null(TR) && !is.null(X)) {
                  B <- start.params$params$B
                }
                fourth <- inter <- NULL
                if (!is.null(TR)) 
                  inter <- start.params$params$fourth  # let's treat this as a vector (vec(B'))'
                vameans <- theta <- lambda <- NULL
                
                if (row.eff) 
                  row.params <- start.params$params$row.params  ## row parameters
                theta <- c(start.params$params$theta)[,1:num.lv]  ## LV coefficients
                if(class(start.params)=="gllvm.quadratic"){
                  theta2 <- start.params$params$theta[,-c(1:num.lv)] 
                }else{
                  if(starting.val!="lingllvm"){
                    if(!is.null(X)){
                      theta2 <- fit$params[,(ncol(fit$params)-num.lv+1):ncol(fit$params)]#need to double check this
                    }else if(is.null(X)){
                      theta2 <- fit$params[,-c(1:(num.lv+1))]  
                    }#this still needs to be implement for traits.
                  }else{
                    theta2<-matrix(0,nrow=p,ncol=num.lv)
                    for(j in 1:p){
                      for(q in 1:num.lv){
                        theta2[j,q]<--.5/(sum((start.params$lvs[,q]-(sum(y[,j]*start.params$lvs[,q])/sum(y[,j])))^2*y[,j])/sum(y[,j]))
                      }
                    }
                    if(any(is.infinite(theta2))){
                      theta2[is.infinite(theta2)]<--0.5
                    }
                  }
                }
                
                vameans <- matrix(start.params$lvs, ncol = num.lv)
                lambda <- start.params$A
            } else {
                stop("Model which is set as starting parameters isn't the suitable you are trying to fit. Check that attributes y, X, TR and row.eff match to each other.")
            }
        }
        if (is.null(offset)) 
            offset <- matrix(0, nrow = n, ncol = p)
        
        phi <- phis <- NULL
        if (family == "negative.binomial") {
            phis <- res$phi
            if (any(phis > 10)) 
                phis[phis > 100] <- 100
            if (any(phis < 0.01)) 
                phis[phis < 0.01] <- 0.01
            res$phi <- phis
            phis <- 1/phis
            
        }
        if(family=="ordinal"){
          K = max(y00)-min(y00)
          if(zeta.struc=="species"){
            zeta <- c(t(fit$zeta[,-1]))
            zeta <- zeta[!is.na(zeta)]
          }else{
            zeta <- fit$zeta[-1]
          }
          
        }else{
          zeta = 0
        }
        
        q <- num.lv
        
        if (!is.null(row.params)) {
            r0 <- row.params
        } else {
            r0 <- rep(0, n)
        }
        a <- c(beta0)
        # diag(theta) <- log(diag(theta)) !!!
        theta2 <- theta[, -c(1:num.lv)]
        theta <- theta[, 1:num.lv][lower.tri(theta[, 1:num.lv], diag = TRUE)]
        
        u <- vameans
        if (!is.null(phis)) {
            phi = (phis)
        } else {
            phi <- rep(1, p)
        }
        q <- num.lv
        sigma <- 1
        
        if (is.null(start.params) || start.params$method != "VA") {
            if (Lambda.struc == "diagonal" || diag.iter > 0) {
                Au <- log(rep(Lambda.start[1], num.lv * n))  #
            } else {
                Au <- c(log(rep(Lambda.start[1], num.lv * n)), rep(0, num.lv * (num.lv - 1)/2 * n))  #1/2, 1
            }
        } else {
            Au = NULL
            for (d in 1:num.lv) {
                if (start.params$Lambda.struc == "unstructured" || length(dim(start.params$A)) == 3) {
                  Au = c(Au, log(start.params$A[, d, d]))
                } else {
                  Au <- c(Au, log(start.params$A[, d]))
                }
            }
            if(Lambda.struc!="diagonal" && diag.iter==0){
              if(start.params$Lambda.struc=="unstructured"){
                for(d1 in 1:num.lv) {
                  for(d2 in 1:num.lv) {
                    if(d1!=d2){
                      Au <- c(Au,start.params$A[,d1,d2] )
                    }
                  }
                }
              }else{
                Au <- c(Au,rep(0,num.lv*(num.lv-1)/2*n))##this is redundant
              }
              
            }
        }
        if (length(Lambda.start) < 2) {
            Ar <- rep(1, n)
        } else {
            Ar <- rep(Lambda.start[2], n)
        }
        optr <- NULL
        timeo <- NULL
        se <- NULL
        
        
        if (row.eff == FALSE) {
            xr <- matrix(0, 1, p)
        } else {
            xr <- matrix(1, 1, p)
        }
        if (family == "poisson") {
            familyn = 0
        }
        if (family == "negative.binomial") {
            familyn = 1
        }
        if (family == "binomial") {
            familyn <- 2
        }
        if (family == "ordinal") {
            familyn = 3
        }
        
        if (row.eff == "random") {
            # || !is.null(randomX)
            
            if (ridge == T) {
                if (ridge.quadraitc == F) {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 1, ridge = 1, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, 
                    lg_Ar = log(Ar), zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, 
                    maxit = maxit), DLL = "gllvm2")
                } else {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 1, ridge = 1, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, 
                    lg_Ar = log(Ar), zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, 
                    maxit = maxit), DLL = "gllvm2")
                }
            } else {
                if (ridge.quadratic == F) {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 1, ridge = 0, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, 
                    lg_Ar = log(Ar), zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, 
                    maxit = maxit), DLL = "gllvm2")
                } else {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 1, ridge = 0, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = log(sigma), Au = Au, 
                    lg_Ar = log(Ar), zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, 
                    maxit = maxit), DLL = "gllvm2")
                }
            }
            
            
        } else {
            if (ridge == T) {
                if (ridge.quadratic == F) {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 0, ridge = 1, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), 
                    zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                    DLL = "gllvm2")
                } else {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 0, ridge = 1, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), 
                    zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                    DLL = "gllvm2")
                  
                }
            } else {
                if (ridge.quadratic == F) {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 0, ridge = 0, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), 
                    zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                    DLL = "gllvm2")
                } else {
                  objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, method = 0, 
                    model = 1, random = 0, ridge = 0, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = matrix(r0), b = rbind(a), 
                    B = matrix(B), lambda = theta, lambda2 = t(theta2), u = u, lg_phi = log(phi), log_sigma = 0, Au = Au, lg_Ar = log(Ar), 
                    zeta = zeta, lg_gamma = rep(0, num.lv), lg_gamma2 = rep(0, num.lv)), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                    DLL = "gllvm2")
                }
            }
        }
        
        if (optimizer == "nlminb") {
            timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, maxit = maxit, trace = trace2)), 
                silent = !trace2))
        }
        if (optimizer == "optim") {
          if(!is.null(par.scale)){
            if(par.scale=="coef"){
              parscale<-abs(objr$par)
              parscale[parscale==0]<-1
            }else if(is.numeric(par.scale)){
              parscale<-rep(par.scale,length(objr$par))
            }
          }else{
            parscale <- rep(1,length(objr$par))
          }
          if(is.null(fns.cale)|!is.numeric(fn.scale)){
            fnscale<-1
          }
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit,parscale=parscale, fnscale=fnscale, trace = trace2),hessian = FALSE),silent = !trace2))
        }
        if (inherits(optr, "try-error")) 
            warning(optr[1])
        
        if (diag.iter > 0 && Lambda.struc == "unstructured" && row.eff == "random" && is.null(randomX)) {
            objr1 <- objr
            optr1 <- optr
            param1 <- optr$par
            nam <- names(param1)
            r1 <- matrix(param1[nam == "r0"])
            b1 <- rbind(param1[nam == "b"])
            B1 <- matrix(param1[nam == "B"])
            
            lambda1 <- param1[nam == "lambda"]
            lambda2 <- matrix(-1*abs(param1[nam=="lambda2"]),nrow=num.lv)
            u1 <- matrix(param1[nam == "u"], n, num.lv)
            lg_phi1 <- param1[nam == "lg_phi"]
            lg_sigma1 <- param1[nam == "log_sigma"]
            #previously  c(pmax(param1[nam=="Au"],rep(log(0.001), num.lv*n)), rep(0.01,num.lv*(num.lv-1)/2*n))
            #this line adds the covariance parameters after diag iter, it didn't start though, this does.
            Au1<- c(rep(0,length(param1[names(param1)=="Au"])), rep(0,num.lv*(num.lv-1)/2*n))
            Ar1 <- param1[nam == "lg_Ar"]
            lg_gamma <- param1[nam == "lg_gamma"]
            lg_gamma2 <- param1[nam == "lg_gamma2"]
            zeta <- param1[nam == "zeta"]
            
            if (row.eff == "random") {
                
                if (ridge == T) {
                  if (ridge.quadratic == F) {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 1, ridge = 1, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = lg_sigma1, Au = Au1, 
                      lg_Ar = Ar1, zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, 
                      maxit = 1000), DLL = "gllvm2")
                  } else {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 1, ridge = 1, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = lg_sigma1, Au = Au1, 
                      lg_Ar = Ar1, zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, 
                      maxit = 1000), DLL = "gllvm2")
                  }
                } else {
                  if (ridge.quadratic == F) {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 1, ridge = 0, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = lg_sigma1, Au = Au1, 
                      lg_Ar = Ar1, zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, 
                      maxit = 1000), DLL = "gllvm2")
                  } else {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 1, ridge = 0, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = lg_sigma1, Au = Au1, 
                      lg_Ar = Ar1, zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, 
                      maxit = 1000), DLL = "gllvm2")
                  }
                }
                
            } else {
                if (ridge == T) {
                  if (ridge.quadratic == F) {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 0, ridge = 1, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = Ar1, 
                      zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                      DLL = "gllvm2")
                  } else {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 0, ridge = 1, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = Ar1, 
                      zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                      DLL = "gllvm2")
                  }
                } else {
                  if (ridge.quadratic == F) {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 0, ridge = 0, ridge_quadratic = 0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = Ar1, 
                      zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                      DLL = "gllvm2")
                  } else {
                    objr <- TMB::MakeADFun(data = list(y = y, x = Xd, xr = xr, offset = offset, num_lv = num.lv, family = familyn, 
                      method = 0, model = 1, random = 0, ridge = 0, ridge_quadratic = 1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent = TRUE, parameters = list(r0 = r1, 
                      b = b1, B = B1, lambda = lambda1, lambda2 = lambda2, u = u1, lg_phi = lg_phi1, log_sigma = 0, Au = Au1, lg_Ar = Ar1, 
                      zeta = zeta, lg_gamma = lg_gamma, lg_gamma2 = lg_gamma2), inner.control = list(mgcmax = 1e+200, maxit = 1000), 
                      DLL = "gllvm2")
                  }
                }
            }
            
            if (optimizer == "nlminb") {
                timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol, trace = trace2)), silent = !trace2))
            }
            if (optimizer == "optim") {
              if(!is.null(par.scale)){
                if(par.scale=="coef"){
                  parscale<-abs(objr$par)
                  parscale[parscale==0]<-1
                }else if(is.numeric(par.scale)){
                  parscale<-rep(par.scale,length(objr$par))
                }
              }else{
                parscale <- rep(1,length(objr$par))
              }
              timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit,parscale=parscale, fnscale=fnscale, trace = trac2),hessian = FALSE),silent = !trace2))
            }
            if (inherits(optr, "try-error")) {
                optr <- optr1
                objr <- objr1
                Lambda.struc <- "diagonal"
            }
            
        }
        if(inherits(optr,"try-error")) warning(optr[1]);
        param <- objr$env$last.par.best
        if (family %in% c("negative.binomial")) {
            phis = exp(param[names(param) == "lg_phi"])
        }
        if(family == "ordinal"){
          zetas <- param[names(param)=="zeta"]
          if(zeta.struc=="species"){
            zetanew <- matrix(NA,nrow=p,ncol=K)
            idx<-0
            for(j in 1:ncol(y)){
              k<-max(y[,j])-2
              if(k>0){
                for(l in 1:k){
                  zetanew[j,l+1]<-zetas[idx+l]
                } 
              }
              idx<-idx+k
            }
            zetanew[,1] <- 0 
            row.names(zetanew) <- colnames(y00); colnames(zetanew) <- paste(min(y):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
          }else{
            zetanew <- c(0,zetas)
            names(zetanew) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
          }
          
          zetas<-zetanew
          out$y<-y00
        }
        if (ridge == T) {
            gi <- names(param) == "lg_gamma"
        }
        if (ridge.quadratic == T) {
            gi2 <- names(param) == "lg_gamma2"
        }
        bi <- names(param) == "b"
        Bi <- names(param) == "B"
        li <- names(param) == "lambda"
        l2i <- names(param) == "lambda2"
        ui <- names(param) == "u"
        if (row.eff != FALSE) {
            ri <- names(param) == "r0"
            if (row.eff == "random") {
                row.params = param[ri]
            } else {
                row.params <- param[ri]
            }
            if (row.eff == "random") 
                sigma <- exp(param["log_sigma"])
        }
        if (!is.null(randomX)) {
            Bri <- names(param) == "Br"
            Br <- matrix(param[Bri], ncol(xb), p)
            Sri <- names(param) == "sigmaB"
            L <- diag(ncol(xb))
            if (ncol(xb) > 1) {
                sigmaB <- diag(exp(param[Sri]))
                Srij <- names(param) == "sigmaij"
                Sr <- param[Srij]
                L[upper.tri(L)] <- Sr
                D <- diag(diag(t(L) %*% L))
            } else {
                D <- 1
                sigmaB <- (exp(param[Sri]))
            }
            sigmaB_ <- solve(sqrt(D)) %*% (t(L) %*% L) %*% solve(sqrt(D))
            sigmaB <- sigmaB %*% sigmaB_ %*% t(sigmaB)
            
        }
        beta0 <- param[bi]
        B <- param[Bi]
        
        lvs <- (matrix(param[ui], n, q))
        theta <- matrix(0, p, num.lv)
        if (p > 1) {
            theta[lower.tri(theta, diag = TRUE)] <- param[li]
            theta <- cbind(theta, -1 * abs(matrix(param[l2i], ncol = num.lv)))
        } else {
            theta <- param[li]
            theta <- c(theta, param[l2i])
        }
        # diag(theta) <- exp(diag(theta)) !!!
        
        new.loglik <- objr$env$value.best[1]
        
        
        if ((n.i == 1 || out$logL > (new.loglik)) && is.finite(new.loglik) && !inherits(optr, "try-error")) {
            objr1 <- objr
            optr1 <- optr
            out$convergence <- optr1$convergence
            out$logL <- new.loglik
            
            out$lvs <- lvs
            out$params$theta <- theta
            rownames(out$lvs) <- rownames(out$y)
            if (ridge == T) {
                out$params$gamma <- exp(param[gi])
                names(out$params$gamma) <- paste("LV", 1:num.lv, sep = "")
            }
            if (ridge.quadratic == T) {
                out$params$gamma2 <- exp(param[gi2])
                names(out$params$gamma2) <- paste("LV", 1:num.lv, "^2", sep = "")
            }
            colnames(out$params$theta) <- c(paste("LV", 1:num.lv, sep = ""), paste("LV", 1:num.lv, "^2", sep = ""))
            rownames(out$params$theta) <- colnames(out$y)
            
            names(beta0) <- colnames(out$y)
            out$params$beta0 <- beta0
            out$params$B <- B
            names(out$params$B) = colnames(Xd)
            
            if (row.eff != FALSE) {
                if (row.eff == "random") {
                  out$params$sigma = sigma
                  names(out$params$sigma) = "sigma"
                }
                out$params$row.params <- row.params
                names(out$params$row.params) <- rownames(out$y)
            }
            if (family %in% c("negative.binomial")) {
                out$params$phi <- 1/phis
                names(out$params$phi) <- colnames(out$y)
                out$params$inv.phi <- phis
                names(out$params$inv.phi) <- colnames(out$y)
            }
            if (family == "ordinal") {
                out$params$zeta <- zetas
            }
            if (!is.null(randomX)) {
                out$params$Br <- Br
                out$params$sigmaB <- sigmaB
                out$corr <- sigmaB_
            }
            if (family == "binomial") 
                out$link <- link
            out$row.eff <- row.eff
            out$time <- timeo
            out$start <- res
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
                    k <- k + 1
                    r <- r + 1
                  }
                }
            }
            out$A <- A
            
            if (row.eff == "random") {
                Ar <- exp(param[names(param) == "lg_Ar"])
                out$Ar <- (Ar)^2
            }
        }
        
        
        
        n.i <- n.i + 1
    }
    tr <- try({
        if (sd.errors && !is.infinite(out$logL)) {
            if (trace) 
                cat("Calculating standard errors for parameters...\n")
            
            sdr <- optimHess(pars, objr$fn, objr$gr, control = list(reltol = reltol, maxit = maxit))
            
            m <- dim(sdr)[1]
            incl <- rep(TRUE, m)
            incld <- rep(FALSE, m)
            incl[names(objr$par) == "lg_Ar"] <- FALSE
            if (row.eff != "random") {
                incl[names(objr$par) == "log_sigma"] <- FALSE
            }
            if (familyn != 6) 
                incl[names(objr$par) == "zeta"] <- FALSE
            
            incl[names(objr$par) == "Au"] <- FALSE
            incld[names(objr$par) == "Au"] <- TRUE
            
            incl[names(objr$par)=="lg_gamma"] <- FALSE;
            incl[names(objr$par)=="lg_gamma2"] <- FALSE;
            
            if (row.eff == "random") {
                incl[names(objr$par) == "lg_Ar"] <- FALSE
                incld[names(objr$par) == "lg_Ar"] <- TRUE
                incl[names(objr$par) == "r0"] <- FALSE
                incld[names(objr$par) == "r0"] <- TRUE
            }
            if (row.eff == FALSE) 
                incl[names(objr$par) == "r0"] <- FALSE
            if (row.eff == "fixed") 
                incl[1] <- FALSE
            incl[names(objr$par) == "u"] <- FALSE
            incld[names(objr$par) == "u"] <- TRUE
            if (familyn == 0 || familyn == 2 || familyn == 6) 
                incl[names(objr$par) == "lg_phi"] <- FALSE
            if (familyn == 6) 
                incl[names(objr$par) == "zeta"] <- TRUE
            
            
            A.mat <- -sdr[incl, incl]  # a x a
            D.mat <- -sdr[incld, incld]  # d x d
            B.mat <- -sdr[incl, incld]  # a x d
            cov.mat.mod <- try(MASS::ginv(A.mat - B.mat %*% solve(D.mat) %*% t(B.mat)),silent=T)
            se <- sqrt(diag(abs(cov.mat.mod)))
            
            incla<-rep(FALSE, length(incl))
            incla[names(objr$par)=="u"] <- TRUE
            out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
            
            
            if (row.eff == "fixed") {
                se.row.params <- c(0, se[1:(n - 1)])
                names(se.row.params) = rownames(out$y)
                se <- se[-(1:(n - 1))]
            }
            se.beta0 <- se[1:p]
            se <- se[-(1:p)]
            se.B <- se[1:length(B)]
            se <- se[-(1:length(B))]
            
            se.theta <- matrix(0, p, num.lv)
            se.theta[lower.tri(se.theta, diag = TRUE)] <- se[1:(p * num.lv - sum(0:(num.lv - 1)))]
            colnames(se.theta) <- paste("LV", 1:num.lv, sep = "")
            rownames(se.theta) <- colnames(out$y)
            out$sd$theta <- se.theta
            se <- se[-(1:(p * num.lv - sum(0:(num.lv - 1))))]
            # diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) !!!
            se.lambdas2 <- matrix(se[1:(p * num.lv)], p, num.lv,byrow=T)
            colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2", sep = "")
            rownames(se.lambdas2) <- colnames(out$y)
            out$sd$theta <- cbind(out$sd$theta, se.lambdas2)
            se <- se[-(1:(p * num.lv))]
            
            out$sd$beta0 <- se.beta0
            names(out$sd$beta0) <- colnames(out$y)
            out$sd$B <- se.B
            names(out$sd$B) <- colnames(Xd)
            if (row.eff == "fixed") {
                out$sd$row.params <- se.row.params
            }
            
            if (family %in% c("negative.binomial")) {
                se.lphis <- se[1:p]
                out$sd$inv.phi <- se.lphis * out$params$inv.phi
                out$sd$phi <- se.lphis * out$params$phi
                names(out$sd$inv.phi) <- names(out$sd$phi) <- colnames(y)
                se <- se[-(1:p)]
            }
            if (row.eff == "random") {
                out$sd$sigma <- se[1] * out$params$sigma
                names(out$sd$sigma) <- "sigma"
                se = se[-1]
            }
            if (!is.null(randomX)) {
                nr <- ncol(xb)
                out$sd$sigmaB <- se * c(diag(out$params$sigmaB), rep(1, nr * (nr - 1)/2))
            }
            if(family %in% c("ordinal")){
              se.zetanew <- se.zetas <- se;
              if(zeta.struc == "species"){
                se.zetanew <- matrix(NA,nrow=p,ncol=K)
                idx<-0
                for(j in 1:ncol(y)){
                  k<-max(y[,j])-2
                  if(k>0){
                    for(l in 1:k){
                      se.zetanew[j,l+1]<-se.zetas[idx+l]
                    } 
                  }
                  idx<-idx+k
                }
                se.zetanew[,1] <- 0
                out$sd$zeta <- se.zetanew
                row.names(out$sd$zeta) <- colnames(y00); colnames(out$sd$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
                
              }else{
                se.zetanew <- c(0, se.zetanew)
                out$sd$zeta <- se.zetanew
                names(out$sd$zeta) <- paste(min(y00):(max(y00)-1),"|",(min(y00)+1):max(y00),sep="")
                
              }
            }
            
        }}, silent=T)
    
    if (inherits(tr, "try-error")) {
        cat("Standard errors for parameters could not be calculated, due to singular fit.\n")
    }
    
    if (is.null(formula1)) {
        out$formula <- formula
    } else {
        out$formula <- formula1
    }
    
    out$D <- Xd
    
    out$TMBfn <- objr1
    out$TMBfn$par <- optr1$par  #ensure params in this fn take final values
    out$logL <- -out$logL
    
    # if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
    if (row.eff == "random") 
        out$logL = out$logL + n * 0.5
    # if(!is.null(randomX)) out$logL = out$logL + p*0.5*ncol(xb)
    
    return(out)
}

