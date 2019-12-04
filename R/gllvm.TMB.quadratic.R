########################################################################################
## GLLVM, with estimation done via Variational approximation using TMB-package
## Original author: Jenni Niku, Bert van der Veen
##########################################################################################
gllvm.TMB.quadratic <- function(y, X = NULL, formula = NULL, num.lv = 2, family = "poisson",
                                Lambda.struc="unstructured", row.eff = FALSE, reltol = reltol, trace = trace, trace2 = trace2,
                                seed = NULL,maxit = maxit, start.lvs = NULL, offset=NULL, sd.errors = TRUE,
                                n.init=n.init,start.params=NULL,
                                optimizer="optim",starting.val="res",diag.iter=diag.iter,
                                Lambda.start=Lambda.start, jitter.var=jitter.var, ridge=ridge, ridge.quadratic = ridge.quadratic, start.method=start.method, par.scale=par.scale, fn.scale=fn.scale, zeta.struc = zeta.struc) {
  n <- dim(y)[1]
  p <- dim(y)[2]
  tr <- NULL
  num.lv <- num.lv
  y <- as.matrix(y)
  formula1 <- formula
  
  if (!is.numeric(y))
    stop( "y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  if (is.null(rownames(y)))
    rownames(y) <- paste("Row", 1:n, sep = "")
  if (is.null(colnames(y)))
    colnames(y) <- paste("Col", 1:p, sep = "")
  if(family=="binomial"){
    link="probit"
  }
  if(num.lv==0){
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
  num.X <- 0;
  
  if(!is.null(X)){
    
    if (!is.null(formula)) {
      xb <- as.matrix(model.matrix(formula, data = data.frame(X)))
      X <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)"))])
      colnames(X) <- colnames(xb)[!(colnames(xb) %in% c("(Intercept)"))]
      Xd <- X1 <- X
      
      num.X <- dim(X)[2]
    } else {
      n1 <- colnames(X)
      formula = paste("~", n1[1], sep = "")
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
      
      for (i in 2:length(nxd))
        formulab <- paste(formulab, nxd[i], sep = "+")
      formula1 <- formulab
    }}
  if (is.null(formula) && is.null(X)) {
    formula = "~ 1"
  }
  
  ## Set initial values for model parameters (including dispersion prm) and latent variables
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  n.i <- 1
  
  out <- list( y = y, X = X, logL = Inf, X.design = X)
  old.logL <- Inf
  if (n.init > 1)
    seed <- sample(1:10000, n.init)
  
  while(n.i <= n.init){
    if(n.init > 1 && trace)
      if(n.i==2|old.logL>out$logL){
        cat("Initial run ", n.i, "LL",out$logL , "\n")
      }else{
        cat("Initial run ", n.i, "\n")
      }
    old.logL <- out$logL

    fit <- start.values.gllvm.TMB.quadratic(y = y, X = X, TR = NULL, family = family, offset= offset, num.lv = num.lv, start.lvs = start.lvs, seed = seed[n.i], starting.val = starting.val, jitter.var = jitter.var, row.eff = row.eff, start.method=start.method, zeta.struc = zeta.struc)
    
    sigma <- 1
    if (is.null(start.params)) {
      beta0 <- fit$params[, 1]
      betas <- NULL
      if (!is.null(X))
        betas <- c(fit$params[, 2:(num.X + 1)])
      lambdas <- NULL
      
        lambdas <- as.matrix(fit$params[,(ncol(fit$params)-num.lv*2+1):(ncol(fit$params)-num.lv)])
          lambdas[upper.tri(lambdas)] <- 0
      
      row.params <- NULL
      
      if (row.eff != FALSE) {
        row.params <- fit$row.params
        if (row.eff == "random") {
          sigma <- sd(row.params)#1;#
        }
      }#rep(0,n)
      lvs <- NULL
        lvs <- matrix(fit$index, ncol = num.lv)
        if(!is.null(X)){
          lambdas <- fit$params[,(ncol(fit$params) - num.lv*2 + 1):(ncol(fit$params)-num.lv)]
          lambda2 <- fit$params[,(ncol(fit$params)-num.lv+1):ncol(fit$params)]
        }else if(is.null(X)){
          lambdas <- fit$params[,(ncol(fit$params) - num.lv*2 + 1):(ncol(fit$params)-num.lv)]
          lambda2 <- fit$params[,-c(1:(num.lv+1))]  
        }
        
      
      #subtract a fraction from the 0 quadratic scores, otherwise the optimization can't get away from the 0s where necessary.
    } else{
      if (dim(start.params$y) == dim(y) &&
          is.null(X) == is.null(start.params$X) &&
          (row.eff == start.params$row.eff)) {
        beta0 <- start.params$params$beta0 ## column intercepts
        betas <- NULL
        if (!is.null(X))
          if(!(dim(X) == dim(start.params$X))) stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that predictors X are the same in both models.")
        betas <- c(start.params$params$Xcoef) ## covariates coefficients
        lambdas <- NULL
        lambda2 <- NULL
          lambdas <- start.params$params$theta[,1:num.lv]
          lambdas[upper.tri(lambdas)] <- 0
          if(class(start.params)=="gllvm.quadratic"){
            lambda2 <- start.params$params$theta[,-c(1:num.lv)] 
          }else{
            if(!is.null(X)){
              lambda2 <- fit$params[,(ncol(fit$params)-num.lv+1):ncol(fit$params)]
            }else if(is.null(X)){
              lambda2 <- fit$params[,-c(1:(num.lv+1))]  
            }#this still needs to be implement for traits.
          }
        row.params <- NULL
        if (start.params$row.eff != FALSE) {
          row.params <- start.params$params$row.params
          if(row.params=="fixed")
            row.params[1] <- 0
          if(row.params=="random")
            sigma <- start.params$params$sigma
        }## row parameters
        lvs <- NULL
          lvs <- matrix(start.params$lvs, ncol = num.lv)
      } else {
        stop( "Model which is set as starting parameters isn't the suitable for the one you are trying to fit. Check that attributes y, X and row.eff match to each other.")
      }
    }
    phis <- NULL
    if (family == "negative.binomial") {
      phis <- fit$phi
      if (any(phis > 20))
        phis[phis > 20] <- 20
      if (any(phis < 0.20))
        phis[phis < 0.2] <- 0.2
      fit$phi <- phis
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
    
    if (is.null(offset))
      offset <- matrix(0, nrow = n, ncol = p)
    
    current.loglik <- -1e6; iter <- 1; err <- 10;
    if(!is.null(row.params)){ r0 <- row.params} else {r0 <- rep(0,n)}
    a <- c(beta0)
    b <- NULL; if(!is.null(X)) b <- matrix(betas, ncol(X), p,byrow = TRUE)
      # diag(lambdas) <- log(diag(lambdas)) !!!
      lambda <- lambdas[lower.tri(lambdas,diag = TRUE)]
      u <- lvs
    
    if(!is.null(phis)) { 
      phi <- phis 
    } else { 
      phi <- rep(1, p); 
      fit$phi <- phi
    }
    
    q <- num.lv
    
    optr<-NULL
    timeo<-NULL
    se <- NULL
    
        if(is.null(start.params)  || start.params$method!="VA"){
          if(Lambda.struc=="diagonal" || diag.iter>0){
            Au <- log(rep(Lambda.start[1],num.lv*n)) #1/2, 1
          } else{
            Au <- c(log(rep(Lambda.start[1],num.lv*n)),rep(0,num.lv*(num.lv-1)/2*n)) #1/2, 1
          }
        } else {
          Au <- NULL
          for(d in 1:num.lv) {
            if(start.params$Lambda.struc=="unstructured" || length(dim(start.params$A))==3){
              Au <- c(Au,log(start.params$A[,d,d]))
            } else {
              Au <- c(Au,log(start.params$A[,d]))
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


      if(length(Lambda.start)<2){ Ar <- rep(1,n)} else {Ar <- rep(Lambda.start[2],n)}
      if(row.eff==FALSE){xr <- matrix(0,1,p)} else {xr <- matrix(1,1,p)}
      if(!is.null(X)){Xd <- cbind(1,X)} else {Xd <- matrix(1,n)}
      extra <- 0
      if(family == "poisson") { familyn <- 0}
      if(family == "negative.binomial") { familyn <- 1}
      if(family == "binomial") { familyn <- 2}
      if(family == "ordinal") { familyn <- 3}
      if(row.eff=="random"){
        if(ridge==T){
          if(ridge.quadratic==F){
          objr <- TMB::MakeADFun(
            data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=1, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
            parameters = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0),lambda = lambda, lambda2 = t(lambda2), u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar),zeta=zeta, lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
            inner.control=list(mgcmax = 1e+200,maxit = maxit),
            DLL = "gllvm2")
          }else{
            objr <- TMB::MakeADFun(
              data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=1, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
              parameters = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0),lambda = lambda, lambda2 = t(lambda2), u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar),zeta=zeta, lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
              inner.control=list(mgcmax = 1e+200,maxit = maxit),
              DLL = "gllvm2")
          }
        }else{
          if(ridge.quadratic==F){
            objr <- TMB::MakeADFun(
              data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=0, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
              parameters = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0),lambda = lambda, lambda2 = t(lambda2), u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar),zeta=zeta, lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
              inner.control=list(mgcmax = 1e+200,maxit = maxit),
              DLL = "gllvm2")
            }else{
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=0, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0 = matrix(r0), b = rbind(a,b), B = matrix(0),lambda = lambda, lambda2 = t(lambda2), u = u,lg_phi=log(phi),log_sigma=log(sigma),Au=Au,lg_Ar=log(Ar),zeta=zeta, lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")   
                
              }
          }
          
      } else {
       if(ridge==T){
         if(ridge.quadratic==F){
           objr <- TMB::MakeADFun(
             data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=1, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
             parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, lambda2=t(lambda2), u = u,lg_phi=log(phi),log_sigma=0,Au=Au,lg_Ar=log(Ar),zeta=zeta,lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
             inner.control=list(mgcmax = 1e+200,maxit = maxit),
             DLL = "gllvm2")##GLLVM
         }else{
           objr <- TMB::MakeADFun(
             data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=1, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
             parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, lambda2=t(lambda2), u = u,lg_phi=log(phi),log_sigma=0,Au=Au,lg_Ar=log(Ar),zeta=zeta,lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
             inner.control=list(mgcmax = 1e+200,maxit = maxit),
             DLL = "gllvm2")##GLLVM
         }
       }else{
         if(ridge.quadratic==F){
           objr <- TMB::MakeADFun(
             data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=0, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
             parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, lambda2=t(lambda2), u = u,lg_phi=log(phi),log_sigma=0,Au=Au,lg_Ar=log(Ar),zeta=zeta,lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
             inner.control=list(mgcmax = 1e+200,maxit = maxit),
             DLL = "gllvm2")##GLLVM
         }else{
           objr <- TMB::MakeADFun(
             data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=0, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
             parameters = list(r0=matrix(r0), b = rbind(a,b),B=matrix(0),lambda = lambda, lambda2=t(lambda2), u = u,lg_phi=log(phi),log_sigma=0,Au=Au,lg_Ar=log(Ar),zeta=zeta,lg_gamma=rep(0,num.lv), lg_gamma2=rep(0,num.lv)),
             inner.control=list(mgcmax = 1e+200,maxit = maxit),
             DLL = "gllvm2")##GLLVM
         }
         }
         }
  
      if(optimizer=="nlminb") {
        timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol, iter.max=maxit, eval.max=maxit,trace=trace2)),silent = !trace2))
      }
      if(optimizer=="optim") {
        if(!is.null(par.scale)){
          if(par.scale=="coef"){
            parscale<-abs(objr$par)#this trick comes from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12044
            parscale[parscale==0]<-1
          }else if(is.numeric(par.scale)){
            parscale<-rep(par.scale,length(objr$par))
          }
        }else{
          parscale <- rep(1,length(objr$par))
        }
        if(is.null(fn.scale)|!is.numeric(fn.scale)){
          fnscale<-1
        }else{
          fnscale<-fn.scale
        }
        timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit,parscale=parscale,fnscale=fnscale, trace=trace2),hessian = FALSE),silent = !trace2))
      }
      if(inherits(optr,"try-error")) warning(optr[1]);
      if(diag.iter>0 && Lambda.struc=="unstructured" && num.lv>1 && !inherits(optr,"try-error")){
        objr1 <- objr
        optr1 <- optr
        param1 <- optr$par
        nam <- names(param1)
        r1 <- matrix(param1[nam=="r0"])
        b1 <- matrix(param1[nam=="b"],num.X+1,p)
        lambda1 <- param1[nam=="lambda"]
        lambda2 <- matrix(-1*abs(param1[nam=="lambda2"]),nrow=num.lv)
        u1 <- matrix(param1[nam=="u"],n,num.lv)
        lg_phi1 <- param1[nam=="lg_phi"]
        log_sigma1 <- param1[nam=="log_sigma"]
        #previously  c(pmax(param1[nam=="Au"],rep(log(0.001), num.lv*n)), rep(0.01,num.lv*(num.lv-1)/2*n))
        #this line adds the covariance parameters after diag iter, it didn't start though, this does, sometimes.
        Au1<- c(rep(0,length(param1[names(param1)=="Au"])), rep(0.01,num.lv*(num.lv-1)/2*n))
        lg_Ar1 <- param1[nam=="lg_Ar"]
        lg_gamma <- param1[nam=="lg_gamma"]
        lg_gamma2 <- param1[nam=="lg_gamma2"]
        zeta <- param1[nam=="zeta"]
        
        if(row.eff == "random"){
          if(ridge==T){
            if(ridge.quadratic==F){
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=1, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=log_sigma1,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")
            }else{
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=1, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=log_sigma1,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")
            }
          }else{
            if(ridge.quadratic==F){
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=0, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=log_sigma1,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")
            }else{
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=1, ridge=0, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=log_sigma1,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")
            }
          }
            
        } else {
          if(ridge==T){
            if(ridge.quadratic==F){
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=1, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=0,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")#GLLVM#
            }else{
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=1, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=0,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")#GLLVM#
            }
          }else{
            if(ridge.quadratic==F){
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=0, ridge_quadratic=0, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=0,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")#GLLVM#
            }else{
              objr <- TMB::MakeADFun(
                data = list(y = y, x = Xd,xr=xr,offset=offset, num_lv = num.lv,family=familyn,extra=extra,model=0,random=0, ridge=0, ridge_quadratic=1, trace = as.integer(trace2), zetastruc = ifelse(zeta.struc=="species",1,0)), silent=TRUE,
                parameters = list(r0=r1, b = b1,B=matrix(0),lambda = lambda1, lambda2 = lambda2, u = u1,lg_phi=lg_phi1,log_sigma=0,Au=Au1,lg_Ar=lg_Ar1,zeta=zeta, lg_gamma=lg_gamma, lg_gamma2=lg_gamma2), #log(phi)
                inner.control=list(mgcmax = 1e+200,maxit = maxit),
                DLL = "gllvm2")#GLLVM#
            }
          }
        }
        if(optimizer=="nlminb") {
          timeo <- system.time(optr <- try(nlminb(objr$par, objr$fn, objr$gr,control = list(rel.tol=reltol, iter.max=maxit, eval.max=maxit, trace=trace2)),silent = !trace2))
        }
        if(optimizer=="optim") {
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
          timeo <- system.time(optr <- try(optim(objr$par, objr$fn, objr$gr,method = "BFGS",control = list(reltol=reltol,maxit=maxit,parscale=parscale, fnscale=fnscale, trace=trace2),hessian = FALSE),silent = !trace2))
        }
        if(optimizer=="optim"){
          if(inherits(optr, "try-error") || is.nan(optr$value) || is.na(optr$value)|| is.infinite(optr$value)){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }else{
          if(inherits(optr, "try-error") || is.nan(optr$objective) || is.na(optr$objective)|| is.infinite(optr$objective)){optr=optr1; objr=objr1; Lambda.struc="diagonal"}
        }
      }
      if(inherits(optr,"try-error")) warning(optr[1]);
      param<-objr$env$last.par.best
      if(family =="negative.binomial") {
        phis <- exp(param[names(param)=="lg_phi"])
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
      if(ridge==T){
        gi <- names(param)=="lg_gamma"
      }
      if(ridge.quadratic==T){
        gi2 <- names(param)=="lg_gamma2"
      }
      
      bi <- names(param)=="b"
      li <- names(param)=="lambda"
      l2i <- names(param)=="lambda2"
      ui <- names(param)=="u"
      if(row.eff!=FALSE) {
        ri <- names(param)=="r0"
        if(row.eff=="fixed") row.params <- param[ri]#c(0,param[ri])
        if(row.eff=="random"){ sigma <- exp(param["log_sigma"]); row.params <- param[ri]}
      }
      betaM <- matrix(param[bi],p,num.X+1,byrow=TRUE)
      beta0 <- betaM[,1]
      if(!is.null(X)) betas <- betaM[,-1]
        lvs<-(matrix(param[ui],n,q))
        theta <- matrix(0,p,num.lv)
        if(p>1) {
          theta[lower.tri(theta,diag=TRUE)] <- param[li];
          theta<-cbind(theta,-1*abs(matrix(param[l2i],ncol=num.lv,byrow=T)))
        } else {theta <- param[li]
        theta<-c(theta,param[l2i])}
        # diag(theta) <- exp(diag(theta)) !!!
        
      
      new.loglik <- objr$env$value.best[1]
      if(family %in% c("negative.binomial","gaussian")) {
        phis <- exp(param[names(param)=="lg_phi"])
      }

    if(((n.i==1 || out$logL > abs(new.loglik)) && new.loglik>0) && !inherits(optr, "try-error")){
      out$start <- fit
      objr1 <- objr; optr1=optr;
      out$convergence <- optr1$convergence
      out$logL <- new.loglik
        out$lvs <- lvs
        out$params$theta <- theta
        rownames(out$lvs) <- rownames(out$y);
        if(ridge==T){
          out$params$gamma <- exp(param[gi])
          names(out$params$gamma) <- paste("LV",1:num.lv,sep="")
        }
        if(ridge.quadratic==T){
          out$params$gamma2 <- exp(param[gi2])
          names(out$params$gamma2) <- paste("LV",1:num.lv,"^2",sep="")
        }
        if(num.lv>1) {
          colnames(out$lvs) <- paste("LV", 1:num.lv, sep="")
          colnames(out$params$theta) <- c(paste("LV", 1:num.lv, sep=""),paste("LV", 1:num.lv, "^2", sep=""))
          rownames(out$params$theta) <- colnames(out$y)}
      
      names(beta0) <- colnames(out$y); out$params$beta0 <- beta0;
      if(!is.null(X)){betas <- matrix(betas,ncol=ncol(X)); out$params$Xcoef <- betas;
      rownames(out$params$Xcoef) <- colnames(out$y); colnames(out$params$Xcoef) <- colnames(X); }
      if(family=="ordinal"){
        out$params$zeta <- zetas
      }
      if(family =="negative.binomial") {
        out$params$inv.phi <- phis; names(out$params$inv.phi) <- colnames(out$y);
        out$params$phi <- 1/phis; names(out$params$phi) <- colnames(out$y);
      }
      if(row.eff!=FALSE) {
        if(row.eff=="random"){ out$params$sigma=sigma; names(out$params$sigma)="sigma"}
        out$params$row.params <- row.params; names(out$params$row.params) <- rownames(out$y)
      }
      if(family == "binomial") out$link <- link;
      
      out$row.eff <- row.eff
      out$time <- timeo
      pars <- optr$par
      
      
      param <- objr$env$last.par.best
        Au <- param[names(param)=="Au"]
        A <- array(0,dim=c(n,num.lv,num.lv))
        for (d in 1:num.lv){
          for(i in 1:n){
            A[i,d,d] <- exp(Au[(d-1)*n+i]);
          }
        }
        if(length(Au)>num.lv*n){
          k <- 0;
          for (c1 in 1:num.lv){
            r <- c1+1;
            while (r <=num.lv){
              for(i in 1:n){
                A[i,r,c1] <- Au[num.lv*n+k*n+i];
                A[i,c1,r] <- A[i,r,c1];
              }
              k <- k; r <- r+1;
            }
          }
        }
        out$A <- A
      
      if(row.eff=="random"){
        Ar <- exp(param[names(param)=="lg_Ar"])
        out$Ar <- Ar^2
      }
    }
    
    n.i <- n.i+1;
  }
  tr<-try({
    if(sd.errors && !is.infinite(out$logL)) {
      if(trace) cat("Calculating standard errors for parameters...\n")
      sdr <- optimHess(pars, objr$fn, objr$gr, control = list(reltol=reltol,maxit=maxit))#maxit=maxit
      m <- dim(sdr)[1]; incl <- rep(TRUE,m); incld <- rep(FALSE,m); inclr <- rep(FALSE,m)
      incl[names(objr$par)=="B"] <- FALSE
      
      incl[names(objr$par)=="lg_Ar"] <- FALSE;
      incl[names(objr$par)=="Au"] <- FALSE;
      
      incl[names(objr$par)=="lg_gamma"] <- FALSE;
      incl[names(objr$par)=="lg_gamma2"] <- FALSE;
      
      if(row.eff=="random") {
        inclr[names(objr$par)=="r0"] <- TRUE;
        incl[names(objr$par)=="lg_Ar"] <- FALSE; incld[names(objr$par)=="lg_Ar"] <- TRUE
        incl[names(objr$par)=="r0"] <- FALSE; incld[names(objr$par)=="r0"] <- TRUE
      }
      if(row.eff=="fixed") {incl[1] <- FALSE; incl[names(objr$par)=="log_sigma"] <- FALSE}
      if(row.eff==FALSE) {incl[names(objr$par)=="r0"] <- FALSE; incl[names(objr$par)=="log_sigma"] <- FALSE}
      
        inclr[names(objr$par)=="u"] <- TRUE;
        incl[names(objr$par)=="u"] <- FALSE;
        incld[names(objr$par)=="u"] <- TRUE;
        incld[names(objr$par)=="Au"] <- TRUE;

      if(familyn!=1) incl[names(objr$par)=="lg_phi"] <- FALSE
      if(familyn!=3) incl[names(objr$par)=="zeta"] <- FALSE
        
      A.mat <- -sdr[incl, incl] # a x a
      D.mat <- -sdr[incld, incld] # d x d
      B.mat <- -sdr[incl, incld] # a x d
      cov.mat.mod <- try(MASS::ginv(A.mat-B.mat%*%solve(D.mat)%*%t(B.mat)))
      se <- sqrt(diag(abs(cov.mat.mod)))
      
      incla<-rep(FALSE, length(incl))
      incla[names(objr$par)=="u"] <- TRUE
      out$Hess <- list(Hess.full=sdr, incla = incla, incl=incl, incld=incld, cov.mat.mod=cov.mat.mod)
      
      if(row.eff=="fixed") { se.row.params <- c(0,se[1:(n-1)]); names(se.row.params) <- rownames(out$y); se <- se[-(1:(n-1))] }
      sebetaM <- matrix(se[1:((num.X+1)*p)],p,num.X+1,byrow=TRUE);  se <- se[-(1:((num.X+1)*p))]
        se.lambdas <- matrix(0,p,num.lv); se.lambdas[lower.tri(se.lambdas, diag = TRUE)] <- se[1:(p * num.lv - sum(0:(num.lv-1)))];
        colnames(se.lambdas) <- paste("LV", 1:num.lv, sep="");
        rownames(se.lambdas) <- colnames(out$y)
        out$sd$theta <- se.lambdas; se <- se[-(1:(p * num.lv - sum(0:(num.lv-1))))];
        # diag(out$sd$theta) <- diag(out$sd$theta)*diag(out$params$theta) !!!
        se.lambdas2 <-  matrix(se[1:(p * num.lv)],p,num.lv,byrow=T);
        colnames(se.lambdas2) <- paste("LV", 1:num.lv, "^2",sep="");
        rownames(se.lambdas2) <- colnames(out$y);
        out$sd$theta <- cbind(out$sd$theta,se.lambdas2); se <- se[-(1:(p * num.lv))]
      
      out$sd$beta0 <- sebetaM[,1]; names(out$sd$beta0) <- colnames(out$y);
      if(!is.null(X)){
        out$sd$Xcoef <- matrix(sebetaM[,-1],nrow = nrow(sebetaM));
        rownames(out$sd$Xcoef) <- colnames(y); colnames(out$sd$Xcoef) <- colnames(X);
      }
      if(row.eff=="fixed") {out$sd$row.params <- se.row.params}
      
      if(family %in% c("negative.binomial")) {
        se.lphis <- se[1:p];  out$sd$inv.phi <- se.lphis*out$params$inv.phi;
        out$sd$phi <- se.lphis*out$params$phi;
        names(out$sd$phi) <- colnames(y);  se <- se[-(1:p)]
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
      if(row.eff=="random") { out$sd$sigma <- se*out$params$sigma; names(out$sd$sigma) <- "sigma" }
      
    }})
  if(inherits(tr, "try-error")) { cat("Standard errors for parameters could not be calculated.\n") }
  
  if(is.null(formula1)){ out$formula <- formula} else {out$formula <- formula1}
  
  out$TMBfn <- objr
  out$TMBfn$par <- objr$env$last.par.best #optr$par #ensure params in this fn take final values
  out$logL <- -out$logL
  
  #if(num.lv > 0) out$logL = out$logL + n*0.5*num.lv
  if(row.eff == "random") out$logL = out$logL + n*0.5
  #if(!is.null(randomX)) out$logL = out$logL + p*0.5*ncol(xb)
  return(out)
}

