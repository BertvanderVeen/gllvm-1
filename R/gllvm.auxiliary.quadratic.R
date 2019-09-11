start.values.gllvm.TMB.quadratic <- function(y, X = NULL, TR=NULL, family, 
                                             offset= NULL, trial.size = 1, num.lv = 0, start.lvs = NULL, 
                                             seed = NULL,power=NULL,starting.val="res",formula=NULL, 
                                             jitter.var=0,yXT=NULL, row.eff=FALSE, 
                                             link = "probit", randomX = NULL) {
  if(!is.null(seed)) set.seed(seed)
  N<-n <- nrow(y); p <- ncol(y); y <- as.matrix(y)
  num.T <- 0; if(!is.null(TR)) num.T <- dim(TR)[2]
  num.X <- 0; if(!is.null(X)) num.X <- dim(X)[2]
  Br <- sigmaB <- sigmaij <- NULL
  mu<-NULL
  out <- list()
  
  
  row.params <- rep(0, n);
  if(starting.val %in% c("res","random") || row.eff == "random"){
    rmeany <- rowMeans(y)
    if(family=="binomial"){
      rmeany=1e-3+0.99*rmeany
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  binomial(link = link)$linkfun(rmeany) - binomial(link = link)$linkfun(rmeany[1])
      } else{
        row.params <-  binomial(link = link)$linkfun(rmeany) - binomial(link = link)$linkfun(mean(rmeany))
      }
    } else if(family=="gaussian"){
      rmeany=1e-3+0.99*rmeany
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  rmeany - rmeany[1]
      } else{
        row.params <-  rmeany - mean(rmeany)
      }
    } else {
      if(row.eff %in% c("fixed",TRUE)) {
        row.params <-  row.params <- log(rmeany)-log(rmeany[1])
      } else{
        row.params <-  row.params <- log(rmeany)-log(mean(y))
      }
    }
    if(any(abs(row.params)>1.5)) row.params[abs(row.params)>1.5] <- 1.5 * sign(row.params[abs(row.params)>1.5])
  }
  sigma=1
  if(!is.numeric(y))
    stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
  
  if(!(family %in% c("poisson","negative.binomial","binomial","ordinal")))
    stop("inputed family not allowed...sorry =(")
  
    unique.ind <- which(!duplicated(y))
    if(is.null(start.lvs)) {
      index <- mvtnorm::rmvnorm(N, rep(0, num.lv));
      unique.index <- as.matrix(index[unique.ind,])
    }
    
    if(!is.null(start.lvs)) {
      index <- as.matrix(start.lvs)
      unique.index <- as.matrix(index[unique.ind,])
    }
  
  y <- as.matrix(y)
  
  if(family == "ordinal") {
    max.levels <- apply(y,2,function(x) length(min(x):max(x)));
    if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If al columns only have two levels, please use family == binomial instead. Thanks")
  }
  
  if(is.null(rownames(y))) rownames(y) <- paste("row",1:N,sep="")
  if(is.null(colnames(y))) colnames(y) <- paste("col",1:p,sep="")
  
  options(warn = -1)
  
  if(family!="ordinal") { ## Using logistic instead of prbit regession here for binomial, but whatever...
    if(starting.val=="res" && is.null(start.lvs) ){# && num.lv>0
      if(is.null(TR)){
        if(family!="gaussian") {
          if(!is.null(X)) fit.mva <- mvabund::manyglm(y ~ X, family = family, K = trial.size)
          if(is.null(X)) fit.mva <- mvabund::manyglm(y ~ 1, family = family, K = trial.size)
          resi <- residuals(fit.mva); resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
          coef <- t(fit.mva$coef)
        } else {
          if(!is.null(X)) fit.mva <- mvabund::manylm(y ~ X)
          if(is.null(X)) fit.mva <- mvabund::manylm(y ~ 1)
          resi <- residuals(fit.mva); resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
          coef <- t(fit.mva$coef)
          fit.mva$phi <- apply(fit.mva$residuals,2,sd)
        }
        gamma=NULL
          lastart <- FAstart(mu=NULL, family=family, y=y, num.lv = num.lv, phis=fit.mva$phi, resi=resi)
          gamma<-lastart$gamma
          index<-lastart$index
      } else {
        n1 <- colnames(X)
        n2 <- colnames(TR)
        form1 <-NULL
        if(is.null(formula)){
          form1 <- "~ 1"
          for(m in 1:length(n2)){
            for(l in 1:length(n1)){
              ni <- paste(n1[l],n2[m],sep = "*")
              form1 <- paste(form1,ni,sep = "+")
            }
          }
          formula=form1
        }
        trait.TMB<-getFromNamespace("trait.TMB","gllvm") ##CHANGE THIS LINE ON MERGE
          fit.mva <- trait.TMB(y, X = X, TR = TR, formula=formula(formula), family = family, num.lv = 0, Lambda.struc = "diagonal", trace = FALSE, sd.errors = FALSE, maxit = 1000, seed=seed,n.init=1,starting.val="zero",yXT = yXT, row.eff = row.eff, diag.iter = 0, optimizer = "nlminb", randomX = randomX);
          fit.mva$coef=fit.mva$params
          if(row.eff=="random") sigma=fit.mva$params$sigma
        out$fitstart <- fit.mva
        if(!is.null(form1)){
          if(!is.null(fit.mva$coef$row.params)) row.params=fit.mva$coef$row.params
          env <- (fit.mva$coef$B)[n1]
          trait <- (fit.mva$coef$B)[n2]
          inter <- (fit.mva$coef$B)[!names(fit.mva$coef$B) %in% c(n1,n2)]
          B <- c(env,trait,inter)
        } else {
          B<-fit.mva$coef$B#<-rep(0, length(fit.mva$coef$B));
          if(!is.null(fit.mva$coef$row.params)) row.params=fit.mva$coef$row.params
        }
        
        fit.mva$phi <- phi <- fit.mva$coef$phi
        ds.res <- matrix(NA, n, p)
        rownames(ds.res) <- rownames(y)
        colnames(ds.res) <- colnames(y)
        mu <- (matrix(fit.mva$X.design%*%fit.mva$coef$B,n,p)+ matrix(fit.mva$coef$beta0,n,p,byrow = TRUE))
        if(row.eff %in% c(TRUE, "random", "fixed")) {mu <- mu + row.params }
        if(!is.null(randomX)) {
          Br <- fit.mva$params$Br
          sigmaB <- fit.mva$params$sigmaB
          if(ncol(fit.mva$Xrandom)>1) sigmaij <- fit.mva$params$sigmaB[lower.tri(fit.mva$params$sigmaB)]
          mu <- mu + fit.mva$Xrandom%*%Br
        }
        if(family %in% c("poisson", "negative.binomial")) {
          mu <- exp(mu)
        }
        if(family == "binomial") {
          mu <-  binomial(link = link)$linkinv(mu)
        }
        
        gamma=NULL
          lastart <- FAstart(mu, family=family, y=y, num.lv = num.lv, phis=fit.mva$phi)
          gamma<-lastart$gamma
          index<-lastart$index
      }
      
      if(is.null(TR)){params <- cbind(coef,gamma)
      } else { params <- cbind((fit.mva$coef$beta0),gamma)}
    } else {
      if(family!="gaussian") {
        if(is.null(TR)){
          if(!is.null(X)) fit.mva <- mvabund::manyglm(y ~ X + index, family = family, K = trial.size)
          if(is.null(X)) fit.mva <- mvabund::manyglm(y ~ index, family = family, K = trial.size)
        } else {
          fit.mva <- mvabund::manyglm(y ~ index, family = family, K = trial.size)
          env  <-  rep(0,num.X)
          trait  <-  rep(0,num.T)
          inter <- rep(0, num.T * num.X)
          B <- c(env,trait,inter)
        }
      } else {
        if(is.null(TR)){
          if(!is.null(X)) fit.mva <- mvabund::manylm(y ~ X + index)
          if(is.null(X)) fit.mva <- mvabund::manylm(y ~ index)
        } else {
          fit.mva <- mvabund::manylm(y ~ index)
          env  <-  rep(0,num.X)
          trait  <-  rep(0,num.T)
          inter <- rep(0, num.T * num.X)
          B <- c(env,trait,inter)
        }
        fit.mva$phi <- apply(fit.mva$residuals,2,sd)
      }
      params <- t(fit.mva$coef)
    }}
  
  
  if(family == "negative.binomial") {
    phi <- fit.mva$phi  + 1e-5
  } else if(family == "gaussian") {
    phi <- fit.mva$phi
  } else { phi <- NULL }
  
  if(family == "ordinal") {
    max.levels <- length(unique(c(y)))
    params <- matrix(NA,p,ncol(cbind(1,X))+num.lv)
    zeta <- matrix(NA,p,max.levels - 1)
    zeta[,1] <- 0 ## polr parameterizes as no intercepts and all cutoffs vary freely. Change this to free intercept and first cutoff to zero
    for(j in 1:p) {
      y.fac <- factor(y[,j])
      if(length(levels(y.fac)) > 2) {
        if(starting.val%in%c("zero","res")){
          if(is.null(X) || !is.null(TR)) cw.fit <- MASS::polr(y.fac ~ 1, method = "probit")
          if(!is.null(X) & is.null(TR) ) cw.fit <- MASS::polr(y.fac ~ X, method = "probit")
        } else {
          if(is.null(X) || !is.null(TR)) cw.fit <- MASS::polr(y.fac ~ index, method = "probit")
        }
        params[j,1:ncol(cbind(1,X))] <- c(cw.fit$zeta[1],-cw.fit$coefficients)
        zeta[j,2:length(cw.fit$zeta)] <- cw.fit$zeta[-1]-cw.fit$zeta[1]
      }
      if(length(levels(y.fac)) == 2) {
        if(starting.val%in%c("zero","res")){
          if(is.null(X) || !is.null(TR)) cw.fit <- glm(y.fac ~ 1, family = binomial(link = "probit"))
          if(!is.null(X) & is.null(TR) ) cw.fit <- glm(y.fac ~ X, family = binomial(link = "probit"))
        } else {
          if(is.null(X) || !is.null(TR)) cw.fit <- glm(y.fac ~ index, family = binomial(link = "probit"))
        }
        params[j,] <- cw.fit$coef
      }
    }
    env <- rep(0,num.X)
    trait <- rep(0,num.T)
    inter <- rep(0, num.T * num.X)
    B=c(env,trait,inter)
  }
  
  if(family!="ordinal" || (family=="ordinal" & starting.val=="res")){
    if(num.lv>1 && p>2){
      gamma<-as.matrix(params[,(ncol(params) - num.lv + 1):ncol(params)])
      qr.gamma <- qr(t(gamma))
      params[,(ncol(params) - num.lv + 1):ncol(params)]<-t(qr.R(qr.gamma))
      index<-(index%*%qr.Q(qr.gamma))
    }}
  if(starting.val=="zero"){
    params=matrix(0,p,1+num.X+num.lv)
    params[,1:(ncol(params) - num.lv)] <- 0
    env <- rep(0,num.X)
    trait <- rep(0,num.T)
    inter <- rep(0, num.T * num.X)
    B=c(env,trait,inter)
      gamma <- matrix(1,p,num.lv)
      gamma[upper.tri(gamma)]=0
      params[,(ncol(params) - num.lv + 1):ncol(params)] <- gamma
      index <- matrix(0,n,num.lv)
    phi <- rep(1,p)
  }
    index <- index+mvtnorm::rmvnorm(n, rep(0, num.lv),diag(num.lv)*jitter.var);
    try({
      gamma.new <- as.matrix(params[,(ncol(params) - num.lv + 1):ncol(params)]);
      sig <- sign(diag(gamma.new));
      params[,(ncol(params) - num.lv + 1):ncol(params)] <- t(t(gamma.new)*sig)
      index <- t(t(index)*sig)}, silent = TRUE)
    #add lambda2 here
    #lambda2
    lambda2<-matrix(-1e-5,nrow=p,ncol=num.lv)
    
    if(!is.null(X) & !is.null(TR)){
      quadratic.start.offset <- cbind(1,index)%*%t(params) + matrix(fit.mva$D%*%matrix(B,ncol=1),n,p)
      
      
    }else if(!is.null(X)&is.null(TR)){
      quadratic.start.offset <- cbind(1,X,index)%*%t(params)
    }else if(is.null(X)&is.null(TR)){
      quadratic.start.offset <- cbind(1,index)%*%t(params)
    }
    if(!is.null(offset)){
      quadratic.start.offset <- quadratic.start.offset + offset
    }
    lambda2<-matrix(0,nrow=p,ncol=num.lv)
    if(family!="ordinal"){
      for(j in 1:p){
        if(family!="negative.binomial"){
          lambda2[j,]<-coef(zetadiv::glm.cons(y[,j]~-1+index^2+offset(quadratic.start.offset[,j]),cons=rep(-1,num.lv-1),cons.inter = -1,family=family))  #first coefficient the packages thinks is the intercept
        }else{
          lambda2[j,]<-coef(zetadiv::glm.cons(y[,j]~-1+index^2+offset(quadratic.start.offset[,j]),cons=rep(-1,num.lv-1),cons.inter = -1,family="poisson"))  
        }
      }
    }
    # }else if(family=="ordinal"){
    #   levels<-c(unique(y))
    #   y.start <- y
    #   y.start[y<mean.level] <- 0
    #   y.start[y>=mean.level] <- 1
    #   for(j in 1:p){
    #       lambda2[j,]<-coef(zetadiv::glm.cons(y.start[,j]~-1+lvs^2+offset(beta0[j]+(lvs%*%t(lambdas))[,j]),cons=rep(-1,num.lv-1),cons.inter = -1,family="binomial"))  
    #     }
    # }
    #subtract a fraction from the 0 quadratic scores, otherwise the optimization can't get away from the 0s where necessary.
    lambda2[lambda2==0]<--1e-5
    params <- cbind(params,lambda2)
  
  out$params <- params
  out$phi <- phi
  out$mu <- mu
  if(!is.null(TR)) { out$B <- B}
  out$index <- index
  if(family == "ordinal") out$zeta <- zeta
  options(warn = 0)
  if(row.eff!=FALSE) {
    out$row.params=row.params
    if(row.eff=="random") out$sigma=sigma
  }
  if(!is.null(randomX)){
    out$Br <- Br
    out$sigmaB <- sigmaB
    out$sigmaij <- sigmaij
  }
  return(out)
}                                        





FAstart <- function(mu, family, y, num.lv, zeta = NULL, phis = NULL, 
                    jitter.var = 0, resi = NULL){
  n<-NROW(y); p <- NCOL(y)
  
  if(is.null(resi)){
    ds.res <- matrix(NA, n, p)
    rownames(ds.res) <- rownames(y)
    colnames(ds.res) <- colnames(y)
    for (i in 1:n) {
      for (j in 1:p) {
        if (family == "poisson") {
          a <- ppois(as.vector(unlist(y[i, j])) - 1, mu[i,j])
          b <- ppois(as.vector(unlist(y[i, j])), mu[i,j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "negative.binomial") {
          phis <- phis + 1e-05
          a <- pnbinom(as.vector(unlist(y[i, j])) - 1, mu = mu[i, j], size = 1/phis[j])
          b <- pnbinom(as.vector(unlist(y[i, j])), mu = mu[i, j], size = 1/phis[j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "binomial") {
          a <- pbinom(as.vector(unlist(y[i, j])) - 1, 1, mu[i, j])
          b <- pbinom(as.vector(unlist(y[i, j])), 1, mu[i, j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "gaussian") {
          a <- pnorm(as.vector(unlist(y[i, j])) - 1, 1, mu[i, j], sd = phis[j])
          b <- pnorm(as.vector(unlist(y[i, j])), 1, mu[i, j], sd = phis[j])
          u <- runif(n = 1, min = a, max = b)
          ds.res[i, j] <- qnorm(u)
        }
        if (family == "ordinal") {
          probK <- NULL
          probK[1] <- pnorm(zeta[j,1]-mu[i,j],log.p = FALSE)
          probK[max(y[,j])] <- 1 - pnorm(zeta[j,max(y[,j]) - 1] - mu[i,j])
          if(max(y[,j]) > 2) {
            j.levels <- 2:(max(y[,j])-1)
            for(k in j.levels) { probK[k] <- pnorm(zeta[j,k] - mu[i,j]) - pnorm(zeta[j,k - 1] - mu[i,j]) }
          }
          probK <- c(0,probK)
          cumsum.b <- sum(probK[1:(y[i,j]+1)])
          cumsum.a <- sum(probK[1:(y[i,j])])
          u <- runif(n = 1, min = cumsum.a, max = cumsum.b)
          if (abs(u - 1) < 1e-05)
            u <- 1
          if (abs(u - 0) < 1e-05)
            u <- 0
          ds.res[i, j] <- qnorm(u)
        }
      }
    }
  } else {
    ds.res <- resi
  }
  resi <- as.matrix(ds.res); resi[is.infinite(resi)] <- 0; resi[is.nan(resi)] <- 0
  
  if(p>2 && n>2){
    if(any(is.nan(resi))){stop("Method 'res' for starting values can not be used, when glms fit too poorly to the data. Try other starting value methods 'zero' or 'random' or change the model.")}
    
    if(n>p){
      fa  <-  try(factanal(resi,factors=num.lv,scores = "regression"))
      if(inherits(fa,"try-error")) stop("Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'.")
      gamma<-matrix(fa$loadings,p,num.lv)
      index <- fa$scores
    } else if(n<p) {
      fa  <-  try(factanal(t(resi),factors=num.lv,scores = "regression"))
      if(inherits(fa,"try-error")) stop("Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'.")
      gamma<-fa$scores
      index <- matrix(fa$loadings,n,num.lv)
    } else {
      tryfit <- TRUE; tryi <- 1
      while(tryfit && tryi<5) {
        fa  <-  try(factanal(rbind(resi,rnorm(p,0,0.01)),factors=num.lv,scores = "regression"), silent = TRUE)
        tryfit <- inherits(fa,"try-error"); tryi <- tryi + 1;
      }
      if(tryfit) stop(attr(fa,"condition")$message, "\n Factor analysis for calculating starting values failed. Maybe too many latent variables. Try smaller 'num.lv' value or change 'starting.val' to 'zero' or 'random'.")
      gamma<-matrix(fa$loadings,p,num.lv)
      index <- fa$scores[1:n,]
    }
  } else {
    gamma <- matrix(1,p,num.lv)
    gamma[upper.tri(gamma)]=0
    index <- matrix(0,n,num.lv)
  }
  
  if(num.lv>1 && p>2){
    qr.gamma <- qr(t(gamma))
    gamma.new<-t(qr.R(qr.gamma))
    sig <- sign(diag(gamma.new));
    gamma <- t(t(gamma.new)*sig)
    index<-(index%*%qr.Q(qr.gamma))
    index <- t(t(index)*sig)
  } else {
    sig <- sign(diag(gamma));
    gamma <- t(t(gamma)*sig)
    index <- t(t(index)*sig)
  }
  if(p>n) {
    sdi <- sqrt(diag(cov(index)))
    sdt <- sqrt(diag(cov(gamma)))
    indexscale <- diag(x = 0.8/sdi, nrow = length(sdi))
    index <- index%*%indexscale
    gammascale <- diag(x = 1/sdt, nrow = length(sdi))
    gamma <- gamma%*%gammascale
  }
  index <- index + mvtnorm::rmvnorm(n, rep(0, num.lv),diag(num.lv)*jitter.var);
  return(list(index = index, gamma = gamma))
}

calc.quad <- function(lambda,theta,Lambda.struc) {
  if(Lambda.struc == "diagonal") out <- 0.5 * (lambda) %*% t(theta^2)
  if(Lambda.struc == "unstructured") {
    if(class(lambda) == "array") { n <- dim(lambda)[1]; num.lv <- dim(lambda)[2] }
    if(class(lambda) == "matrix") { num.lv <- dim(lambda)[2]; n <- 1 }
    if(class(theta) == "matrix") { p <- dim(theta)[1] }
    if(class(theta) == "numeric") { p <- 1; theta <- matrix(theta,1) }
    
    out <- matrix(NA,n,p)
    
    if(n == 1) {
      for(j in 1:p) { out[1,j] <- 0.5 * t(theta[j,]) %*% lambda %*% theta[j,] }
    }
    
    if(n > 1) {
      if(n <= p) out <- t(sapply(1:n, function(x) 0.5 * rowSums(theta * (theta %*% lambda[x,,]))))
      if(n > p) {
        lambda.mat <- aperm(lambda,c(3,2,1)); dim(lambda.mat) <- c(num.lv,num.lv * n)
        f <- function(x) 0.5 * rowSums(matrix((t(lambda.mat) * theta[x,]) %*% theta[x,],ncol=num.lv,byrow=TRUE))
        out <- sapply(1:p,f)
      }
    }
  }
  return(list(mat = out, mat.sum = sum(out)))
}


inf.criteria <- function(fit)
{
  family=fit$family
  abund=fit$y
  num.lv=fit$num.lv
  n <- dim(abund)[1]
  k<-attributes(logLik(fit))$df
  
  BIC <- -2*fit$logL + (k) * log(n)
  # AIC
  AIC <- -2*fit$logL + (k) * 2
  # AICc
  AICc <- AIC + 2*k*(k+1)/(n-k-1)
  list(BIC = BIC, AIC = AIC, AICc = AICc, k = k)
}

# Creates matrix of fourth corner terms from a vector
getFourthCorner<- function(object){
  if(is.null(object$X) || is.null(object$TR)) stop();
  
  n1=colnames(object$X)
  n2=colnames(object$TR)
  
  nams=names(object$params$B)
  fx<-cbind(apply(sapply(n1,function(x){grepl(x, nams)}),1,any), apply(sapply(n2,function(x){grepl(x, nams)}),1,any))
  fourth.index<-rowSums(fx)>1
  nams2=nams[fourth.index]
  fourth.corner=object$params$B[fourth.index]
  
  i=1; j=1;
  fourth<-matrix(0,length(n1),length(n2))
  for (i in 1:length(n1)) {
    for (j in 1:length(n2)) {
      fur=(grepl(n1[i], nams2)+grepl(n2[j], nams2))>1
      if(any(fur)){ fourth[i,j]=fourth.corner[fur]}
    }
  }
  colnames(fourth)=n2
  rownames(fourth)=n1
  return(fourth)
}


# Calculates standard errors for random effects
sdrandom<-function(obj, Vtheta, incl, ignore.u = FALSE){
  r <- obj$env$random
  par = obj$env$last.par.best
  hessian.random <- obj$env$spHess(par, random = TRUE)
  L <- obj$env$L.created.by.newton
  if (ignore.u) {
    diag.term2 <- 0
  } else {
    f <- obj$env$f
    w <- rep(0, length(par))
    reverse.sweep <- function(i) {
      w[i] <- 1
      f(par, order = 1, type = "ADGrad", rangeweight = w, doforward = 0)[r]
    }
    nonr <- setdiff(seq_along(par), r)
    tmp <- sapply(nonr, reverse.sweep)
    if (!is.matrix(tmp))
      tmp <- matrix(tmp, ncol = length(nonr))
    A <- solve(hessian.random, tmp[, incl])
    diag.term2 <- diag(rowSums((A %*% Vtheta) * A))
  }
  diag.term1 <- Matrix::chol2inv(L)
  diag.cov.random <- diag.term1 + diag.term2
  return(diag.cov.random)
}


# draw an ellipse
ellipse<-function(center, covM, rad){
  seg <- 51
  Qc <- chol(covM, pivot = TRUE)
  angles <- (0:seg) * 2 * pi / seg
  unit.circ <- cbind(cos(angles), sin(angles))
  order <- order(attr(Qc, "pivot"))
  ellips <- t(center + rad * t(unit.circ %*% Qc[, order]))
  lines(ellips, col = 4)
}

gamEnvelope <- function(x, y,line.col = "red", envelope.col = c("blue","lightblue"), col = 1, envelopes = TRUE, ...){
  xSort <- sort(x, index.return = TRUE)
  gam.yx <- gam(y[xSort$ix] ~ xSort$x)
  pr.y <- predict.gam(gam.yx, se.fit = TRUE)
  n.obs <- length(xSort$ix)
  prHi <- pr.y$fit + 1.96*pr.y$se.fit
  prLow <- pr.y$fit - 1.96*pr.y$se.fit
  if(envelopes) polygon(xSort$x[c(1:n.obs,n.obs:1)], c(prHi,prLow[n.obs:1]), col = envelope.col[2], border = NA)
  lines(xSort$x, pr.y$fit, col = envelope.col[1])
  abline(h = 0, col = 1)
  points(x, y, col = col, ...)
}

