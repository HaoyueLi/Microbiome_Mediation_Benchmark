

extract_sparseDOSSA <- function(sparseDOSSA_fit,dataset,set) {
  
  # metadata + feature data
  if(dataset=="basis"){
    sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_basis)
  }else if(dataset=="count"){
    sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_count)
  }else if(dataset=="norm"){
    sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_norm)
  }else{
    print("No dataset match")
  }
  rownames(sparsedossa_results) <- sparsedossa_results$X1
  nMetadata <- sum(grepl("Metadata", sparsedossa_results$X1, fixed = TRUE))
  nMicrobes <- sum(grepl("Feature_spike", sparsedossa_results$X1, fixed = TRUE))
  nSamples <- ncol(sparsedossa_results) - 1
  sparsedossa_results <- sparsedossa_results[-1, -1]
  sparsedossa_results<-sparsedossa_results[-c(1:nMetadata),]
  colnames(sparsedossa_results) <- paste('Sample', 1:nSamples, sep='')
  
  # data <- as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)), ])
  start_index<-(set-1)*nMicrobes+1
  end_index<-set*nMicrobes
  data<-as.matrix(sparsedossa_results[c(start_index:end_index),])
  data <- data.matrix(data)
  class(data) <- "numeric"
  
  # Extract Features
  data <- as.data.frame(data)
  
  # Rename Features and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
  # wh.TP <- colnames(features) %in% significant_features
  rownames(data) <- paste("Feature", 1:nMicrobes, sep = "")
  # newname <- paste0(colnames(features)[wh.TP], "_TP")
  # colnames(features)[wh.TP] <- newname
  
  
  
  # Return as list
  return(list(data=data))
}

#MedTest
MedOmniTest <-
  function(x, y, m.list, z=NULL, nperm=9999) {
    
    n <- length(y)
    
    save.seed <- get(".Random.seed", .GlobalEnv)
    perm.mat <-  matrix(NA, nrow=nperm, ncol=length(m.list))
    
    obs.vec <- numeric(length(m.list))
    
    if (is.null(z)) {
      x.adj <- resid(lm(x ~ 1))
      z1 <- rep(1, n)
      Px <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
    } else {
      x.adj <- resid(lm(x ~ z))
      z1 <- cbind(1, z)
      Px <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
    }
    
    if (is.null(z)) {
      y.adj <- resid(lm(y ~ x)) 
      z1 <- cbind(1, x)
      Py <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
      
    } else {
      y.adj <- resid(lm(y ~ x + z))
      z1 <- cbind(1, x, z)
      Py <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
    }
    for (j in 1:length(m.list)) {
      assign(".Random.seed", save.seed, .GlobalEnv)
      
      dmat <- as.matrix(m.list[[j]])
      n <- nrow(dmat)
      dmat <- -dmat^2/2
      
      G <- mean(dmat) + dmat - rowMeans(dmat) - matrix(rep(1, n), ncol=1) %*% colMeans(dmat)
      
      obj <- eigen(G, symmetric=TRUE)
      lambda <- obj$values
      u <- obj$vectors
      u <- u[, lambda > 1E-8, drop = FALSE]
      lambda <- lambda[lambda > 1E-8]
      if (is.null(z)) {
        u.adj <- apply(u, 2, function (xx) resid(lm(xx ~ 1)))
        u.x.adj <- apply(u, 2, function (xx) resid(lm(xx ~ x)))
      } else {
        u.adj <- apply(u, 2, function (xx) resid(lm(xx ~ z)))
        u.x.adj <- apply(u, 2, function (xx) resid(lm(xx ~ z + x)))
      }
      
      obs.vec[j] <- sum(lambda * abs(t(u.adj) %*% x.adj) * abs(t(u.x.adj) %*% y.adj))
      perm.mat[, j] <- sapply(1:nperm, function(i) {
        x.adj.p <- Px %*% x.adj[sample(n)]
        y.adj.p <- Py %*% y.adj[sample(n)]
        f1.f2.p <- sum(lambda * abs(t(u.adj) %*% x.adj) * abs(t(u.x.adj) %*% y.adj.p))
        f2.f1.p <- sum(lambda * abs(t(u.adj) %*% x.adj.p) * abs(t(u.x.adj) %*% y.adj))
        f1.p.f2.p <- sum(lambda * abs(t(u.adj) %*% x.adj.p) * 
                           abs(t(u.x.adj) %*% y.adj.p))
        max(f1.f2.p, f2.f1.p, f1.p.f2.p)
      })
    }
    # marginal p values
    margPs <- colMeans(rbind(sweep(perm.mat, 2, obs.vec, '>'), rep(1, length(m.list)))) 
    
    
    perm.mat.p <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
    
    
    return(list(margPs = margPs, 
         #margP = min(min(margPs * length(m.list)), 1),
         permP = mean(c(rowMins(perm.mat.p) <= min(margPs), 1))))
    
  }

#Huang and Pan
## Test for mediation effect for continuous outcome ##
IEtest.omnibus<-function(
  G=gset,			# mediators (n-by-p matrix)
  S=snp,			# exposure (n-by-1 vector or matrix)
  Y=diz,			# outcome (n-by-1 vector or matrix)
  X=conf,			# covariates (n-by-(q-1) matrix, not including intercept)
  n.draw=1000,		# no. of Monte-Carlo resampling
  adaptive=FALSE,		# TRUE if selecting the number of transformed mediators up to 80% variability
  small=NA			# TRUE: small sample method based on transformation models; FALSE: full results (both original and transformation models)
){
  ########################################################################
  
  n<-dim(G)[1]
  p<-dim(G)[2]
  if (sum(is.na(X))) {
    q<-1
  } else if (is.null(dim(X))) {
    q<-2
  } else q<-dim(X)[2]+1
  if (is.na(small)) small<-n<5*p
  
  ### Estimate U ###
  
  fit.g<-lm(G~X+S)
  fit.y<-lm(Y~X+S+G+S*G)
  if (n-2*p-q-1>=5) sig.y<-sum(fit.y$res^2)/(n-2*p-q-1)
  Sigma<-cov(fit.g$residual)
  svds<-svd(Sigma)
  if (adaptive){
    n.factor<-sum(cumsum(svds$d)/sum(svds$d)<0.8)+1
    U<-svds$u[,1:n.factor]
    p.prime<-n.factor
  } else {U<-svds$u; p.prime=p}
  P<-G%*%U
  
  Ut<-t(U)
  UU<-bdiag(Ut, Ut, Ut)
  
  ### Calculate component-wise S-P-Y p values using delta method ###
  
  fit.g2<-lm(P~X+S)
  
  s.delta<-diag(vcov(fit.g2))[(1:p.prime)*3]
  thetas.p<-NULL; s.thetas.p<-matrix(NA); mse2<-NULL
  for (j in 1:p.prime){
    Pj<-P[,j]
    fit.y2pj<-lm(Y~X+S+X*S+I(S^2)+Pj+S*Pj)
    thetas.p<-rbind(thetas.p, fit.y2pj$coef[c("Pj", "S:Pj")])
    mse2<-c(mse2, summary(fit.y2pj)$sigma^2)
  }
  
  cov1.p<-diag(s.delta)
  Z<-cbind(1, X, S, P, S*P)
  info<-min(mse2)*(t(Z)%*%Z)
  fit.y<-lm(Y~X+S+G+S*G)
  if (!small) info<-(summary(fit.y)$sigma)^2*(t(Z)%*%Z)			## revision 2 ##
  cov2.p<-solve(info)[q+1+1:(2*p.prime), q+1+1:(2*p.prime)]
  
  bcov.p<-bdiag(cov1.p, cov2.p)
  theta.sum.p<-apply(thetas.p, 1, sum)
  beta.p<-fit.g2$coef["S",]
  
  pam.i.p<-rbind(diag(theta.sum.p), diag(beta.p)[rep(1:p.prime, 2), 1:p.prime])
  ie.p<-theta.sum.p*beta.p
  vie.p<-t(pam.i.p)%*%bcov.p%*%pam.i.p
  stat.p<-(t(ie.p)%*%solve(vie.p)%*%ie.p)[1]
  pval.p<-pchisq(stat.p, df=p.prime, lower.tail=FALSE)
  
  pam.p<-c(theta.sum.p, rep(beta.p, 2))
  vie.p.all<-(t(pam.p)%*%bcov.p%*%pam.p)[1]
  ie.p.all<-sum(ie.p)
  pval.p.sum<-pchisq(ie.p.all^2/vie.p.all, df=1, lower.tail=FALSE)
  
  ### Calculate component-wise and total S-G-Y p values using delta method ###
  
  if (!small){
    n<-dim(G)[1]
    p<-dim(G)[2]
    p.ori<-p
    fit.g<-lm(G~X+S)
    fit.y<-lm(Y~X+S+G+S*G)
    cov1<-vcov(fit.g)[(1:p)*3, (1:p)*3]
    cov2<-vcov(fit.y)[(1:(2*p))+3, (1:(2*p))+3]
    bcov<-bdiag(cov1, cov2)
    betga<-fit.y$coef[(1:(2*p))+3]
    theta.sum<-apply(matrix(fit.y$coef[(1:(2*p))+3], nr=p), 1, sum)
    beta<-fit.g$coef["S",]
    
    pam.i<-rbind(diag(theta.sum), diag(beta), diag(beta))
    ie<-theta.sum*beta
    vie<-t(pam.i)%*%bcov%*%pam.i
    stat.comp<-(t(ie)%*%solve(vie)%*%ie)[1]
    pval.comp<-pchisq(stat.comp, df=p, lower.tail=FALSE)
    
    pam<-c(theta.sum, rep(beta, 2))
    vie.all<-(t(pam)%*%bcov%*%pam)[1]
    ie.all<-sum(ie)
    pval.sum<-pchisq(ie.all^2/vie.all, df=1, lower.tail=FALSE)
    
    ### Calculate S-G-Y p values using resampling method ###
    
    set.seed(137)
    s.theta.mc<-as.matrix(bcov)
    theta.mc<-rmvnorm(n.draw, mean=c(beta, betga), sigma=s.theta.mc)	## sample from the distributions of regression coefficients
    delta.mc<-theta.mc[,1:p]
    betag.mc<-theta.mc[,p+1:p]
    gamma.mc<-theta.mc[,2*p+1:p]
    ies.mc<-((betag.mc+gamma.mc)*delta.mc)
    
    stat.mc<-t(ie-0)%*%solve(cov(ies.mc))%*%(ie-0)
    pval.mc<-pchisq(stat.mc, df=p.ori, lower.tail=FALSE)
    
    ies.mc.c<-t(t(ies.mc)-apply(ies.mc, 2, mean))
    vie.mc<-apply(ies.mc.c^2, 1, sum)	
    pval.v.mc<-mean(vie.mc>sum(ie^2))
    
    tie.mc<-apply(ies.mc, 1, sum)
    tie.mc.c<-tie.mc-mean(tie.mc)
    pval.sum.mc<-2*mean(tie.mc.c>abs(sum(ie)))
    
  }
  
  ### Calculate S-P-Y p values using resampling method ###
  
  theta.p.mc<-rmvnorm(n.draw, mean=c(beta.p, thetas.p), sigma=as.matrix(bcov.p))	## sample from the distribution of regression coefficients
  delta.p.mc<-theta.p.mc[,1:p.prime]
  betag.p.mc<-theta.p.mc[,p.prime+1:p.prime]
  gamma.p.mc<-theta.p.mc[,2*p.prime+1:p.prime]
  ies.p.mc<-((betag.p.mc+gamma.p.mc)*delta.p.mc)
  
  ies.p.mc.s<-cov(ies.p.mc)
  stat.p.mc<-t(ie.p-0)%*%solve(ies.p.mc.s)%*%(ie.p-0)
  pval.p.mc<-pchisq(stat.p.mc, df=p.prime, lower.tail=FALSE)
  
  ies.p.mc.c<-t(t(ies.p.mc)-apply(ies.p.mc, 2, mean))
  vie.p.mc<-apply(ies.p.mc.c^2, 1, sum)	
  pval.v.p.mc<-mean(vie.p.mc>sum(ie.p^2))
  
  tie.p.mc<-apply(ies.p.mc, 1, sum)
  tie.p.mc.c<-tie.p.mc-mean(tie.p.mc)
  pval.p.sum.mc<-2*mean(tie.p.mc.c>abs(sum(ie.p)))
  
  ### Calculate omnibus p value ###
  
  if (!small){
    
    cov.1<-cov(ies.mc)
    var.1sum<-var(apply(ies.mc, 1, sum))
    cov.1.inv<-solve(cov.1)
    set.seed(37)
    ie.null.1<-ies.mc.c
    ie.null.1sum<-apply(ie.null.1, 1, sum)
    vie.null.1<-apply(ie.null.1^2, 1, sum)
    
    cov.inv.1.svd<-svd(cov.1.inv)
    cov.inv.1.5<-cov.inv.1.svd$u%*%diag(sqrt(cov.inv.1.svd$d))
    stat.1.mc.null<-apply((ie.null.1%*%cov.inv.1.5)^2, 1, sum)
    pval.1.mc.null<-pchisq(stat.1.mc.null, df=p.ori, lower.tail=FALSE)
    stat.1sum.mc.null<-ie.null.1sum^2/var.1sum
    pval.1sum.mc.null<-pchisq(stat.1sum.mc.null, df=1, lower.tail=FALSE)
    pval.1.vie.mc.null<-(rank(vie.null.1)-0.5)/(n.draw)
  }
  
  cov.2<-cov(ies.p.mc)
  var.2sum<-var(apply(ies.p.mc, 1, sum))
  cov.2.inv<-solve(cov.2)
  set.seed(37)
  ie.null.2<-ies.p.mc.c
  ie.null.2sum<-apply(ie.null.2, 1, sum)
  vie.null.2<-apply(ie.null.2^2, 1, sum)
  
  cov.inv.2.svd<-svd(cov.2.inv)
  cov.inv.2.5<-cov.inv.2.svd$u%*%diag(sqrt(cov.inv.2.svd$d))
  stat.2.mc.null<-apply((ie.null.2%*%cov.inv.2.5)^2, 1, sum)
  pval.2.mc.null<-pchisq(stat.2.mc.null, df=p.prime, lower.tail=FALSE)
  stat.2sum.mc.null<-ie.null.2sum^2/var.2sum
  pval.2sum.mc.null<-pchisq(stat.2sum.mc.null, df=1, lower.tail=FALSE)
  pval.2.vie.mc.null<-(rank(vie.null.2)-0.5)/(n.draw)
  
  if (!small){
    pval.min.null2<-apply(cbind(pval.1sum.mc.null, pval.2sum.mc.null, pval.2.vie.mc.null, pval.1.vie.mc.null), 1, min)
    pval.omnibus2<-mean(pval.min.null2<min(pval.sum.mc, pval.p.sum.mc, pval.v.mc, pval.v.p.mc))
    pval<-c(pval.sum, pval.sum.mc, pval.p.sum, pval.p.sum.mc, pval.comp, pval.mc, pval.p, pval.p.mc, pval.v.mc, pval.v.p.mc,pval.omnibus2)
    names(pval)<-c("pval.Delta (delta mthd)", "pval.Delta (MC)", "pval.Delta.trsfrm (delta mthd)", "pval.Delta.trsfrm (MC)", 
                   "pval.deltas (delta mthd)", "pval.deltas (MC)", "pval.deltas.trsfrm (delta mthd)", "pval.deltas.trsfrm (MC)", 
                   "pval.tau (MC)", "pval.tau.p (MC)", "pval.omnibus")
    
  }
  
  if (small){
    pval.min.null1<-apply(cbind(pval.2sum.mc.null, pval.2.vie.mc.null), 1, min)
    pval.omnibus1<-mean(pval.min.null1<=min(pval.p.sum.mc, pval.v.p.mc))
    pval<-c(pval.p.sum, pval.p.sum.mc, pval.p, pval.p.mc, pval.v.p.mc, pval.omnibus1)
    names(pval)<-c("pval.Delta.trsfrm (delta mthd)", "pval.Delta.trsfrm (MC)", 
                   "pval.deltas.trsfrm (delta mthd)", "pval.deltas.trsfrm (MC)", "pval.tau.p (MC)", "pval.omnibus")
    
  }
  
  return(pval)
  
}

IEtest.omnibus.boot<-function(
  G0=gset,			# mediators (n-by-p matrix)
  S0=snp,			# exposure (n-by-1 vector or matrix)
  Y0=diz,			# outcome (n-by-1 vector or matrix)
  X0=conf,			# covariates (n-by-(q-1) matrix, not including intercept
  n.draw=1000,		# no. of bootstrapping
  small=NA			# TRUE: small sample method based on transformation model; FALSE: full results (both original and transformation models)
){
  ########################################################################
  
  n<-dim(G0)[1]
  p<-dim(G0)[2]
  if (sum(is.na(X0))) {
    q<-1
  } else if (is.null(dim(X0))) {
    q<-2
  } else q<-dim(X0)[2]+1
  if (is.na(small)) small<-n<5*p
  
  ### Estimate U ###
  
  fit.g0<-lm(G0~X0+S0)
  Sigma<-cov(fit.g0$residual)
  svds<-svd(Sigma)
  U<-svds$u; p.prime=p
  
  theta.p.mc<-matrix(0, nc=3*p, nr=n.draw+1)
  theta.mc<-matrix(0, nc=3*p, nr=n.draw+1)
  for (bb in 0:n.draw){
    
    set.seed(bb*13)
    if (bb==0) samp<-1:n
    if (bb>0) samp<-sample(1:n, re=TRUE)
    X<-X0[samp]
    S<-S0[samp]
    G<-G0[samp,]
    Y<-Y0[samp,]
    
    fit.y<-lm(Y~X+S+G+S*G)
    P<-G%*%U
    Ut<-t(U)
    UU<-bdiag(Ut, Ut, Ut)
    
    ### Fit S-Pj-Y ###
    
    fit.g2<-lm(P~X+S)
    thetas.p<-NULL; s.thetas.p<-matrix(NA); mse2<-NULL
    for (j in 1:p.prime){
      Pj<-P[,j]
      fit.y2pj<-lm(Y~X+S+X*S+I(S^2)+Pj+S*Pj)
      thetas.p<-rbind(thetas.p, fit.y2pj$coef[c("Pj", "S:Pj")])
      mse2<-c(mse2, summary(fit.y2pj)$sigma^2)
    }
    beta.p<-fit.g2$coef["S",]
    theta.p.mc[bb+1, 1:(p.prime*3)]<-c(beta.p, thetas.p)
    
    ### Fit S-G-Y ###
    
    if (!small){
      n<-dim(G)[1]
      p<-dim(G)[2]
      p.ori<-p
      fit.g<-lm(G~X+S)
      fit.y<-lm(Y~X+S+G+S*G)
      betga<-fit.y$coef[(1:(2*p))+3]
      theta.sum<-apply(matrix(fit.y$coef[(1:(2*p))+3], nr=p), 1, sum)
      beta<-fit.g$coef["S",]
      theta.mc[bb+1, ]<-c(beta, betga)
    }
    print(bb)
    flush.console()
  }
  
  ### Calculate S-G-Y p values using resampling method ###
  
  if (!small){
    
    delta.mc<-theta.mc[,1:p]
    betag.mc<-theta.mc[,p+1:p]
    gamma.mc<-theta.mc[,2*p+1:p]
    ies.mc<-((betag.mc+gamma.mc)*delta.mc)
    ie<-ies.mc[1,]
    
    ies.mc.c<-t(t(ies.mc)-apply(ies.mc, 2, mean))
    vie.mc<-apply(ies.mc.c^2, 1, sum)	
    pval.v.mc<-mean(vie.mc>sum(ie^2))
    
    tie.mc<-apply(ies.mc, 1, sum)
    tie.mc.c<-tie.mc-mean(tie.mc)
    pval.sum.mc<-2*mean(tie.mc.c>abs(sum(ie)))
    
  }
  
  ### Calculate S-P-Y p values using resampling method ###
  
  delta.p.mc<-theta.p.mc[,1:p.prime]
  betag.p.mc<-theta.p.mc[,p.prime+1:p.prime]
  gamma.p.mc<-theta.p.mc[,2*p.prime+1:p.prime]
  ies.p.mc<-((betag.p.mc+gamma.p.mc)*delta.p.mc)
  ie.p<-ies.p.mc[1,]
  
  ies.p.mc.c<-t(t(ies.p.mc)-apply(ies.p.mc, 2, mean))
  vie.p.mc<-apply(ies.p.mc.c^2, 1, sum)	
  pval.v.p.mc<-mean(vie.p.mc>sum(ie.p^2))
  
  tie.p.mc<-apply(ies.p.mc, 1, sum)
  tie.p.mc.c<-tie.p.mc-mean(tie.p.mc)
  pval.p.sum.mc<-2*mean(tie.p.mc.c>abs(sum(ie.p)))
  
  
  if (!small){ 
    pval<-c(pval.sum.mc, pval.p.sum.mc, pval.v.mc, pval.v.p.mc)
    names(pval)<-c("pval.Delta (bootstrap)", "pval.Delta.trsfrm (bootstrap)", 
                   "pval.deltas (bootstrap)", "pval.deltas.trsfrm (bootstrap)")
  }
  if (small){
    pval<-c(pval.p.sum.mc, pval.v.p.mc)
    names(pval)<-c("pval.Delta.trsfrm (bootstrap)", "pval.deltas.trsfrm (bootstrap)")
  }
  return(pval)
  
}

#Helper functions for HDMA: https://www.frontiersin.org/articles/10.3389/fgene.2019.01195/full#T2
# Internal function: parallel computing check
checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel & (ncore > 1)) {
    if (ncore > parallel::detectCores()) {
      message("You requested ", ncore, " cores. There are only ", 
              parallel::detectCores(), " in your machine!")
      ncore <- parallel::detectCores()
    }
    if (verbose) 
      message("    Running ", program.name, " with ", ncore, " cores in parallel...   (", 
              Sys.time(), ")")
    doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              Sys.time(), ")")
    registerDoSEQ()
  }
}

## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", 
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ", 
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text, 
                   "]};invisible(b)}")
  return(script)
}

## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}

## Internal function: scale data (obsolete function)

scaleto <- function(dat) {
  if (is.null(dat)) 
    return(list(dn = NULL, d = NULL, ds = NULL))
  dat_scale <- scale(dat)
  dat_names <- names(dat)
  if (any(class(dat) %in% c("matrix", "data.frame", "data.table"))) {
    dat_names <- colnames(dat)
    dat <- as.matrix(data.frame(dat_scale))
  } else {
    dat_names <- names(dat)
    dat <- as.numeric(dat_scale)
  }
  dat_scale <- as.numeric(attributes(dat_scale)[["scaled:scale"]])
  return(list(dn = dat_names, d = dat, ds = dat_scale))
}

# Internal function: Sure Independent Screening
# Global variables:
globalVariables("n")
globalVariables("M_chunk")

himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel, ncore, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }
  
  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                             glm.family, ", data = datarun))"), "c(1,4)")))
  
  checkParallel(tag, parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}
  
  colnames(results) <- M.names
  return(results)
}


#HDMA: https://www.frontiersin.org/articles/10.3389/fgene.2019.01195/full#T2
####### Load three internal functions before using hdma function to do mediation analysis.
####### Three internal functions are established by YinanZheng (2016) and posted in the website of 
####### https://github.com/YinanZheng/HIMA/blob/master/R/utils.R.  
####### There is need to download the script of utils.R from the website and then load into your workplace of R. 
####### The hdma function is aimed at doing high-dimensional mediation analysis.
####### R version: R/3.5.1.
hdma <- function (X, Y, M, COV.XM = NULL, COV.MY = COV.XM, family = c("gaussian","binomial"), method = c("lasso", "ridge"), topN = NULL,
                  parallel = FALSE, ncore = 1, verbose = FALSE, ...){	
  ####################################################################################################################################
  #########################################                   Function body                ###########################################
  ####################################################################################################################################
  ####### INPUT
  ####### X : Independent variable that is a vector
  ####### Y : Dependent variable that is a vector and can be either continuous or binary variable
  ####### M : High-dimensional mediators that can be either data.frame or matrix. Rows represent samples, columns represent variables
  ####### COV.XM : a data.frame or matrix of covariates dataset for testing the association X ~ M. Default = NULL. 
  #######          If the covariates contain mixed types, please make sure all categorical variables are properly transformed into factor
  #######          type.
  ####### COV.MY : a data.frame or matrix of covariates dataset for testing the association Y ~ M. Using covariates should be careful.
  #######          If the cavariables are not specified, the covariates for Y ~ M are the same with that of M ~ X.
  ####### family : either 'gaussian' or 'binomial', relying on the type of outcome (Y). See hdi package.
  ####### method : either "lasso" or "ridge" to estimate the effect of M -> Y.
  ####### topN : an integer can be used to set the number of top markers by the method of sure independent screening. Default = NULL.
  #######        If topN is NULL, it will be either ceiling(n/log(n)) if family = 'gaussian', or ceiling(n/(2*log(n))) if family = 
  #######	      'binomial', where n is the sample size. If the sample size is greater than topN (pre-specified or calculated), all
  #######        markers will be harbored in the test.
  ####### parallel : logical parameter. The parameter can be employed to enable your computer to do parallel calculation. Default = FALSE.
  ####### ncore : the parameter can be used to set the number of cores to run parallel computing when parallel == TRUE. By default max
  ########	number of cores available in the machine will be utilized.
  ####### verbose : logical. Default = FALSE.
  ####### ... : other arguments passed to hdi.
  ####################################################################################################################################
  ####### Values 
  ####### alpha : the coefficient can reflect the association of X –> M, note that the effect is adjusted by covariables when covariables
  #######	are not NULL.
  ####### beta : the coefficient can reflect the association of M –> Y, note that the effect is adjusted by X. When covariables are not
  #######	NULL, the effect is adjusted by X and covariables. 
  ####### gamma : the coefficient can reflect the linkage of X –> Y, it can represent the total effect.
  ####### alpha*beta : the estimator of mediation effect.
  ####### %total effect : alpha*beta/gamma*100. The proportion of the mediation effect is out of the total effect.
  ####### p-values : joint significant test for mediators.
  ####################################################################################################################################
  ####### checking the necessary packages 
  pkgs <- list("hdi","MASS","doParallel", "foreach","iterators")
  checking<-unlist(lapply(pkgs, require, character.only = T))
  if(any(checking==F))
    stop("Please install the necessary packages first!")	
  family <- match.arg(family)
  method <- match.arg(method)
  if (parallel & (ncore == 1)) ncore <- parallel::detectCores()
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2*n/log(n)) 
  } else {
    d <- topN      # the number of top mediators that associated with independent variable (X)
  }
  d <- min(p, d)   # if d > p select all mediators
  #############################################################################################################################
  ################################           Step-1 Sure Independent Screening (SIS)          #################################
  #############################################################################################################################
  message("Step 1: Sure Independent Screening ...", " (", Sys.time(), ")")
  if(family == "binomial") 
  {
    if(verbose) message("Screening M using the association between X and M: ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", parallel = parallel, 
                                   ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", parallel = parallel,
                           ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)						 
  #################################################################################################################################
  ####################################          Step-2  High-dimensional Inference (HDI)         ##################################
  #################################################################################################################################
  message("Step 2: High-dimensional inference (", method, ") ...", "     (", Sys.time(), ")")
  ## Based on the SIS results in step 1. We will find the most influential M on Y.	
  if (is.null(COV.MY)) {
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM, Y, family = family) else fit <- ridge.proj(XM, Y, family = family)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    set.seed(1029)
    if (method == "lasso") fit <- lasso.proj(XM_COV, Y, family = family) else fit <- ridge.proj(XM_COV, Y, family = family)
  }
  P_hdi<-fit$pval[1:length(ID)]
  index<-which(P_hdi<=0.05)
  if(verbose)  message("Non-zero ",method, " beta estimate(s) of mediator(s) found: ", paste0(names(index), collapse = ","))
  if(length(index)==0)
    return(NULL)
  ID_test <- ID[index]  
  if(family == "binomial")
  {
    ## This has been done in step 1 (when Y is binary)
    alpha_est <- alpha[,ID_test, drop = FALSE]
  } else {
    if(verbose) message(" Estimating alpha (effect of X on M): ", appendLF = FALSE)
    alpha_est <- himasis(NA, M[,ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                         parallel = FALSE, ncore = ncore, verbose, tag = "site-by-site ordinary least squares estimation")
  }   
  beta_P<-P_hdi[index]
  beta_hat<-fit$bhat[index]              # the estimator for beta
  alpha_hat<-as.numeric(alpha_est[1, ])
  ab_est<-beta_hat*alpha_hat
  alpha_P<-alpha_est[2,]
  PA <- rbind(beta_P,alpha_P)
  P_value <- apply(PA,2,max)
  ###################################################################################################################################    
  ###############################          STEP 3   Computing the propotion of mediation effect          ############################
  ###################################################################################################################################
  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
  results <- data.frame(alpha = alpha_hat, beta = beta_hat, gamma = gamma_est, `alpha*beta` = ab_est, `%total effect` 
                        =ab_est/gamma_est*100, `P.value` = P_value, check.names = FALSE)
  message("Done!", " (", Sys.time(), ")")
  doParallel::stopImplicitCluster()
  return(results)
}

TSSnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # TSS Normalizing the Data
  features_TSS <-
    vegan::decostand(
      features_norm,
      method = "total",
      MARGIN = 1,
      na.rm = TRUE)
  
  # Convert back to data frame
  features_TSS <- as.data.frame(features_TSS)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TSS) <- dd
  
  
  # Return
  return(features_TSS)
}


######################
## CLR Normalization #
######################

# Apply CLR Normalization To A Dataset

CLRnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # CLR Normalizing the Data
  features_CLR <- chemometrics::clr(features_norm + 1)
  
  # Convert back to data frame
  features_CLR <- as.data.frame(features_CLR)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CLR) <- dd
  
  
  # Return
  return(features_CLR)
}

######################
## CSS Normalization #
######################

# Apply CSS Normalization To A Dataset

CSSnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # CSS Normalizing the Data
  # Create the metagenomeSeq object
  MGS = metagenomeSeq::newMRexperiment(
    t(features_norm),
    featureData = NULL,
    libSize = NULL,
    normFactors = NULL
  )
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = metagenomeSeq::cumNorm(MGS, p = metagenomeSeq::cumNormStat(MGS))
  # Save the normalized data as data.frame
  features_CSS = as.data.frame(t(
    metagenomeSeq::MRcounts(MGS, norm = TRUE, log = FALSE)))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CSS) <- dd
  
  
  # Return as list
  return(features_CSS)
}


#' High-dimensional Mediation Analysis
#' 
#' \code{hima} is used to estimate and test high-dimensional mediation effects.
#' 
#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a \code{data.frame} or \code{matrix} of high-dimensional mediators. Rows represent samples, columns 
#' represent variables.
#' @param COV.XM a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{M ~ X}. 
#' Covariates specified here will not participate penalization. Default = \code{NULL}. If the covariates 
#' contain mixed data types, please make sure all categorical variables are properly formatted as \code{factor} type.
#' @param COV.MY a \code{data.frame} or \code{matrix} of covariates dataset for testing the association \code{Y ~ M}. 
#' Covariates specified here will not participate penalization. If not specified, the same set of covariates for 
#' \code{M ~ X} will be applied. Using different sets of covariates is allowed but this needs to be handled carefully.
#' @param family either 'gaussian' or 'binomial', depending on the data type of outcome (\code{Y}). See 
#' \code{\link{ncvreg}}
#' @param penalty the penalty to be applied to the model. Either 'MCP' (the default), 'SCAD', or 
#' 'lasso'. See \code{\link{ncvreg}}.
#' @param topN an integer specifying the number of top markers from sure independent screening. 
#' Default = \code{NULL}. If \code{NULL}, \code{topN} will be either \code{ceiling(n/log(n))} if 
#' \code{family = 'gaussian'}, or \code{ceiling(n/(2*log(n)))} if \code{family = 'binomial'}, 
#' where \code{n} is the sample size. If the sample size is greater than topN (pre-specified or calculated), all 
#' mediators will be included in the test (i.e. low-dimensional scenario).
#' @param parallel logical. Enable parallel computing feature? Default = \code{TRUE}.
#' @param ncore number of cores to run parallel computing Valid when \code{parallel == TRUE}. 
#' By default max number of cores available in the machine will be utilized.
#' @param verbose logical. Should the function be verbose? Default = \code{FALSE}.
#' @param ... other arguments passed to \code{\link{ncvreg}}.
#' 
#' @return A data.frame containing mediation testing results of selected mediators. 
#' \itemize{
#'     \item{alpha: }{coefficient estimates of exposure (X) --> mediators (M).}
#'     \item{beta: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for exposure).}
#'     \item{gamma: }{coefficient estimates of exposure (X) --> outcome (Y) (total effect).}
#'     \item{alpha*beta: }{mediation effect.}
#'     \item{\% total effect: }{alpha*beta / gamma. Percentage of the mediation effect out of the total effect.}
#'     \item{adjusted.p: }{statistical significance of the mediator (Bonferroni procedure).}
#'     \item{BH.FDR: }{statistical significance of the mediator (Benjamini-Hochberg procedure).}
#' }
#'
#' @examples
#' n <- 100  # sample size
#' p <- 500 # the dimension of covariates
#' 
#' # the regression coefficients alpha (exposure --> mediators)
#' alpha <- rep(0, p) 
#' 
#' # the regression coefficients beta (mediators --> outcome)
#' beta1 <- rep(0, p) # for continuous outcome
#' beta2 <- rep(0, p) # for binary outcome
#' 
#' # the first four markers are true mediators
#' alpha[1:4] <- c(0.45,0.5,0.6,0.7)
#' beta1[1:4] <- c(0.55,0.6,0.65,0.7)
#' beta2[1:4] <- c(1.45,1.5,1.55,1.6)
#'
#' # these are not true mediators
#' alpha[7:8] <- 0.5
#' beta1[5:6] <- 0.8
#' beta2[5:6] <- 1.7
#' 
#' # Generate simulation data
#' simdat_cont = simHIMA(n, p, alpha, beta1, seed=1029) 
#' simdat_bin = simHIMA(n, p, alpha, beta2, binaryOutcome = TRUE, seed=1029) 
#' 
#' # Run HIMA with MCP penalty by default
#' # When Y is continuous (default)
#' hima.fit <- hima(simdat_cont$X, simdat_cont$Y, simdat_cont$M, verbose = TRUE) 
#' head(hima.fit)
#' 
#' # When Y is binary (should specify family)
#' hima.logistic.fit <- hima(simdat_bin$X, simdat_bin$Y, simdat_bin$M, 
#' family = "binomial", verbose = TRUE) 
#' head(hima.logistic.fit)
#' 
#' @export
my_hima <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM, 
                 family = c("gaussian", "binomial"), 
                 penalty = c("MCP", "SCAD", "lasso"), 
                 topN = NULL, 
                 parallel = FALSE, 
                 ncore = 1, 
                 verbose = FALSE, 
                 ...) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  
  if (parallel & (ncore == 1)) ncore <- parallel::detectCores()
  
  n <- nrow(M)
  p <- ncol(M)
  
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
  
  d <- min(p, d) # if d > p select all mediators
  
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Sure Independent Screening ...", "     (", Sys.time(), ")")
  
  if(family == "binomial")
  {
    # Screen M using X given the limited information provided by Y (binary)
    # Therefore the family is still gaussian
    if(verbose) message("    Screening M using the association between X and M: ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                                   parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("    Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", 
                           parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)
  
  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  message("Step 2: Penalized estimate (", penalty, ") ...", "     (", Sys.time(), ")")
  
  ## Based on the screening results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    fit <- ncvreg(XM, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV.MY))), ...)
  }
  # plot(fit)
  
  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) message("    Tuning parameter lambda selected: ", lam)
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  if(length(ID_1_non) == 0)
  {
    if(verbose) message("    All ", penalty, " beta estimates of the ", length(ID), " mediators are zero. Please check ncvreg package for more options (?ncvreg).")
  } else {
    if(verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
    beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
    ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M
    ## 
    
    if(family == "binomial")
    {
      ## This has been done in step 1 (when Y is binary)
      alpha <- alpha[,ID_test, drop = FALSE]
    } else {
      if(verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
      alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", 
                       modelstatement = "Mone ~ X", parallel = FALSE, ncore = ncore, 
                       verbose, tag = "site-by-site ordinary least squares estimation")
    }
    
    #########################################################################
    ################################ STEP 3 #################################
    #########################################################################
    if(verbose) message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")
    
    alpha_est_ID_test <- as.numeric(alpha[1, ])  #  the estimator for alpha
    P_adjust_alpha <- alpha[2, ]  # the adjusted p-value for alpha (bonferroni)
    # P_adjust_alpha[P_adjust_alpha > 1] <- 1
    P_fdr_alpha <- p.adjust(alpha[2, ], "fdr")  # the adjusted p-value for alpha (FDR)
    
    alpha_est <- alpha_est_ID_test
    
    ## Post-test based on the oracle property of the MCP penalty
    if (is.null(COV.MY)) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
    }
    
    res <- summary(glm(Y ~ ., family = family, data = YMX))$coefficients
    est <- res[2:(length(ID_test) + 1), 1]  # the estimator for beta
    P_adjust_beta <- res[2:(length(ID_test) + 1), 4]  # the adjused p-value for beta (bonferroni)
    # P_adjust_beta[P_adjust_beta > 1] <- 1
    P_fdr_beta <- p.adjust(res[2:(length(ID_test) + 1), 4], "fdr")  # the adjusted p-value for beta (FDR)
    
    ab_est <- alpha_est * beta_est
    
    ## Use the maximum value as p value 
    PA <- rbind(P_adjust_beta, P_adjust_alpha)
    P_value <- apply(PA, 2, max)
    
    FDRA <- rbind(P_fdr_beta, P_fdr_alpha)
    FDR <- apply(FDRA, 2, max)
    
    # Total effect
    if (is.null(COV.MY)) {
      YX <- data.frame(Y = Y, X = X)
    } else {
      YX <- data.frame(Y = Y, X = X, COV.MY)
    }
    
    gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
    
    results <- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                          `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                          `adjusted.p` = P_value, `BH.FDR` = FDR, check.names = FALSE)
    
    message("Done!", "     (", Sys.time(), ")")
    
    doParallel::stopImplicitCluster()
    
    return(results)
  }
}
