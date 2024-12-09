## NPMLE

phfit_d = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = as.matrix(as.matrix(X)[A==a,])
  tt = sort(unique(Td[Dd==1]))
  l = length(tt)
  k = length(Td)
  p = ncol(X)
  if (!is.null(par)){
    beta = beta0 = par[2+1:p]
    delta_g = delta_g0 = par[1]
    delta_r = delta_r0 = par[2]
  } else {
    beta = beta0 = rep(0, p)
    delta_g = delta_r = delta_g0 = delta_r0 = 0
  }
  lam = rep(1/l,l)
  iter = 0
  while(iter<100){
    iter = iter+1
    Xb = as.numeric(X%*%beta)
    S0 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))))
    S1 = t(sapply(Td, function(l) colSums((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*X)))
    S2 = t(sapply(Td, function(l) t(X)%*%diag((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l)))%*%X))
    dbeta = as.numeric(t(X-S1/S0)%*%Dd)
    ddbeta = t(S1/S0)%*%diag(Dd)%*%(S1/S0) - matrix(colSums(Dd*S2/S0),p,p)
    S1 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*Dg*(Tg<l)))
    ddelta_g = sum(Dd*((Td>Tg)*Dg-S1/S0))
    dddelta_g = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0)
    S1 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*Dr*(Tr<l)))
    ddelta_r = sum(Dd*((Td>Tr)*Dr-S1/S0))
    dddelta_r = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0)
    beta = beta - ginv(ddbeta) %*% dbeta
    if (sum(Dd*(Tg<Td))==0) {
      delta_g = 0
    } else {
      delta_g = delta_g - ddelta_g/dddelta_g
    }
    if (sum(Dd*(Tr<Td))==0) {
      delta_r = 0
    } else {
      delta_r = delta_r - ddelta_r/dddelta_r
    }
    delta_g = sign(delta_g)*min(abs(delta_g),4.2)
    delta_r = sign(delta_r)*min(abs(delta_r),4.2)
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*
         exp(Xb+(Tg<t)*delta_g+(Tr<t)*delta_r)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_g-delta_g0,delta_r-delta_r0)))
    #print(c(beta,delta_g,delta_r,tol))
    if (tol<0.00001) break
    beta0 = beta; delta_g0 = delta_g; delta_r0 = delta_r
  }
  return(list(beta=beta,delta_g=delta_g,delta_r=delta_r,tt=tt,lam=lam))
}

phfit_g = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = as.matrix(as.matrix(X)[A==a,])
  tt = sort(unique(Tg[Dg==1]))
  l = length(tt)
  k = length(Tg)
  p = ncol(X)
  if (!is.null(par)){
    beta = beta0 = par[1+1:p]
    delta_r = delta_r0 = par[1]
  } else {
    beta = beta0 = rep(0, p)
    delta_r = delta_r0 = 0
  }
  lam = rep(1/l,l)
  iter = 0
  while(iter<100){
    iter = iter+1
    Xb = as.numeric(X%*%beta)
    S0 = sapply(Tg, function(l) sum((Tg>=l)*exp(Xb+delta_r*Dr*(Tr<l))))
    S1 = t(sapply(Tg, function(l) colSums((Tg>=l)*exp(Xb+delta_r*Dr*(Tr<l))*X)))
    S2 = t(sapply(Tg, function(l) t(X)%*%diag((Tg>=l)*exp(Xb+delta_r*Dr*(Tr<l)))%*%X))
    dbeta = as.numeric(t(X-S1/S0)%*%Dg)
    ddbeta = t(S1/S0)%*%diag(Dg)%*%(S1/S0) - matrix(colSums(Dg*S2/S0),p,p)
    S1 = sapply(Tg, function(l) sum((Tg>=l)*exp(Xb+delta_r*Dr*(Tr<l))*Dr*(Tr<l)))
    ddelta_r = sum(Dg*((Tg>Tr)*Dr-S1/S0))
    dddelta_r = sum(Dg*(S1/S0)^2) - sum(Dg*S1/S0)
    beta = beta - ginv(ddbeta) %*% dbeta
    if (sum(Dg*(Tr<Tg))==0) {
      delta_r = 0
    } else {
      delta_r = delta_r - ddelta_r/dddelta_r
    }
    delta_r = sign(delta_r)*min(abs(delta_r),4.2)
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dg*(Tg==t))/sum((Tg>=t)*
                                                       exp(Xb+(Tr<t)*delta_r)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_r-delta_r0)))
    #print(c(beta,delta_r,tol))
    if (tol<0.00001) break
    beta0 = beta; delta_r0 = delta_r
  }
  return(list(beta=beta,delta_r=delta_r,tt=tt,lam=lam))
}

phfit_r = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = as.matrix(as.matrix(X)[A==a,])
  tt = sort(unique(Tr[Dr==1]))
  l = length(tt)
  k = length(Tr)
  p = ncol(X)
  if (!is.null(par)){
    beta = beta0 = par[1+1:p]
    delta_g = delta_g0 = par[1]
  } else {
    beta = beta0 = rep(0, p)
    delta_g = delta_g0 = 0
  }
  lam = rep(1/l,l)
  iter = 0
  while(iter<100){
    iter = iter+1
    Xb = as.numeric(X%*%beta)
    S0 = sapply(Tr, function(l) sum((Tr>=l)*exp(Xb+delta_g*Dg*(Tg<l))))
    S1 = t(sapply(Tr, function(l) colSums((Tr>=l)*exp(Xb+delta_g*Dg*(Tg<l))*X)))
    S2 = t(sapply(Tr, function(l) t(X)%*%diag((Tr>=l)*exp(Xb+delta_g*Dg*(Tg<l)))%*%X))
    dbeta = as.numeric(t(X-S1/S0)%*%Dr)
    ddbeta = t(S1/S0)%*%diag(Dr)%*%(S1/S0) - matrix(colSums(Dr*S2/S0),p,p)
    S1 = sapply(Tr, function(l) sum((Tr>=l)*exp(Xb+delta_g*Dg*(Tg<l))*Dg*(Tg<l)))
    ddelta_g = sum(Dr*((Tr>Tg)*Dg-S1/S0))
    dddelta_g = sum(Dr*(S1/S0)^2) - sum(Dr*S1/S0)
    beta = beta - ginv(ddbeta) %*% dbeta
    if (sum(Dr*(Tg<Tr))==0) {
      delta_g = 0
    } else {
      delta_g = delta_g - ddelta_g/dddelta_g
    }
    delta_g = sign(delta_g)*min(abs(delta_g),4.2)
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*
                                                       exp(Xb+(Tg<t)*delta_g)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_g-delta_g0)))
    #print(c(beta,delta_g,tol))
    if (tol<0.00001) break
    beta0 = beta; delta_g0 = delta_g
  }
  return(list(beta=beta,delta_g=delta_g,tt=tt,lam=lam))
}

phfit_c = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Td = Td[A==a]
  Dd = Dd[A==a]
  X = as.matrix(as.matrix(X)[A==a,])
  Tc = Td
  Dc = 1 - Dd
  tt = sort(unique(Td[Dd==0]))
  l = length(tt)
  k = length(Tc)
  if (!is.null(X)) {
    p = ncol(X)
    beta = beta0 = rep(0, p)
  } else {
    beta = beta0 = 0
  }
  if (!is.null(par)){
    beta = beta0 = par
  }
  lam = rep(1/l,l)
  iter = 0
  while(iter<100){
    iter = iter+1
    if (!is.null(X)) {
      Xb = as.numeric(X%*%beta)
    } else {
      Xb = 1
    }
    if (!is.null(X)) {
      S0 = sapply(Tc, function(l) sum((Tc>=l)*exp(Xb)))
      S1 = t(sapply(Tc, function(l) colSums((Tc>=l)*exp(Xb)*X)))
      S2 = t(sapply(Tc, function(l) t(X)%*%diag((Tc>=l)*exp(Xb))%*%X))
      dbeta = as.numeric(t(X-S1/S0)%*%Dc)
      ddbeta = t(S1/S0)%*%diag(Dc)%*%(S1/S0) - matrix(colSums(Dc*S2/S0),p,p)
      beta = beta - ginv(ddbeta) %*% dbeta
      Xb = as.numeric(X%*%beta)
    }
    lam = sapply(tt, function(t) sum(Dc*(Tc==t))/sum((Tc>=t)*exp(Xb)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0)))
    #print(c(beta,tol))
    if (tol<0.00001) break
    beta0 = beta
  }
  return(list(beta=beta,tt=tt,lam=lam))
}


testcif <- function(fit1,fit0,maxt=NULL){
  tt = fit1$tt
  F1 = colMeans(fit1$Feff, na.rm=TRUE)
  F0 = colMeans(fit0$Feff, na.rm=TRUE)
  IF1 = fit1$EIF
  IF0 = fit0$EIF
  Ti = rep(TRUE, length(tt))
  if (!is.null(maxt)) Ti = (tt<=maxt)
  Tt = sum((F1*(1-F1)-F0*(1-F0))*diff(c(0,tt))*Ti)
  IFt = colSums(((1-2*F1)*t(IF1)-(1-2*F0)*t(IF0))*diff(c(0,tt))*Ti)
  Vt = sd(IFt, na.rm=TRUE)/sqrt(nrow(IF1))
  p1 = 2*pnorm(-abs(Tt/Vt))
  Tt = sum((F1-F0)*diff(c(0,F1+F0))*Ti)
  V1 = colSums(t(IF1-IF0)*diff(c(0,F1+F0))*Ti)
  V2 = colSums((F1-F0)*apply(cbind(0,IF1+IF0),1,diff)*Ti)
  Vt = sd(V1+V2)/sqrt(nrow(IF1))
  p2 = 2*pnorm(-abs(Tt/Vt))
  return(list(p1=p1,p2=p2))
}

matchy = function(x,y,newx,exact=TRUE){
  if (exact) {
    newy = sapply(newx, function(t) y[max(which(x==t))])
  } else {
    newy = sapply(newx, function(t) y[max(which(x<=t))])
  }
  newy = as.numeric(newy)
  newy[is.na(newy)] = 0
  return(newy)
}