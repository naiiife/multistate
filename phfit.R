library(MASS)

phfit_d = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = X[A==a,]
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
  Xb = as.numeric(X%*%beta)
  lam = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*
                       exp(Xb+(Tg<t)*delta_g+(Tr<t)*delta_r)))
  lam[is.nan(lam)] = 0
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                      exp((tt>Tg[i])*delta_g+(tt>Tr[i])*delta_r)))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dd-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    #S0 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))))
    #S1 = sapply(Td, function(l) colSums((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*X))
    #S2 = sapply(Td, function(l) t(X)%*%diag((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l)))%*%X)
    #S1 = t(S1)
    #S2 = t(S2)
    #dbeta = as.numeric(t(X-S1/S0)%*%Dd)
    #ddbeta = t(S1/S0)%*%diag(Dd)%*%(S1/S0) - matrix(colSums(Dd*S2/S0),p,p)
    dddelta_g = - sum(sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                         (tt>Tg[i])*exp(Xb[i]+delta_g+(tt>Tr[i])*delta_r))))
    dddelta_r = - sum(sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                         (tt>Tr[i])*exp(Xb[i]+delta_r+(tt>Tg[i])*delta_g))))
    ddelta_g = dddelta_g + sum(Dd*(Tg<Td))
    ddelta_r = dddelta_r + sum(Dd*(Tr<Td))
    #S1 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*Dg*(Tg<l)))
    #ddelta_g = sum(Dd*(Dg-S1/S0))
    #dddelta_g = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0^2)
    #S1 = sapply(Td, function(l) sum((Td>=l)*exp(Xb+delta_g*Dg*(Tg<l)+delta_r*Dr*(Tr<l))*Dr*(Tr<l)))
    #ddelta_r = sum(Dd*(Dr-S1/S0))
    #dddelta_r = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0^2)
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
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*
                                exp(Xb+(Tg<t)*delta_g+(Tr<t)*delta_r)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_g-delta_g0,delta_r-delta_r0)))
    #print(c(delta_g,delta_r,beta,tol))
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
  X = X[A==a,]
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
  Xb = as.numeric(X%*%beta)
  lam = sapply(tt, function(t) sum(Dg*(Tg==t))/sum((Tg>=t)*
                       exp(Xb+(Tr<t)*delta_r)))
  lam[is.nan(lam)] = 0
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Tg[i])*lam*exp((tt>Tr[i])*delta_r)))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dg-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    dddelta_r = - sum(sapply(1:k, function(i) sum((tt<=Tg[i])*lam*
                      (tt>Tr[i])*exp(Xb[i]+delta_r))))
    ddelta_r = dddelta_r + sum(Dg*(Tr<Tg))
    beta = beta - ginv(ddbeta) %*% dbeta
    if (sum(Dg*(Tr<Tg))==0) {
      delta_r = 0
    } else {
      delta_r = delta_r - ddelta_r/dddelta_r
    }
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dg*(Tg==t))/sum((Tg>=t)*
                               exp(Xb+(Tr<t)*delta_r)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_r-delta_r0)))
    #print(c(delta_r,beta,tol))
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
  X = X[A==a,]
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
  Xb = as.numeric(X%*%beta)
  lam = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*
                               exp(Xb+(Tg<t)*delta_g)))
  lam[is.nan(lam)] = 0
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Tr[i])*lam*exp((tt>Tg[i])*delta_g)))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dr-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    dddelta_g = - sum(sapply(1:k, function(i) sum((tt<=Tr[i])*lam*
                         (tt>Tg[i])*exp(Xb[i]+delta_g))))
    ddelta_g = dddelta_g + sum(Dr*(Tg<Tr))
    beta = beta - ginv(ddbeta) %*% dbeta
    if (sum(Dr*(Tg<Tr))==0) {
      delta_g = 0
    } else {
      delta_g = delta_g - ddelta_g/dddelta_g
    }
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*
                                                       exp(Xb+(Tg<t)*delta_g)))
    lam[is.nan(lam)] = 0
    tol = max(abs(c(beta-beta0,delta_g-delta_g0)))
    #print(c(delta_g,beta,tol))
    if (tol<0.00001) break
    beta0 = beta; delta_g0 = delta_g
  }
  return(list(beta=beta,delta_g=delta_g,tt=tt,lam=lam))
}

phfit_c = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL){
  Td = Td[A==a]
  Dd = Dd[A==a]
  X = X[A==a,]
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
  Xb = as.numeric(X%*%beta)
  lam = sapply(tt, function(t) sum(Dc*(Tc==t))/sum((Tc>=t)*exp(Xb)))
  lam[is.nan(lam)] = 0
  while(TRUE){
    if (!is.null(X)) {
      Xb = as.numeric(X%*%beta)
    } else {
      Xb = 1
    }
    Lam = sapply(1:k, function(i) sum((tt<=Tc[i])*lam)) * exp(Xb)
    if (!is.null(X)) {
    dbeta = t(X)%*%(Dc-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
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