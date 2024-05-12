library(MASS)

phfit_d = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a){
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
  lam = rep(1/l, l)
  p = ncol(X)
  beta = beta0 = rep(0, p)
  delta_g = delta_r = delta_g0 = delta_r0 = 0
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                      exp((tt>Tg[i])*delta_g+(tt>Tr[i])*delta_r)))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dd-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    dddelta_g = - sum(sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                         (tt>Tg[i])*exp(Xb[i]+delta_g+(tt>Tr[i])*delta_r))))
    dddelta_r = - sum(sapply(1:k, function(i) sum((tt<=Td[i])*lam*
                         (tt>Tr[i])*exp(Xb[i]+delta_r+(tt>Tg[i])*delta_g))))
    ddelta_g = dddelta_g + sum(Dd*(Tg<Td))
    ddelta_r = dddelta_r + sum(Dd*(Tr<Td))
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
    print(c(delta_g,delta_r,beta,tol))
    if (tol<0.0001) break
    beta0 = beta; delta_g0 = delta_g; delta_r0 = delta_r
  }
  return(list(beta=beta,delta_g=delta_g,delta_r=delta_r,tt=tt,lam=lam))
}

phfit_g = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a){
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
  lam = rep(1/l, l)
  p = ncol(X)
  beta = beta0 = rep(0, p)
  delta_r = delta_r0 = 0
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
    print(c(delta_r,beta,tol))
    if (tol<0.0001) break
    beta0 = beta; delta_r0 = delta_r
  }
  return(list(beta=beta,delta_r=delta_r,tt=tt,lam=lam))
}

phfit_r = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a){
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
  lam = rep(1/l, l)
  p = ncol(X)
  beta = beta0 = rep(0, p)
  delta_g = delta_g0 = 0
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
    print(c(delta_g,beta,tol))
    if (tol<0.0001) break
    beta0 = beta; delta_g0 = delta_g
  }
  return(list(beta=beta,delta_g=delta_g,tt=tt,lam=lam))
}

phfit_c = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = X[A==a,]
  Tc = Td
  Dc = 1 - Dd
  tt = sort(unique(Td[Dd==0]))
  l = length(tt)
  k = length(Tg)
  lam = rep(1/l, l)
  if (!is.null(X)) {
    p = ncol(X)
    beta = beta0 = rep(0, p)
  } else {
    beta = beta0 = 0
  }
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
    print(c(beta,tol))
    if (tol<0.0001) break
    beta0 = beta
  }
  return(list(beta=beta,tt=tt,lam=lam))
}

phfit_3 = function(Tr,Dr,Td,Dd,A,X,a){
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = X[A==a,]
  tt = sort(unique(Td[Dd==1]))
  l = length(tt)
  k = length(Td)
  lam = rep(1/l, l)
  p = ncol(X)
  beta = beta0 = rep(0, p)
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Td[i])*(tt>=Tr[i])*(Dr[i]==1)*lam))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dd*Dr-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    beta = beta - ginv(ddbeta) %*% dbeta
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dd*(Td==t)*Dr)/sum((Td>=t)*Dr*(Tr<=t)*
                                     exp(Xb)))
    lam[is.nan(lam)] = 0
    tol = max(abs(beta-beta0))
    if (tol<0.0001) break
    beta0 = beta
  }
  return(list(beta=beta,tt=tt,lam=lam))
}
