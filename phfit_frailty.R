library(statmod)
ngrid = 8
quad = gauss.quad(ngrid, "hermite")
b_list = quad$nodes
w_list = quad$weights

phfit = function(Tg,Dg,Tr,Dr,Td,Dd,A,X,a,par=NULL,se=FALSE){
  Tg = Tg[A==a]
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dg = Dg[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = as.matrix(as.matrix(X)[A==a,])
  tt = sort(unique(c(Td[Dd==1],Tg[Dg==1],Tr[Dr==1])))
  l = length(tt)
  k = length(Td)
  p = ncol(X)
  betad = betag = betar = rep(0, p)
  delta_gd = delta_rd = delta_gr = delta_rg = 0
  lamd = lamg = lamr = rep(1/l,l); sigma = 0.1
  est = c(betad,betag,betar,delta_gd,delta_rd,delta_gr,delta_rg,sigma)
  iter = 0
  while(iter<1000){
    iter = iter+1; est0 = est
    Xbd = as.numeric(X%*%betad)
    Xbg = as.numeric(X%*%betag)
    Xbr = as.numeric(X%*%betar)
    Ee0 = Ee1 = Ee2 = 0
    for (b in 1:ngrid){
      e = sqrt(2*sigma)*b_list[b]
      Lamd = sapply(1:k, function(i) sum(lamd*(tt<=Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
      dlamd = sapply(1:k, function(i) sum(lamd*(tt==Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
      Lamg = sapply(1:k, function(i) sum(lamg*(tt<=Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
      dlamg = sapply(1:k, function(i) sum(lamg*(tt==Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
      Lamr = sapply(1:k, function(i) sum(lamr*(tt<=Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
      dlamr = sapply(1:k, function(i) sum(lamr*(tt==Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
      dLamd = Dd*log(dlamd); dLamd[dlamd==0] = 0
      dLamg = Dg*log(dlamg); dLamg[dlamg==0] = 0
      dLamr = Dr*log(dlamr); dLamr[dlamr==0] = 0
      lik = exp(-Lamd+dLamd-Lamg+dLamg-Lamr+dLamr)*w_list[b]
      Ee0 = Ee0 + lik
      Ee1 = Ee1 + lik*e*e
      Ee2 = Ee2 + lik*exp(e)
    }
    sigma = min(mean(Ee1/Ee0),40)
    Eee = Ee2/Ee0
    # d
    S0 = sapply(Td, function(l) sum((Td>=l)*Eee*exp(Xbd+delta_gd*(Tg<l)+delta_rd*(Tr<l))))
    S1 = t(sapply(Td, function(l) colSums((Td>=l)*Eee*exp(Xbd+delta_gd*(Tg<l)+delta_rd*(Tr<l))*X)))
    S2 = t(sapply(Td, function(l) t(X)%*%diag((Td>=l)*Eee*exp(Xbd+delta_gd*(Tg<l)+delta_rd*(Tr<l)))%*%X))
    dbetad = as.numeric(t(X-S1/S0)%*%Dd)
    ddbetad = t(S1/S0)%*%diag(Dd)%*%(S1/S0) - matrix(colSums(Dd*S2/S0),p,p)
    S1 = sapply(Td, function(l) sum((Td>=l)*Eee*exp(Xbd+delta_gd*(Tg<l)+delta_rd*(Tr<l))*(Tg<l)))
    ddelta_gd = sum(Dd*((Td>Tg)*Dg-S1/S0))
    dddelta_gd = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0)
    S1 = sapply(Td, function(l) sum((Td>=l)*Eee*exp(Xbd+delta_gd*(Tg<l)+delta_rd*(Tr<l))*(Tr<l)))
    ddelta_rd = sum(Dd*((Td>Tr)*Dr-S1/S0))
    dddelta_rd = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0)
    betad = betad - ginv(ddbetad) %*% dbetad
    # g
    S0 = sapply(Tg, function(l) sum((Tg>=l)*Eee*exp(Xbg+delta_rg*(Tr<l))))
    S1 = t(sapply(Tg, function(l) colSums((Tg>=l)*Eee*exp(Xbg+delta_rg*(Tr<l))*X)))
    S2 = t(sapply(Tg, function(l) t(X)%*%diag((Tg>=l)*Eee*exp(Xbg+delta_rg*(Tr<l)))%*%X))
    dbetag = as.numeric(t(X-S1/S0)%*%Dg)
    ddbetag = t(S1/S0)%*%diag(Dg)%*%(S1/S0) - matrix(colSums(Dg*S2/S0),p,p)
    S1 = sapply(Tg, function(l) sum((Tg>=l)*Eee*exp(Xbg+delta_rg*(Tr<l))*(Tr<l)))
    ddelta_rg = sum(Dg*((Tg>Tr)*Dr-S1/S0))
    dddelta_rg = sum(Dg*(S1/S0)^2) - sum(Dg*S1/S0)
    betag = betag - ginv(ddbetag) %*% dbetag
    # r
    S0 = sapply(Tr, function(l) sum((Tr>=l)*Eee*exp(Xbr+delta_gr*(Tg<l))))
    S1 = t(sapply(Tr, function(l) colSums((Tr>=l)*Eee*exp(Xbr+delta_gr*(Tg<l))*X)))
    S2 = t(sapply(Tr, function(l) t(X)%*%diag((Tr>=l)*Eee*exp(Xbr+delta_gr*(Tg<l)))%*%X))
    dbetar = as.numeric(t(X-S1/S0)%*%Dr)
    ddbetar = t(S1/S0)%*%diag(Dr)%*%(S1/S0) - matrix(colSums(Dr*S2/S0),p,p)
    S1 = sapply(Tr, function(l) sum((Tr>=l)*Eee*exp(Xbr+delta_gr*(Tg<l))*(Tg<l)))
    ddelta_gr = sum(Dr*((Tr>Tg)*Dg-S1/S0))
    dddelta_gr = sum(Dr*(S1/S0)^2) - sum(Dr*S1/S0)
    betar = betar - ginv(ddbetar) %*% dbetar
    if (sum(Dd*(Tg<Td))==0) {
      delta_gd = 0
    } else {
      delta_gd = delta_gd - ddelta_gd/dddelta_gd
    }
    if (sum(Dd*(Tr<Td))==0) {
      delta_rd = 0
    } else {
      delta_rd = delta_rd - ddelta_rd/dddelta_rd
    }
    if (sum(Dg*(Tr<Tg))==0) {
      delta_rg = 0
    } else {
      delta_rg = delta_rg - ddelta_rg/dddelta_rg
    }
    if (sum(Dr*(Tg<Tr))==0) {
      delta_gr = 0
    } else {
      delta_gr = delta_gr - ddelta_gr/dddelta_gr
    }
    delta_gr = sign(delta_gr)*min(abs(delta_gr),4.5)
    delta_gd = sign(delta_gd)*min(abs(delta_gd),4.5)
    delta_rd = sign(delta_rd)*min(abs(delta_rd),4.5)
    delta_rg = sign(delta_rg)*min(abs(delta_rg),4.5)
    Xbd = as.numeric(X%*%betad)
    Xbg = as.numeric(X%*%betag)
    Xbr = as.numeric(X%*%betar)
    lamd = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*
                     Eee*exp(Xbd+(Tg<t)*delta_gd+(Tr<t)*delta_rd)))
    lamd[is.nan(lamd)] = 0
    lamg = sapply(tt, function(t) sum(Dg*(Tg==t))/sum((Tg>=t)*
                     Eee*exp(Xbg+(Tr<t)*delta_rg)))
    lamg[is.nan(lamg)] = 0
    lamr = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*
                     Eee*exp(Xbr+(Tg<t)*delta_gr)))
    lamr[is.nan(lamr)] = 0
    est = c(betad,betag,betar,delta_gd,delta_rd,delta_gr,delta_rg,sigma)
    tol = max(abs(est-est0))
    #print(c(delta_gd,delta_rd,sigma,tol))
    if (tol<0.0001) break
  }
  if (se){
    lamd0 = lamd
    lamg0 = lamg
    lamr0 = lamr
    Xbd = as.numeric(X%*%betad)
    Xbg = as.numeric(X%*%betag)
    Xbr = as.numeric(X%*%betar)
    Ee0 = 0
    for (b in 1:ngrid){
      e = sqrt(2*sigma)*b_list[b]
      Lamd = sapply(1:k, function(i) sum(lamd*(tt<=Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
      dlamd = sapply(1:k, function(i) sum(lamd*(tt==Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
      Lamg = sapply(1:k, function(i) sum(lamg*(tt<=Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
      dlamg = sapply(1:k, function(i) sum(lamg*(tt==Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
      Lamr = sapply(1:k, function(i) sum(lamr*(tt<=Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
      dlamr = sapply(1:k, function(i) sum(lamr*(tt==Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
      dLamd = Dd*log(dlamd); dLamd[dlamd==0] = 0
      dLamg = Dg*log(dlamg); dLamg[dlamg==0] = 0
      dLamr = Dr*log(dlamr); dLamr[dlamr==0] = 0
      lik = exp(-Lamd+dLamd-Lamg+dLamg-Lamr+dLamr)*w_list[b]
      Ee0 = Ee0 + lik
    }
    loglik0 = log(Ee0)
    loglik = NULL
    for (q in 1:length(est)){
      est1 = est
      est1[q] = est[q] + 1/sqrt(k)
      betad = est1[1:p]
      betag = est1[p+(1:p)]
      betar = est1[2*p+(1:p)]
      delta_gd = est1[3*p+1]
      delta_rd = est1[3*p+2]
      delta_gr = est1[3*p+3]
      delta_rg = est1[3*p+4]
      sigma = est1[3*p+5]
      Xbd = as.numeric(X%*%betad)
      Xbg = as.numeric(X%*%betag)
      Xbr = as.numeric(X%*%betar)
      lamd = lamd0
      lamg = lamg0
      lamr = lamr0
      for (iterh in 1:2){
      Ee0 = 0; Ee2 = 0
      for (b in 1:ngrid){
        e = sqrt(2*sigma)*b_list[b]
        lamd = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*
                          Eee*exp(Xbd+(Tg<t)*delta_gd+(Tr<t)*delta_rd)))
        lamd[is.nan(lamd)] = 0
        lamg = sapply(tt, function(t) sum(Dg*(Tg==t))/sum((Tg>=t)*
                          Eee*exp(Xbg+(Tr<t)*delta_rg)))
        lamg[is.nan(lamg)] = 0
        lamr = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*
                          Eee*exp(Xbr+(Tg<t)*delta_gr)))
        lamr[is.nan(lamr)] = 0
        Lamd = sapply(1:k, function(i) sum(lamd*(tt<=Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
        dlamd = sapply(1:k, function(i) sum(lamd*(tt==Td[i])*exp(e+Xbd[i]+delta_gd*(Tg[i]<tt)+delta_rd*(Tr[i]<tt))))
        Lamg = sapply(1:k, function(i) sum(lamg*(tt<=Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
        dlamg = sapply(1:k, function(i) sum(lamg*(tt==Tg[i])*exp(e+Xbg[i]+delta_rg*(Tr[i]<tt))))
        Lamr = sapply(1:k, function(i) sum(lamr*(tt<=Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
        dlamr = sapply(1:k, function(i) sum(lamr*(tt==Tr[i])*exp(e+Xbr[i]+delta_gr*(Tg[i]<tt))))
        dLamd = Dd*log(dlamd); dLamd[dlamd==0] = 0
        dLamg = Dg*log(dlamg); dLamg[dlamg==0] = 0
        dLamr = Dr*log(dlamr); dLamr[dlamr==0] = 0
        lik = exp(-Lamd+dLamd-Lamg+dLamg-Lamr+dLamr)*w_list[b]
        Ee0 = Ee0 + lik
        Ee2 = Ee2 + lik*exp(e)
      }
      Eee = Ee2/Ee0
      }
      loglik = cbind(loglik,log(Ee0))
    }
    dloglik = (loglik-loglik0)/(1/sqrt(k))
    V = ginv(t(dloglik)%*%dloglik)
    se = sqrt(diag(V))
    lamd = lamd0
    lamg = lamg0
    lamr = lamr0
  }
  betad = est[1:p]
  betag = est[p+(1:p)]
  betar = est[2*p+(1:p)]
  delta_gd = est[3*p+1]
  delta_rd = est[3*p+2]
  delta_gr = est[3*p+3]
  delta_rg = est[3*p+4]
  sigma = est[3*p+5]
  return(list(betad=betad,betag=betag,betar=betar,
              delta_gd=delta_gd,delta_rd=delta_rd,
              delta_gr=delta_gr,delta_rg=delta_rg,
              tt=tt,lamd=lamd,lamg=lamg,lamr=lamr,
              sigma=sigma,iter=iter,se=se))
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