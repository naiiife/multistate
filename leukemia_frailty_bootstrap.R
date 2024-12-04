set.seed(2024)
cif1 = cif0 = cif1o = cif0o = cif_spe = cif_cin = cif_int = NULL
sig1 = sig0 = NULL

for (bs in 1:B){
  cat('Bootstrap round',bs,'of',B,'\n')
  ss1 = sample(which(dat$A==1),replace=TRUE)
  ss0 = sample(which(dat$A==0),replace=TRUE)
  ss = c(ss1,ss0)
  A = dat$A[ss]
  X = as.matrix(dat$X[ss,])
  Tg=dat$Tg[ss];Dg=dat$Dg[ss];Tr=dat$Tr[ss]
  Dr=dat$Dr[ss];Td=dat$Td[ss];Dd=dat$Dd[ss]
  
  fit1 = try(phfit(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1))
  fit0 = try(phfit(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0))
  if ('try-error' %in% class(fit1)) next
  if ('try-error' %in% class(fit0)) next
  sig1 = append(sig1,fit1$sigma)
  sig0 = append(sig0,fit0$sigma)
  
  F_1.b = F_0.b = F_1o.b = F_0o.b = Fd_spe.b = Fd_int.b = Fd_cin.b = 0
  for (b in 1:ngrid){
    e = sqrt(2*fit1$sigma)*b_list[b]
    Xb = as.numeric(X%*%fit1$betad)+e
    delta_g = fit1$delta_gd
    delta_r = fit1$delta_rd
    lam_d = matchy(fit1$tt, fit1$lamd, tt)
    lam_od1 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
    lam_ogd1 = lam_od1 * exp(delta_g)
    lam_ord1 = lam_od1 * exp(delta_r)
    lam_ogrd1 = lam_orgd1 = lam_od1 * exp(delta_g+delta_r)
    Xb = as.numeric(X%*%fit1$betag)+e
    delta_r = fit1$delta_rg
    lam_g = matchy(fit1$tt, fit1$lamg, tt)
    lam_og1 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
    lam_org1 = lam_og1 *exp(delta_r)
    Xb = as.numeric(X%*%fit1$betar)+e
    delta_g = fit1$delta_gr
    lam_r = matchy(fit1$tt, fit1$lamr, tt)
    lam_or1 = sapply(1:l, function(t) lam_r[t]*exp(Xb))
    lam_ogr1 = lam_or1 * exp(delta_g)
    
    e = sqrt(2*fit0$sigma)*b_list[b]
    Xb = as.numeric(X%*%fit0$betad)+e
    delta_g = fit0$delta_gd
    delta_r = fit0$delta_rd
    lam_d = matchy(fit0$tt, fit0$lamd, tt)
    lam_od0 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
    lam_ogd0 = lam_od0 * exp(delta_g)
    lam_ord0 = lam_od0 * exp(delta_r)
    lam_ogrd0 = lam_orgd0 = lam_od0 * exp(delta_g+delta_r)
    Xb = as.numeric(X%*%fit0$betag)+e
    delta_r = fit0$delta_rg
    lam_g = matchy(fit0$tt, fit0$lamg, tt)
    lam_og0 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
    lam_org0 = lam_og0 *exp(delta_r)
    Xb = as.numeric(X%*%fit0$betar)+e
    delta_g = fit0$delta_gr
    lam_r = matchy(fit0$tt, fit0$lamr, tt)
    lam_or0 = sapply(1:l, function(t) lam_r[t]*exp(Xb))
    lam_ogr0 = lam_or0 * exp(delta_g)
    
    # observable incidence
    lam_og_A = A*lam_og1 + (1-A)*lam_og0
    lam_od_A = A*lam_od1 + (1-A)*lam_od0
    lam_or_A = A*lam_or1 + (1-A)*lam_or0
    lam_ogr_A = A*lam_ogr1 + (1-A)*lam_ogr0
    lam_ogd_A = A*lam_ogd1 + (1-A)*lam_ogd0
    lam_org_A = A*lam_org1 + (1-A)*lam_org0
    lam_ord_A = A*lam_ord1 + (1-A)*lam_ord0
    lam_ogrd_A = A*lam_ogrd1 + (1-A)*lam_ogrd0
    lam_orgd_A = A*lam_orgd1 + (1-A)*lam_orgd0
    Lam_o_A = t(apply(lam_od_A+lam_og_A+lam_or_A, 1, cumsum))
    Lam_og_A = t(apply(lam_ogd_A+lam_ogr_A, 1, cumsum))
    Lam_or_A = t(apply(lam_ord_A+lam_org_A, 1, cumsum))
    Lam_ogr_A = t(apply(lam_ogrd_A, 1, cumsum))
    Lam_org_A = t(apply(lam_orgd_A, 1, cumsum))
    
    dF_or_A = exp(-Lam_o_A)*lam_or_A
    dF_og_A = exp(-Lam_o_A)*lam_og_A
    dF_od_A = exp(-Lam_o_A)*lam_od_A
    F_or_A = t(apply(dF_or_A, 1, cumsum))
    F_og_A = t(apply(dF_og_A, 1, cumsum))
    F_od_A = t(apply(dF_od_A, 1, cumsum))
    dF_og__A = t(apply(dF_og_A*exp(Lam_og_A), 1, cumsum))*exp(-Lam_og_A)
    dF_or__A = t(apply(dF_or_A*exp(Lam_or_A), 1, cumsum))*exp(-Lam_or_A)
    dF_ogr_A = dF_og__A*lam_ogr_A
    F_ogr_A = t(apply(dF_ogr_A, 1, cumsum))
    dF_org_A = dF_or__A*lam_org_A
    F_org_A = t(apply(dF_org_A, 1, cumsum))
    dF_ord_A = dF_or__A*lam_ord_A
    F_ord_A = t(apply(dF_ord_A, 1, cumsum))
    dF_ogd_A = dF_og__A*lam_ogd_A
    F_ogd_A = t(apply(dF_ogd_A, 1, cumsum))
    dF_org__A = t(apply(dF_org_A*exp(Lam_org_A), 1, cumsum))*exp(-Lam_org_A)
    dF_ogr__A = t(apply(dF_ogr_A*exp(Lam_ogr_A), 1, cumsum))*exp(-Lam_ogr_A)
    dF_orgd_A = dF_org__A*lam_orgd_A
    F_orgd_A = t(apply(dF_orgd_A, 1, cumsum))
    dF_ogrd_A = dF_ogr__A*lam_ogrd_A
    F_ogrd_A = t(apply(dF_ogrd_A, 1, cumsum))
    F_A = colMeans(F_od_A+F_ord_A+F_ogd_A+F_orgd_A+F_ogrd_A)
    
    # potential incidence
    #prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
    #lam_or1 = lam_ogr1 = prev_g*lam_ogr1+(1-prev_g)*lam_or1
    Lam_o_1 = t(apply(lam_od1+lam_og1+lam_or1, 1, cumsum))
    Lam_og_1 = t(apply(lam_ogd1+lam_ogr1, 1, cumsum))
    Lam_or_1 = t(apply(lam_ord1+lam_org1, 1, cumsum))
    Lam_ogr_1 = t(apply(lam_ogrd1, 1, cumsum))
    Lam_org_1 = t(apply(lam_orgd1, 1, cumsum))
    dF_or_1 = exp(-Lam_o_1)*lam_or1
    dF_og_1 = exp(-Lam_o_1)*lam_og1
    dF_od_1 = exp(-Lam_o_1)*lam_od1
    F_or_1 = t(apply(dF_or_1, 1, cumsum))
    F_og_1 = t(apply(dF_og_1, 1, cumsum))
    F_od_1 = t(apply(dF_od_1, 1, cumsum))
    dF_og__1 = t(apply(dF_og_1*exp(Lam_og_1), 1, cumsum))*exp(-Lam_og_1)
    dF_or__1 = t(apply(dF_or_1*exp(Lam_or_1), 1, cumsum))*exp(-Lam_or_1)
    dF_ogr_1 = dF_og__1*lam_ogr1
    F_ogr_1 = t(apply(dF_ogr_1, 1, cumsum))
    dF_org_1 = dF_or__1*lam_org1
    F_org_1 = t(apply(dF_org_1, 1, cumsum))
    dF_ord_1 = dF_or__1*lam_ord1
    F_ord_1 = t(apply(dF_ord_1, 1, cumsum))
    dF_ogd_1 = dF_og__1*lam_ogd1
    F_ogd_1 = t(apply(dF_ogd_1, 1, cumsum))
    dF_org__1 = t(apply(dF_org_1*exp(Lam_org_1), 1, cumsum))*exp(-Lam_org_1)
    dF_ogr__1 = t(apply(dF_ogr_1*exp(Lam_ogr_1), 1, cumsum))*exp(-Lam_ogr_1)
    dF_orgd_1 = dF_org__1*lam_orgd1
    F_orgd_1 = t(apply(dF_orgd_1, 1, cumsum))
    dF_ogrd_1 = dF_ogr__1*lam_ogrd1
    F_ogrd_1 = t(apply(dF_ogrd_1, 1, cumsum))
    F_1.b = F_1.b + colMeans(F_od_1+F_ord_1+F_ogd_1+F_orgd_1+F_ogrd_1)*w_list[b]/sqrt(pi)
    
    Lam_o_0 = t(apply(lam_od0+lam_og0+lam_or0, 1, cumsum))
    Lam_og_0 = t(apply(lam_ogd0+lam_ogr0, 1, cumsum))
    Lam_or_0 = t(apply(lam_ord0+lam_org0, 1, cumsum))
    Lam_ogr_0 = t(apply(lam_ogrd0, 1, cumsum))
    Lam_org_0 = t(apply(lam_orgd0, 1, cumsum))
    dF_or_0 = exp(-Lam_o_0)*lam_or0
    dF_og_0 = exp(-Lam_o_0)*lam_og0
    dF_od_0 = exp(-Lam_o_0)*lam_od0
    F_or_0 = t(apply(dF_or_0, 1, cumsum))
    F_og_0 = t(apply(dF_og_0, 1, cumsum))
    F_od_0 = t(apply(dF_od_0, 1, cumsum))
    dF_og__0 = t(apply(dF_og_0*exp(Lam_og_0), 1, cumsum))*exp(-Lam_og_0)
    dF_or__0 = t(apply(dF_or_0*exp(Lam_or_0), 1, cumsum))*exp(-Lam_or_0)
    dF_ogr_0 = dF_og__0*lam_ogr0
    F_ogr_0 = t(apply(dF_ogr_0, 1, cumsum))
    dF_org_0 = dF_or__0*lam_org0
    F_org_0 = t(apply(dF_org_0, 1, cumsum))
    dF_ord_0 = dF_or__0*lam_ord0
    F_ord_0 = t(apply(dF_ord_0, 1, cumsum))
    dF_ogd_0 = dF_og__0*lam_ogd0
    F_ogd_0 = t(apply(dF_ogd_0, 1, cumsum))
    dF_org__0 = t(apply(dF_org_0*exp(Lam_org_0), 1, cumsum))*exp(-Lam_org_0)
    dF_ogr__0 = t(apply(dF_ogr_0*exp(Lam_ogr_0), 1, cumsum))*exp(-Lam_ogr_0)
    dF_orgd_0 = dF_org__0*lam_orgd0
    F_orgd_0 = t(apply(dF_orgd_0, 1, cumsum))
    dF_ogrd_0 = dF_ogr__0*lam_ogrd0
    F_ogrd_0 = t(apply(dF_ogrd_0, 1, cumsum))
    F_0.b = F_0.b + colMeans(F_od_0+F_ord_0+F_ogd_0+F_orgd_0+F_ogrd_0)*w_list[b]/sqrt(pi)
    
    # counterfactual incidence, separable
    # a_or=a_gr=a_og=a_rg=0,a_od=a_rd=a_gd=1
    lam_og_a = lam_og0
    lam_od_a = lam_od1
    lam_or_a = lam_or0
    lam_ogr_a = lam_ogr0
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org0
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_spe.b = Fd_spe.b + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # counterfactual incidence, conditionally interventional
    # a_or=a_gr=0,a_od=a_og=a_rd=a_gd=a_rg=1
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_or_a = lam_or0
    lam_ogr_a = lam_ogr0
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_cin.b = Fd_cin.b + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # potential incidence, overall
    # lam_or = lam_ogr
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr1+(1-prev_g)*lam_or1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    F_1o.b = F_1o.b + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # lam_or = lam_ogr
    lam_og_a = lam_og0
    lam_od_a = lam_od0
    lam_ogd_a = lam_ogd0
    lam_org_a = lam_org0
    lam_ord_a = lam_ord0
    lam_ogrd_a = lam_ogrd0
    lam_orgd_a = lam_orgd0
    prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    F_0o.b = F_0o.b + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # counterfactual incidence, interventional
    # lam_or = lam_ogr
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_int.b = Fd_int.b + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
  }
  
  cif1 = rbind(cif1, F_1.b)
  cif0 = rbind(cif0, F_0.b)
  cif1o = rbind(cif1o, F_1o.b)
  cif0o = rbind(cif0o, F_0o.b)
  cif_spe = rbind(cif_spe, Fd_spe.b)
  cif_cin = rbind(cif_cin, Fd_cin.b)
  cif_int = rbind(cif_int, Fd_int.b)
}
