setwd('D:/Papers/Multistate/codes/')
source('phfit_npmle.R')
dat = read.csv('leukemia.csv')
dat = transform(dat, MRD=as.numeric(MRDPRET>0), SEX=SEX1-1,
                CR=DISEASESTATUS-1, DIAGNOSIS=TALLORBALL-1,
                TRANSPLANT=YZLX-1)
dat = dat[dat$MRD==1,]
A = dat$TRANSPLANT # 1 Haplo-SCT 0 MSDT
Tc = dat$AGVHDT
Dc = as.numeric(dat$AGVHD>0)
Tc[Dc==0] = dat$CGVHDT[Dc==0]
Dc[Dc==0] = as.numeric(dat$CGVHD[Dc==0]>0)
Tm = dat$RELAPSET
Dm = dat$RELAPSE
Td = dat$OST
Dd = dat$OS
#Tc[764] = Td[764]; Dc[764] = 0  #record error
Tc[Tc==Tm & Dc==1] = Tc[Tc==Tm & Dc==1] - 0.01
Td[Tc==Td & Dc==1] = Td[Tc==Td & Dc==1] + 0.01
Td[Tm==Td & Dm==1] = Td[Tm==Td & Dm==1] + 0.01
Tg = Tc; Dg = Dc
Tr = Tm; Dr = Dm 
n = nrow(dat)
dat$Tr=Tr;dat$Td=Td;dat$Tg=Tg
dat$Dr=Dr;dat$Dd=Dd;dat$Dg=Dg
dat$A = A

X = dat[,c('AGE','SEX','CR','DIAGNOSIS')]
dat$X = X
X = as.matrix(X)
X[,1] = (X[,1]-mean(X[,1]))/sd(X[,1])
dat$X = X
p = ncol(X)
tt = sort(unique(c(Td,Tr,Tg)))
l = length(tt)


# hazard of d
fit_d = phfit_d(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1,
                par=c(0.216,3.401,-0.114,0.249,0.263,0.015,0.216))
Xb = as.numeric(X%*%fit_d$beta)
delta_g = fit_d$delta_g
delta_r = fit_d$delta_r
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_od1 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
lam_ogd1 = lam_od1 * exp(delta_g)
lam_ord1 = lam_od1 * exp(delta_r)
lam_ogrd1 = lam_orgd1 = lam_od1 * exp(delta_g+delta_r)
fit_d = phfit_d(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0,
                par=c(1.147,4.087,-0.460,-0.014,0.209,-0.545))
Xb = as.numeric(X%*%fit_d$beta)
delta_g = fit_d$delta_g
delta_r = fit_d$delta_r
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_od0 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
lam_ogd0 = lam_od0 * exp(delta_g)
lam_ord0 = lam_od0 * exp(delta_r)
lam_ogrd0 = lam_orgd0 = lam_od0 * exp(delta_g+delta_r)
# hazard of g
fit_g = phfit_g(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1,
                par=c(-0.032,-0.021,-0.010,0.017,0.183))
Xb = as.numeric(X%*%fit_g$beta)
delta_r = fit_g$delta_r
lam_g = matchy(fit_g$tt, fit_g$lam, tt)
lam_og1 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
lam_org1 = lam_og1 *exp(delta_r)
fit_g = phfit_g(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0,
                par=c(0,0.299,-0.181,1.074,-0.896))
Xb = as.numeric(X%*%fit_g$beta)
delta_r = fit_g$delta_r
lam_g = matchy(fit_g$tt, fit_g$lam, tt)
lam_og0 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
lam_org0 = lam_og0 *exp(delta_r)
# hazard of r
fit_r = phfit_r(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1,
                par=c(-0.044,-0.021,-0.139,0.845,0.875))
Xb = as.numeric(X%*%fit_r$beta)
delta_g = fit_r$delta_g
lam_r = matchy(fit_r$tt, fit_r$lam, tt)
lam_or1 = sapply(1:l, function(t) lam_r[t]*exp(Xb))
lam_ogr1 = lam_or1 * exp(delta_g)
fit_r = phfit_r(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0,
                par=c(0.192,-0.205,-0.200,0.767,-0.029))
Xb = as.numeric(X%*%fit_r$beta)
delta_g = fit_r$delta_g
lam_r = matchy(fit_r$tt, fit_r$lam, tt)
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
F_1 = colMeans(F_od_1+F_ord_1+F_ogd_1+F_orgd_1+F_ogrd_1)

#prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
#lam_or0 = lam_ogr0 = prev_g*lam_ogr0+(1-prev_g)*lam_or0
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
F_0 = colMeans(F_od_0+F_ord_0+F_ogd_0+F_orgd_0+F_ogrd_0)

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
Fd_spe = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)

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
Fd_cin = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)

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
F_1o = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)

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
F_0o = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)

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
Fd_int = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)

B = 200
source('leukemia_bootstrap.R')

quantile(exp(delta1[,2]),c(.025,.975))
quantile(exp(delta1[,4]),c(.025,.975))
quantile(exp(delta1[,1]),c(.025,.975))

index = which(tt<=2500)
tti = tt[index]
F1 = F_1[index]
F0 = F_0[index]
F1o = F_1o[index]
F0o = F_0o[index]
F2 = Fd_spe[index]
F3 = Fd_int[index]
F4 = Fd_cin[index]
sde = F2 - F0
sie = F1 - F2
ide = F3 - F0o
iie = F1o - F3
cde = F4 - F0
cie = F1 - F4
se1 = apply(cif1,2,sd,na.rm=TRUE)[index]
se0 = apply(cif0,2,sd,na.rm=TRUE)[index]
se2 = apply(cif_spe,2,sd,na.rm=TRUE)[index]
se3 = apply(cif_int,2,sd,na.rm=TRUE)[index]
se4 = apply(cif_cin,2,sd,na.rm=TRUE)[index]
se.sde = apply(cif_spe-cif0,2,sd,na.rm=TRUE)[index]
se.sie = apply(cif1-cif_spe,2,sd,na.rm=TRUE)[index]
se.ide = apply(cif_int-cif0o,2,sd,na.rm=TRUE)[index]
se.iie = apply(cif1o-cif_int,2,sd,na.rm=TRUE)[index]
se.cde = apply(cif_cin-cif0,2,sd,na.rm=TRUE)[index]
se.cie = apply(cif1-cif_cin,2,sd,na.rm=TRUE)[index]
sde.l = apply(cif_spe-cif0,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
sde.u = apply(cif_spe-cif0,2,function(l) quantile(l,.975,na.rm=TRUE))[index]
sie.l = apply(cif1-cif_spe,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
sie.u = apply(cif1-cif_spe,2,function(l) quantile(l,.975,na.rm=TRUE))[index]
ide.l = apply(cif_int-cif0o,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
ide.u = apply(cif_int-cif0o,2,function(l) quantile(l,.975,na.rm=TRUE))[index]
iie.l = apply(cif1o-cif_int,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
iie.u = apply(cif1o-cif_int,2,function(l) quantile(l,.975,na.rm=TRUE))[index]
cde.l = apply(cif_cin-cif0,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
cde.u = apply(cif_cin-cif0,2,function(l) quantile(l,.975,na.rm=TRUE))[index]
cie.l = apply(cif1-cif_cin,2,function(l) quantile(l,.025,na.rm=TRUE))[index]
cie.u = apply(cif1-cif_cin,2,function(l) quantile(l,.975,na.rm=TRUE))[index]

# CIF plots
plot(tti, F1, type='s', lwd=1.5, lty=1, ylim=c(0,.8), col='brown',
     ylab='Cumulative incidence', xlab='Days since transplantation', 
     main='Cumulative incidence of death')
#points(tti, F1+1.96*se1, type='s', lty=2, col='brown')
#points(tti, F1-1.96*se1, type='s', lty=2, col='brown')
points(tti, F0, type='s', lwd=1.5, lty=1, col='purple')
#points(tti, F0+1.96*se0, type='s', lty=2, col='purple')
#points(tti, F0-1.96*se0, type='s', lty=2, col='purple')
points(tti, F1o, type='s', lwd=1.5, lty=4, col='brown')
points(tti, F0o, type='s', lwd=1.5, lty=4, col='purple')
points(tti, F2, type='s', lwd=1.5, col='darkorange')
#points(tti, F2+1.96*se2, type='s', lty=2, col='darkorange')
#points(tti, F2-1.96*se2, type='s', lty=2, col='darkorange')
points(tti, F3, type='s', lwd=1.5, col='darkblue')
#points(tti, F3+1.96*se3, type='s', lty=2, col='darkblue')
#points(tti, F3-1.96*se3, type='s', lty=2, col='darkblue')
points(tti, F4, type='s', lwd=1.5, col='darkgreen')
#points(tti, F4+1.96*se4, type='s', lty=2, col='darkblue')
#points(tti, F4-1.96*se4, type='s', lty=2, col='darkblue')
legend('topleft', cex=0.9, lty=c(1,1,1,1,1), lwd=rep(2,5), 
       col=c('brown','purple','darkorange','darkblue','darkgreen'),
       legend=c('Haplo-SCT','MSDT','Separable','Interventional','Cond. Inter.'))

plot(tti, sde, type='s', lwd=1.5, ylim=c(-.4,.4), col='darkorange',
     ylab='Diff in cumulative incidence', xlab='Days since transplantation', main='Direct effect')
abline(h=0, lty=5)
points(tti, sde+1.96*se.sde, type='s', lty=2, col='darkorange')
points(tti, sde-1.96*se.sde, type='s', lty=2, col='darkorange')
#points(tti, sde.l, type='s', lty=2, col='darkorange')
#points(tti, sde.u, type='s', lty=2, col='darkorange')
points(tti, ide, type='s', lwd=1.5, col='darkblue')
points(tti, ide+1.96*se.ide, type='s', lty=2, col='darkblue')
points(tti, ide-1.96*se.ide, type='s', lty=2, col='darkblue')
#points(tti, ide.l, type='s', lty=2, col='darkblue')
#points(tti, ide.u, type='s', lty=2, col='darkblue')
points(tti, cde, type='s', lwd=1.5, col='darkgreen')
points(tti, cde+1.96*se.cde, type='s', lty=2, col='darkgreen')
points(tti, cde-1.96*se.cde, type='s', lty=2, col='darkgreen')
#points(tti, cde.l, type='s', lty=2, col='darkgreen')
#points(tti, cde.u, type='s', lty=2, col='darkgreen')
legend('topleft', cex=0.9, lty=c(1,1,1), lwd=c(2,2,2), 
       col=c('darkorange','darkblue','darkgreen'),
       legend=c('Separbale effect','Interventional effect','Cond. inter. effect'))

plot(tti, sie, type='s', lwd=1.5, ylim=c(-.5,.3), col='darkorange',
     ylab='Diff in cumulative incidence', xlab='Days since transplantation', main='Indirect effect')
abline(h=0, lty=5)
points(tti, sie+1.96*se.sie, type='s', lty=2, col='darkorange')
points(tti, sie-1.96*se.sie, type='s', lty=2, col='darkorange')
#points(tti, sie.l, type='s', lty=2, col='darkorange')
#points(tti, sie.u, type='s', lty=2, col='darkorange')
points(tti, iie, type='s', lwd=1.5, col='darkblue')
points(tti, iie+1.96*se.iie, type='s', lty=2, col='darkblue')
points(tti, iie-1.96*se.iie, type='s', lty=2, col='darkblue')
#points(tti, iie.l, type='s', lty=2, col='darkblue')
#points(tti, iie.u, type='s', lty=2, col='darkblue')
points(tti, cie, type='s', lwd=1.5, col='darkgreen')
points(tti, cie+1.96*se.cie, type='s', lty=2, col='darkgreen')
points(tti, cie-1.96*se.cie, type='s', lty=2, col='darkgreen')
#points(tti, cie.l, type='s', lty=2, col='darkgreen')
#points(tti, cie.u, type='s', lty=2, col='darkgreen')
legend('topleft', cex=0.9, lty=c(1,1,1), lwd=c(2,2,2), 
       col=c('darkorange','darkblue','darkgreen'),
       legend=c('Separbale effect','Interventional effect','Cond. inter. effect'))
