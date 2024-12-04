
source('phfit_npmle.R')
dat = read.csv('leukemia.csv')
dat = transform(dat, MRD=as.numeric(MRDPRET>0), SEX=SEX1-1,
                CR=DISEASESTATUS-1, DIAGNOSIS=TALLORBALL-1,
                TRANSPLANT=YZLX-1)
dat = dat[dat$MRD==1,]
A = 1-dat$TRANSPLANT
Tc = dat$AGVHDT
Dc = as.numeric(dat$AGVHD>0)
Tc[Dc==0] = dat$CGVHDT[Dc==0]
Dc[Dc==0] = as.numeric(dat$CGVHD[Dc==0]>0)
Tm = dat$RELAPSET
Dm = dat$RELAPSE
Td = dat$OST
Dd = dat$OS
#Tc[764] = Td[764]; Dc[764] = 0  #record error
Tc[Tc==Tm & Dc==1] = Tc[Tc==Tm & Dc==1] - 0.1
Td[Tc==Td & Dc==1] = Td[Tc==Td & Dc==1] + 0.1
Td[Tm==Td & Dm==1] = Td[Tm==Td & Dm==1] + 0.1
Tg = Tc; Dg = Dc
Tr = Tm; Dr = Dm 
n = nrow(dat)
dat$Tr=Tr;dat$Td=Td;dat$Tg=Tg
dat$Dr=Dr;dat$Dd=Dd;dat$Dg=Dg
dat$A = A

X = dat[,c('AGE','SEX','CR','DIAGNOSIS')]
dat$X = X
X = as.matrix(X)
p = ncol(X)
tt = sort(unique(c(Td,Tr,Tg)))
l = length(tt)

matchy = function(x,y,newx){
  newy = sapply(newx, function(t) y[which(x==t)])
  newy = as.numeric(newy)
  newy[is.na(newy)] = 0
  return(newy)
}

cif11.bs = cif00.bs = cif01.bs = cif10.bs = NULL
#for (bs in 1:200){
#ss = sample(n,n,replace=TRUE)
#A = dat$A[ss]
#X = as.matrix(dat$X[ss,])
#Tg=dat$Tg[ss];Dg=dat$Dg[ss];Tr=dat$Tr[ss]
#Dr=dat$Dr[ss];Td=dat$Td[ss];Dd=dat$Dd[ss]

# hazard of d
fit_d = phfit_d(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1)
Xb = as.numeric(X%*%fit_d$beta)
delta_g = fit_d$delta_g
delta_r = fit_d$delta_r
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_od1 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
lam_ogd1 = lam_od1 * exp(delta_g)
lam_ord1 = lam_od1 * exp(delta_r)
lam_ogrd1 = lam_orgd1 = lam_od1 * exp(delta_g+delta_r)
fit_d = phfit_d(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0)
Xb = as.numeric(X%*%fit_d$beta)
delta_g = fit_d$delta_g
delta_r = fit_d$delta_r
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_od0 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
lam_ogd0 = lam_od0 * exp(delta_g)
lam_ord0 = lam_od0 * exp(delta_r)
lam_ogrd0 = lam_orgd0 = lam_od0 * exp(delta_g+delta_r)
# hazard of g
fit_g = phfit_g(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1)
Xb = as.numeric(X%*%fit_g$beta)
delta_r = fit_g$delta_r
lam_g = matchy(fit_g$tt, fit_g$lam, tt)
lam_og1 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
lam_org1 = lam_og1 *exp(delta_r)
fit_g = phfit_g(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0)
Xb = as.numeric(X%*%fit_g$beta)
delta_r = fit_g$delta_r
lam_g = matchy(fit_g$tt, fit_g$lam, tt)
lam_og0 = sapply(1:l, function(t) lam_g[t]*exp(Xb))
lam_org0 = lam_og0 *exp(delta_r)
# hazard of r
fit_r = phfit_r(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1)
Xb = as.numeric(X%*%fit_r$beta)
delta_g = fit_r$delta_g
lam_r = matchy(fit_r$tt, fit_r$lam, tt)
lam_or1 = sapply(1:l, function(t) lam_r[t]*exp(Xb))
lam_ogr1 = lam_or1 * exp(delta_g)
fit_r = phfit_r(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0)
Xb = as.numeric(X%*%fit_r$beta)
delta_g = fit_r$delta_g
lam_r = matchy(fit_r$tt, fit_r$lam, tt)
lam_or0 = sapply(1:l, function(t) lam_r[t]*exp(Xb))
lam_ogr0 = lam_or0 * exp(delta_g)
# hazard of c
fit_c = phfit_c(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=1)
Xb = as.numeric(X%*%fit_c$beta)
lam_c = matchy(fit_c$tt, fit_c$lam, tt)
lam_c1 = sapply(1:l, function(t) lam_c[t]*exp(Xb))
fit_c = phfit_c(Tg,Dg,Tr,Dr,Td,Dd,A,X,a=0)
Xb = as.numeric(X%*%fit_c$beta)
lam_c = matchy(fit_c$tt, fit_c$lam, tt)
lam_c0 = sapply(1:l, function(t) lam_c[t]*exp(Xb))

fit = glm(A~X, family='binomial')
ps = matrix(1,nrow=n,ncol=2)
pscore = predict(fit, type='response')
ps[,1] = (1 - pscore) #* mean((1-A)/(1-pscore))
ps[,2] = pscore #* mean(A/pscore)

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
lam_c_A = A*lam_c1 + (1-A)*lam_c0
Lam_c_A = t(apply(lam_c_A, 1, cumsum))
SC = 1 - t(apply(exp(-Lam_c_A)*lam_c_A, 1, cumsum))

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

EIF = function(a){
  # counterfactual incidence
  lam_og_a = a['og']*lam_og1 + (1-a['og'])*lam_og0
  lam_od_a = a['od']*lam_od1 + (1-a['od'])*lam_od0
  lam_or_a = a['or']*lam_or1 + (1-a['or'])*lam_or0
  lam_ogr_a = a['gr']*lam_ogr1 + (1-a['gr'])*lam_ogr0
  lam_ogd_a = a['gd']*lam_ogd1 + (1-a['gd'])*lam_ogd0
  lam_org_a = a['rg']*lam_org1 + (1-a['rg'])*lam_org0
  lam_ord_a = a['rd']*lam_ord1 + (1-a['rd'])*lam_ord0
  lam_ogrd_a = a['rd']*lam_ogrd1 + (1-a['rd'])*lam_ogrd0
  lam_orgd_a = a['gd']*lam_orgd1 + (1-a['gd'])*lam_orgd0
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
  Fd = F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a
  
  return(list(Freg=Fd))
}

a1 = rep(1,7)
a0 = rep(0,7)
names(a1) = names(a0) = c('og','or','od','gr','rg','gd','rd')
index = which(tt<=2500)
tti = tt[index]
a01 = a10 = a0
a01['od'] = a01['rd'] = a01['gd'] = 1
a10['od'] = a10['rd'] = a10['gd'] = a10['og'] = a10['rg'] = 1

fit11 = EIF(a1)
fit00 = EIF(a0)
fit01 = EIF(a01)
fit10 = EIF(a10)

cif11.bs = rbind(cif11.bs, colMeans(fit11$Freg))
cif00.bs = rbind(cif00.bs, colMeans(fit00$Freg))
cif01.bs = rbind(cif01.bs, colMeans(fit01$Freg))
cif10.bs = rbind(cif10.bs, colMeans(fit10$Freg))
#}

cif1 = colMeans(fit11$Freg)[index]
cif0 = colMeans(fit00$Freg)[index]
cif2 = colMeans(fit01$Freg)[index]
se1 = apply(cif11.bs,2,sd,na.rm=TRUE)[index]
se0 = apply(cif00.bs,2,sd,na.rm=TRUE)[index]
se2 = apply(cif01.bs,2,sd,na.rm=TRUE)[index]

plot(tti, cif1, type='s', lwd=2, ylim=c(0,1), col='brown',
     ylab='Cumulative incidence', xlab='Time', main='Incidence of the primary event')
points(tti, cif0, type='s', lwd=2, col='darkcyan')
points(tti, cif1+1.96*se1, type='s', lty=2, col='brown')
points(tti, cif1-1.96*se1, type='s', lty=2, col='brown')
points(tti, cif0+1.96*se0, type='s', lty=2, col='darkcyan')
points(tti, cif0-1.96*se0, type='s', lty=2, col='darkcyan')
points(tti, cif2, type='s', lwd=2, col='darkblue')
points(tti, cif2+1.96*se2, type='s', lty=2, col='darkblue')
points(tti, cif2-1.96*se2, type='s', lty=2, col='darkblue')
legend('topleft', cex=0.8, lty=c(1,1,1), lwd=c(2,2,2), 
       col=c('brown','darkcyan','darkblue'),
       legend=c('E{Y(1,1)} MSDT','E{Y(0,0)} Haplo-SCT','E{Y(0,1)} Cross-world'))

nde = cif2 - cif0
nie = cif1 - cif2
se.nde = apply(cif01.bs-cif00.bs,2,sd,na.rm=TRUE)[index]
se.nie = apply(cif11.bs-cif01.bs,2,sd,na.rm=TRUE)[index]
plot(tti, nde, type='s', lwd=2, ylim=c(-.4,.4), col='brown',
     ylab='Diff in cumulative incidence', xlab='Time', main='Treatment effect')
abline(h=0, lty=5)
points(tti, nde+1.96*se.nde, type='s', lty=2, col='brown')
points(tti, nde-1.96*se.nde, type='s', lty=2, col='brown')
points(tti, nie, type='s', lwd=2, col='darkcyan')
points(tti, nie+1.96*se.nie, type='s', lty=2, col='darkcyan')
points(tti, nie-1.96*se.nie, type='s', lty=2, col='darkcyan')
legend('topleft', cex=0.8, lty=c(1,1), lwd=c(2,2), 
       col=c('brown','darkcyan'),
       legend=c('NDE','NIE'))
