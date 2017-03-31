
#
# This file contains R functions implementing the algorithms described  
# in the paper: Lednicki B., Wieczorkowski R.:
# OPTIMAL STRATIFICATION AND SAMPLE ALLOCATION BETWEEN SUBPOPULATIONS
# AND STRATA, Statistics in Transition (2003), Vol. 6, No. 2, 287-305.
#
# actualization date: 31.03.2017
#


cumf<-function(x,L,cumtype="Singh")
# Initial univariate stratifications using cum-rules
# after log transformation of the variable
# Arguments:
#   x - vector with stratification variable values
#   L - number of strata
#   cumtype - "D&H" (cum-square root rule) or
#             "Singh" (cum-cube root rule)
# Value:
#   vector of the initial stratification boundaries
#

{
  eps<-1e-3
  hst<-hist(log(eps+x),plot=FALSE,breaks="Scott")
  mids<-hst$mids
  fy<-hst$dens
  if (cumtype=="D&H") cf<-cumsum(sqrt(fy))
  else if (cumtype=="Singh") cf<-cumsum((fy)^(1/3))
  else { cat("cumtype error ! ","\n"); break }

  min<-min(cf)
  max<-max(cf)
  delta<-(max-min)/L
  
  gr<-double(L-1)
  xi<-double(L-1)
  for (i in 1:(L-1))
  {
     gr[i]<-min+i*delta
	 xi<-min(which(gr[i]<=cf))
	 gr[i]<-mids[xi]
  }

  return(exp(gr)-eps)
}




cumf.power<-function(x,L,r)
# Initial univariate stratifications using cum-power rule 
# (generalized cum-square root rule),
# square root (f)^(1/2) generalize to power (f)^r, where 0<r<1 
# Arguments:
#   x - vector with stratification variable values
#   L - number of strata
#   r - real number in (0,1) interval
#
# Value:
#   vector of the initial stratification boundaries
#

{
  hst<-hist(x,plot=FALSE,nclass=nclass.FD)
  mids<-hst$mids
  fy<-hst$dens
  
  cf<-cumsum((fy)^(r))
  
  min<-min(cf)
  max<-max(cf)
  delta<-(max-min)/L
  
  gr<-double(L-1)
  xi<-double(L-1)
  for (i in 1:(L-1))
  {
    gr[i]<-min+i*delta
    xi<-min(which(gr[i]<=cf))
    gr[i]<-mids[xi]
  }
  
  return(gr)
  
}





# 1d section
# functions for univariate stratifications


nh<-function(gr,L,c,x,type="Neyman",p=0.7,cv=FALSE)
# Function giving sample allocation to strata 
# Arguments:
#    gr - vector with stratum boundaries
#    L  - the number of strata
#    c  - target coefficient of variation
#    x  - vector with stratification variable
#    type - type of allocation rule: "Neyman" or "power"
#    p  - parameter for power allocation rule
#    cv - if TRUE then we compute final coefficient of variation 
# Values:
#    vector with sample sizes allocated to strata
#

{

if (!identical(gr,unique(gr))) gr<-jitter(gr)
  
h<-cut(x,c(min(x),sort(gr),max(x)),include.lowest=TRUE,labels=FALSE)

Ybar<-mean(x)
N<-length(x)
Nh<-table(h)

Wh<-Nh/N
S2h<-tapply(x,h,var)
S2h[is.na(S2h)]<-0
Yh<-tapply(x,h,mean)

ah<-double(L)
nh<-double(L)

if (type=="Neyman")
{
   ah<-Wh*sqrt(S2h)/sum(Wh[1:(L-1)]*sqrt(S2h[1:(L-1)]))
}
else
{
  if (type=="power")
       ah<-(Wh*Yh)^p/sum((Wh[1:(L-1)]*Yh[1:(L-1)])^p)
  else
  cat("Bad type in nc","\n")

}
n<-Nh[L]+sum( Wh[1:(L-1)]^2*S2h[1:(L-1)]/ah[1:(L-1)],na.rm=TRUE )/((c^2)*(Ybar^2)+sum(Wh[1:(L-1)]*S2h[1:(L-1)],na.rm=TRUE)/N)

nh<-(n-Nh[L])*ah
nh[L]<-Nh[L]
nh<-round(nh)

nh<-pmax(2,nh)
nh<-pmin(nh,Nh)

if (any(is.na(nh))) 
  if (cv==TRUE) return(list(nh=99999,CV=NA))
  else return(99999)

if ( any(nh>Nh) ) cat("nh > Nh !!! ","\n")

if (cv==TRUE)
{
  cv<-sqrt(sum((Nh/nh)*(Nh-nh)*S2h))/(N*Ybar)
  cat("cv = ",cv,"\n")
  return(list(nh=nh,CV=cv))
}

else return(nh)

}




bh<-function(x,L,c,type="Neyman",p=0.7,NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
# Function giving numerical solution of the univariate
# stratification and allocation problem
# Arguments:
#    x  - vector with stratification variable
#    L  - the number of strata
#    c  - target coefficient of variation
#    type - type of allocation rule: "Neyman" or "power"
#    p  - parameter for power allocation rule
#    NM_maxit1 - max iteration for Nelder-Mead method in optimization
#                of generalized cum-square-root rule parameter
#    NM_maxit2 - max iteration for Nelder-Mead method in main optimization phase
#    rel_tol - stopping rule, relative tolerance for sequential optimization
#                
# Values:
#    list with vectors containing optimal strata boundaries
#    and optimal sample allocation
#
{

  r0<-0.5
  cat("Initial power parameter in cum rule: ",r0,"\n")
  
  r0<-optim(r0,
            function(r) {sum( nh( cumf.power(x,L,r[1]),L,c,x,type=type,p=p))},
            method="Nelder-Mead", control=list(maxit=NM_maxit1))$par
  cat("Optimised power parameter: ",r0,"\n")
  gr0<-cumf.power(x,L,r0) 
  nh0<-nh(gr0,L,c,x,type=type,p=p)
  cat("Sample size with inital cum-power rule: ",sum(nh0),"\n")
  
    
  sumpop<-sum(nh0)
  while (1)
  {

  gropt<-optim(gr0,function(z) {sum(nh(z,L,c,x,type=type,p=p))},
               method="Nelder-Mead",control=list(maxit=NM_maxit2))$par
 	
	gropt<-sort(gropt)

	pom<-nh(gropt,L,c,x,type=type,p=p,cv=TRUE)
  nhopt<-pom$nh
	#break
	cat("Sequential optimization: sum(nh) = ",sum(nhopt),"\n")
	#if (sumpop==sum(nhopt)) { break }
	if (abs(sumpop-sum(nhopt))/sumpop<=rel_tol) { break }
	gr0<-gropt
	sumpop<-sum(nhopt)
  }	
  cat("Optimal sample size = ",sum(nhopt),"\n")
  
  return(list(bh=c(gropt,max(x)),nh=nhopt,CV=pom$CV))

}






# 2d section
# functions for bi-variate stratifications
# these functions are generalizations of the corresponding
# univariate versions



nh2d<-function(gr,Lx,Ly,cx,cy,x,y,type="Neyman",p=0.7,cv=FALSE)
{

grx<-gr[1:(Lx-1)]
gry<-gr[Lx:(Lx+Ly-2)]

if (!identical(grx,unique(grx))) grx<-jitter(grx)
if (!identical(gry,unique(gry))) gry<-jitter(gry)

  hx<-cut(x,c(min(x),sort(grx),max(x)),include.lowest=TRUE,labels=FALSE)
  hy<-cut(y,c(min(y),sort(gry),max(y)),include.lowest=TRUE,labels=FALSE)
  h<-pmax(hx,hy)
  
  Xbar<-mean(x)
  Ybar<-mean(y)
  N<-length(x)
  Nh<-table(h)
  ##if (length(Nh) != L) cat("Too many strata L !","\n")

  Wh<-Nh/N
  S2hx<-abs(tapply(x,h,var))
  S2hx[is.na(S2hx)]<-0
  Xh<-tapply(x,h,mean)
  S2hy<-abs(tapply(y,h,var))
  S2hy[is.na(S2hy)]<-0
  Yh<-tapply(y,h,mean)
  
  LL<-length(Nh)
  
  ahx<-double(LL)
  ahy<-double(LL)
  
  if (type=="Neyman")
  {
     ahx<-Wh*sqrt(S2hx)/sum(Wh[1:(LL-1)]*sqrt(S2hx[1:(LL-1)]))
     ahy<-Wh*sqrt(S2hy)/sum(Wh[1:(LL-1)]*sqrt(S2hy[1:(LL-1)]))
  }
  else
  {
    if (type=="power")
	{
	   ahx<-(Wh*Xh)^p/sum((Wh[1:(LL-1)]*Xh[1:(LL-1)])^p)
	   ahy<-(Wh*Yh)^p/sum((Wh[1:(LL-1)]*Yh[1:(LL-1)])^p)
	}
	else
	cat("Bad type in nc","\n")	
  }
  nx<-Nh[LL]+sum((Wh[1:(LL-1)]^2*S2hx[1:(LL-1)]/ahx[1:(LL-1)])/((cx^2)*(Xbar^2)+sum(Wh[1:(LL-1)]*S2hx[1:(LL-1)])/N),na.rm=TRUE)
  ny<-Nh[LL]+sum((Wh[1:(LL-1)]^2*S2hy[1:(LL-1)]/ahy[1:(LL-1)])/((cy^2)*(Ybar^2)+sum(Wh[1:(LL-1)]*S2hy[1:(LL-1)])/N),na.rm=TRUE)
  n<-max(nx,ny)
  nh<-pmax((n-Nh[LL])*ahx,(n-Nh[LL])*ahy)
  nh[LL]<-Nh[LL]
nh<-round(nh)

nh<-pmax(2,nh)
nh<-pmin(nh,Nh)
if ( any(nh>Nh) ) cat("nh > Nh !!! ","\n")

if (cv==TRUE)
{
  cvx<-sqrt(sum((Nh/nh)*(Nh-nh)*S2hx))/(N*Xbar)
  cvy<-sqrt(sum((Nh/nh)*(Nh-nh)*S2hy))/(N*Ybar)
  cat("cvx = ",cvx,"\n")
  cat("cvy = ",cvy,"\n")
  return(list(nh=nh,CVx=cvx,CVy=cvy))
}
else return(nh)

}




bh2d<-function(x,y,Lx,Ly,cx,cy,type="Neyman",p=0.7,
               QM_niter=10,NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
{
  
  r0<-c(0.5,0.5)
  cat("Initial power parameters in cum rule: ",r0,"\n")
  
  r0<-optim(r0,
            function(r) {sum( nh2d( c(cumf.power(x,Lx,r[1]),cumf.power(y,Ly,r[2])),
                                    Lx,Ly,cx,cy,x,y,type=type,p=p)  )},
            method="Nelder-Mead", control=list(maxit=NM_maxit1))$par
  cat("Optimised power parameters: ",r0,"\n")
  gr0<-c(cumf.power(x,Lx,r0[1]),cumf.power(y,Ly,r0[2])) 
  nh0<-nh2d(gr0,Lx,Ly,cx,cy,x,y,type=type,p=p)
  cat("Sample size with initial cum-power rule: ",sum(nh0),"\n")
  

  sumpop<-sum(nh0)

  while (1)
  {
    gropt<-optim(gr0,function(z) {sum(nh2d(z,Lx,Ly,cx,cy,x,y,type=type,p=p))},
                 method="Nelder-Mead")$par

	pom<-nh2d(gropt,Lx,Ly,cx,cy,x,y,type=type,p=p,cv=TRUE)
	nhopt<-pom$nh
	##break
	cat("Sequential optimization: sum(nh) = ",sum(nhopt),"\n")

	#if (sum(nhopt)==sumpop) { break }
	if (abs(sum(nhopt)-sumpop)/sumpop<=rel_tol) { break }
	gr0<-gropt
	sumpop<-sum(nhopt)
  }
  
  cat("Optimal sample size = ",sum(nhopt),"\n")
  groptx<-sort(gropt[1:(Lx-1)])
  gropty<-sort(gropt[Lx:(Lx+Ly-2)])
 
  return(list(bhx=c(groptx,max(x)),bhy=c(gropty,max(y)),nh=nhopt,CVx=pom$CVx,CVy=pom$CVy))

}






# 3d section
# functions for tri-variate stratifications
# these functions are generalizations of the corresponding
# univariate versions


nh3d<-function(gr,Lx,Ly,Lz,cx,cy,cz,x,y,z,type="Neyman",p=0.7,cv=FALSE)
{

grx<-gr[1:(Lx-1)]
gry<-gr[Lx:(Lx+Ly-2)]
grz<-gr[(Lx+Ly-1):(Lx+Ly+Lz-3)]

if (!identical(grx,unique(grx))) grx<-jitter(grx)
if (!identical(gry,unique(gry))) gry<-jitter(gry)
if (!identical(grz,unique(grz))) grz<-jitter(grz)


  hx<-cut(x,c(min(x),sort(grx),max(x)),include.lowest=TRUE,labels=FALSE)
  hy<-cut(y,c(min(y),sort(gry),max(y)),include.lowest=TRUE,labels=FALSE)
  hz<-cut(z,c(min(z),sort(grz),max(z)),include.lowest=TRUE,labels=FALSE)
  ##h<-paste(hx,hy,sep="")
  h<-pmax(hx,hy,hz)
  
  Xbar<-mean(x)
  Ybar<-mean(y)
  Zbar<-mean(z)
  N<-length(x)
  Nh<-table(h)
  ##if (length(Nh) != L) cat("Too many strata L !","\n")

  Wh<-Nh/N
  S2hx<-abs(tapply(x,h,var))
  Xh<-tapply(x,h,mean)
  S2hy<-abs(tapply(y,h,var))
  Yh<-tapply(y,h,mean)
  S2hz<-abs(tapply(z,h,var))
  Zh<-tapply(z,h,mean)

  #cat("S2hx ",S2hx,"\n")
  #cat("S2hy ",S2hy,"\n")
  
  LL<-length(table(h))
  #cat("LL ",LL,"\n")
  
  ahx<-double(LL)
  ahy<-double(LL)
  ahz<-double(LL)
   
  if (type=="Neyman")
  {
     ahx<-Wh*sqrt(S2hx)/sum(Wh[1:(LL-1)]*sqrt(S2hx[1:(LL-1)]))
     ahy<-Wh*sqrt(S2hy)/sum(Wh[1:(LL-1)]*sqrt(S2hy[1:(LL-1)]))
     ahz<-Wh*sqrt(S2hz)/sum(Wh[1:(LL-1)]*sqrt(S2hz[1:(LL-1)]))
  }
  else
  {
    if (type=="power")
	{
	   ahx<-(Wh*Xh)^p/sum((Wh[1:(LL-1)]*Xh[1:(LL-1)])^p)
	   ahy<-(Wh*Yh)^p/sum((Wh[1:(LL-1)]*Yh[1:(LL-1)])^p)
   	   ahz<-(Wh*Zh)^p/sum((Wh[1:(LL-1)]*Zh[1:(LL-1)])^p)
	}
	else
	cat("Bad type in nh3d","\n")	
  }
  nx<-Nh[LL]+sum((Wh[1:(LL-1)]^2*S2hx[1:(LL-1)]/ahx[1:(LL-1)])/((cx^2)*(Xbar^2)+sum(Wh[1:(LL-1)]*S2hx[1:(LL-1)])/N),na.rm=TRUE)
  ny<-Nh[LL]+sum((Wh[1:(LL-1)]^2*S2hy[1:(LL-1)]/ahy[1:(LL-1)])/((cy^2)*(Ybar^2)+sum(Wh[1:(LL-1)]*S2hy[1:(LL-1)])/N),na.rm=TRUE)
  nz<-Nh[LL]+sum((Wh[1:(LL-1)]^2*S2hz[1:(LL-1)]/ahz[1:(LL-1)])/((cz^2)*(Zbar^2)+sum(Wh[1:(LL-1)]*S2hz[1:(LL-1)])/N),na.rm=TRUE)
  n<-max(nx,ny,nz)
  nh<-pmax((n-Nh[LL])*ahx,(n-Nh[LL])*ahy,(n-Nh[LL])*ahz)
  nh[LL]<-Nh[LL]
nh<-round(nh)
#cat("nh ",nh,"\n")
#cat("Nh ",Nh,"\n")
#cat("n ",sum(nh),"\n")

nh<-pmax(2,nh)
nh<-pmin(nh,Nh)

if (any(is.na(nh))) 
  if (cv==TRUE) return(list(nh=99999,CV=NA))
else return(99999)

if ( any(nh>Nh) ) cat("nh > Nh !!! ","\n")

if (cv==TRUE)
{
  cvx<-sqrt(sum((Nh/nh)*(Nh-nh)*S2hx))/(N*Xbar)
  cvy<-sqrt(sum((Nh/nh)*(Nh-nh)*S2hy))/(N*Ybar)
  cvz<-sqrt(sum((Nh/nh)*(Nh-nh)*S2hz))/(N*Zbar)
  cat("cvx = ",cvx,"\n")
  cat("cvy = ",cvy,"\n")
  cat("cvz = ",cvz,"\n")
  return(list(nh=nh,CVx=cvx,CVy=cvy,CVz=cvz))
  
}
else return(nh)

}




bh3d<-function(x,y,z,Lx,Ly,Lz,cx,cy,cz,type="Neyman",p=0.7,
               QM_niter=10,NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
{
  
  r0<-c(0.5,0.5,0.5)
  cat("Initial power parameters in cum rule: ",r0,"\n")

  r0<-optim(r0,
            function(r) {sum( nh3d( c(cumf.power(x,Lx,r[1]),cumf.power(y,Ly,r[2]),
                                      cumf.power(z,Lz,r[3])),
                                    Lx,Ly,Lz,cx,cy,cz,x,y,z,type=type,p=p)  )},
            method="Nelder-Mead", control=list(maxit=NM_maxit1))$par
  cat("Optimised powe parameters:",r0,"\n")
  gr0<-c(cumf.power(x,Lx,r0[1]),cumf.power(y,Ly,r0[2]),cumf.power(z,Lz,r0[3])) 
  nh0<-nh3d(gr0,Lx,Ly,Lz,cx,cy,cz,x,y,z,type=type,p=p)
  cat("Sample size with initial cum-power rule: ",sum(nh0),"\n")
  
  
  
  sumpop<-sum(nh0)
  while (1)
  {	
    gropt<-optim(gr0,function(w) {sum(nh3d(w,Lx,Ly,Lz,cx,cy,cz,x,y,z,type=type,p=p))},
                 method="Nelder-Mead",
                 control=list(maxit=NM_maxit2))$par
    
    pom<-nh3d(gropt,Lx,Ly,Lz,cx,cy,cz,x,y,z,type=type,p=p,cv=TRUE)
    nhopt<-pom$nh
    
    # break
    cat("Sequential optimization: sum(nh) = ",sum(nhopt),"\n")
    if (abs(sum(nhopt)-sumpop)/sumpop<=rel_tol) { break }
    gr0<-gropt
    sumpop<-sum(nhopt)
  }	
  cat("Optimal sample size = ",sum(nhopt),"\n")
  groptx<-sort(gropt[1:(Lx-1)])
  gropty<-sort(gropt[Lx:(Lx+Ly-2)])
  groptz<-sort(gropt[(Lx+Ly-1):(Lx+Ly+Lz-3)])
  
  return(list(bhx=c(groptx,max(x)),bhy=c(gropty,max(y)),bhz=c(groptz,max(z)),nh=nhopt,
                      CVx=pom$CVx,CVy=pom$CVy,CVz=pom$CVz))
  
}





# # example of solving univariate stratification and
# # allocation problem
# ex1<-bh(x,5,0.01,type="Neyman",NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
# ex1<-as.data.frame(ex1)
# print(ex1)
# 
# 
# # example of solving bi-variate stratification and
# # allocation problem
# ex2<-bh2d(x,y,7,7,0.01,0.01,type="Neyman",NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
# ex2<-as.data.frame(ex2)
# print(ex2)
# 
# 
# # example of solving tri-variate stratification and
# # allocation problem
# ex3<-bh3d(x,y,z,7,7,7,0.01,0.01,0.01,type="Neyman",NM_maxit1=10,NM_maxit2=100,rel_tol=0.01)
# ex3<-as.data.frame(ex3)
# print(ex3)


