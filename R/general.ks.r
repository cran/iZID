##compute KS statistic for possion,nb,bb,bnb and related zero-inflated or hurdle models with specific parameters
general.ks<-function(N,lambda,r,p,n,alpha1,alpha2,distri,type=c('general','zi','h'),phi=0){
  if((type=='general')|(phi==0)){
    x_new=switch(distri,poisson=stats::rpois(N,lambda=lambda),nb=stats::rnbinom(N,size=ceiling(r),prob=p),bb=extraDistr::rbbinom(N,size=ceiling(n),alpha=alpha1,beta=alpha2),
       bnb=extraDistr::rbnbinom(N,size=ceiling(r),alpha=alpha1,beta=alpha2))
    probs_new=switch(distri,poisson=stats::ppois(0:max(x_new),lambda=lambda),nb=stats::pnbinom(0:max(x_new),size=ceiling(r),prob=p),
       bb=extraDistr::pbbinom(0:max(x_new),size=ceiling(n),alpha=alpha1,beta=alpha2),bnb=extraDistr::pbnbinom(0:max(x_new),size=ceiling(r),alpha=alpha1,beta=alpha2))
  }else{
    if(type=='zi'){
	  x_new=switch(distri,poisson=sample.zi(N,phi,distri='poisson',lambda),nb=sample.zi(N,phi,distri='nb',r=r,p=p),
        bb=sample.zi(N,phi,distri='bb',alpha1=alpha1,alpha2=alpha2,n=n),bnb=sample.zi(N,phi,distri='bnb',r=r,alpha1=alpha1,alpha2=alpha2))
      probs_new=switch(distri,poisson=stats::ppois(0:max(x_new),lambda=lambda),nb=stats::pnbinom(0:max(x_new),size=ceiling(r),prob=p),
        bb=extraDistr::pbbinom(0:max(x_new), size=ceiling(n),alpha=alpha1,beta=alpha2),bnb=extraDistr::pbnbinom(0:max(x_new),size=ceiling(r),alpha=alpha1,beta=alpha2))	   
      probs_new=phi+(1-phi)*probs_new
	}else{
	  x_new=switch(distri,poisson=sample.h(N,phi,distri='poisson',lambda=lambda),nb=sample.h(N,phi,distri='nb',r=r,p=p),
        bb=sample.h(N,phi,distri='bb',alpha1=alpha1,alpha2=alpha2,n=n),bnb=sample.h(N,phi,distri='bnb',r=r,alpha1=alpha1,alpha2=alpha2))
      probs_new=switch(distri,poisson=stats::ppois(0:max(x_new),lambda=lambda),nb=stats::pnbinom(0:max(x_new),size=ceiling(r),prob=p),
        bb=extraDistr::pbbinom(0:max(x_new), size=ceiling(n),alpha=alpha1,beta=alpha2),bnb=extraDistr::pbnbinom(0:max(x_new),size=ceiling(r),alpha=alpha1,beta=alpha2))	   
      probs_new=phi+(1-phi)*(probs_new-probs_new[1])/(1-probs_new[1])
	}
  }
  
 step_new=stats::stepfun(0:max(x_new),c(0,probs_new))  ##step function for the new samples
 z=stats::knots(step_new)  ##Line 25-27 are taken from R code of function disc_ks_test in gsgeneral package
 dev=c(0, stats::ecdf(x_new)(z)-step_new(z))
 Dn_new=max(abs(dev)) ##KS statistic for the new samples
 return(Dn_new)
 }
    
 
#examples
# N=2000
# tol=1e-6
# r=5
# p=0.3
# n=10
# alpha1=2
# alpha2=3
# lambda=2

 # general.ks(N,tol,lambda,distri='poisson',type='general',phi=0,sample.zi,sample.h)  ##0.01185173 
 # general.ks(N,tol,lambda,distri='poisson',type='h',phi=0,sample.zi,sample.h)    #0.02080992
 # general.ks(N,tol,lambda,distri='poisson',type='zi',phi=0,sample.zi,sample.h)   #0.005787068
 # general.ks(N,tol,lambda,distri='poisson',type='h',phi=0.3,sample.zi=sample.zi,sample.h=sample.h)  #0.008422593
 # general.ks(N,tol,lambda,distri='poisson',type='zi',phi=0.3,sample.zi=sample.zi,sample.h=sample.h) #0.01076694
 # general.ks(N,tol,r=r,p=p,distri='nb',type='general',phi=0,sample.zi,sample.h)  #0.01850778 
 # general.ks(N,tol,r=r,p=p,distri='nb',type='h',phi=0.2,sample.zi=sample.zi,sample.h=sample.h) #0.01405876
 # general.ks(N,tol,r=r,p=p,distri='nb',type='zi',phi=0.2,sample.zi=sample.zi,sample.h=sample.h) #0.01145171
 # general.ks(N,tol,n=n,alpha1=alpha1,alpha2=alpha2,distri='bb',type='general',phi=0,sample.zi,sample.h) #0.02068931
 # general.ks(N,tol,n=n,alpha1=alpha1,alpha2=alpha2,distri='bb',type='zi',phi=0,sample.zi,sample.h) #0.01181069
 # general.ks(N,tol,n=n,alpha1=alpha1,alpha2=alpha2,distri='bb',type='zi',phi=0.4,sample.zi=sample.zi,sample.h=sample.h) #0.008104895
 # general.ks(N,tol,n=n,alpha1=alpha1,alpha2=alpha2,distri='bb',type='h',phi=0.4,sample.zi=sample.zi,sample.h=sample.h) #0.01878075
 # general.ks(N,tol,r=r,alpha1=alpha1,alpha2=alpha2,distri='bnb',type='general',phi=0,sample.zi,sample.h)  #0.01346866
 # general.ks(N,tol,r=r,alpha1=alpha1,alpha2=alpha2,distri='bnb',type='zi',phi=0,sample.zi,sample.h)  #0.02045238
 # general.ks(N,tol,r=r,alpha1=alpha1,alpha2=alpha2,distri='bnb',type='zi',phi=0.3,sample.zi=sample.zi,sample.h=sample.h) #0.01808392
 # general.ks(N,tol,r=r,alpha1=alpha1,alpha2=alpha2,distri='bnb',type='h',phi=0.3,sample.zi=sample.zi,sample.h=sample.h)  #0.01417337
