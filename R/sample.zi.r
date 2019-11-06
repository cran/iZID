#'Generate random deviates for zero-inflated models
#' @export
#'
#' @rdname sample.h

sample.zi<-function(N,phi,distri='poisson',lambda=NA,r=NA,p=NA,alpha1=NA,alpha2=NA,n=NA) {
  ##check the validity of some inputs
  distri=tolower(distri)[1]  ##so that all letters are in lower case. only use the first value if user inputs a vector but not a scalar.
  N=ceiling(N)[1]
  phi=phi[1]
  if(N<=0)
    stop('The sample size N is too small.')
  if((phi>=1)|(phi<=0))
    stop('phi is not between 0 and 1 (not including 0 and 1).')
  if(!(distri%in%c('poisson','nb','bb','bnb')))
   stop('please input a distribution name among poisson,nb,bb,bnb.')
  if(distri=="poisson"){
    lambda=lambda[1]
	if(lambda<=0)
	   stop('lambda is less or equal to 0.')
  }
  if(distri=="nb"){
    r=ceiling(r[1])
	p=p[1]
	if(r<=0)
	   stop('r is too small.')
	if((p<=0)|(p>=1))
	   stop('p is not between 0 and 1 (not including 0 and 1).')
  }
  if((distri=="bb")|(distri=="bnb")){
    alpha1=alpha1[1]
	alpha2=alpha2[1]
	if(alpha1<=0)
	   stop('alpha1 is less or equal to 0.')
	if(alpha2<=0)
	   stop('alpha2 is less or equal to 0.')
  }
  if(distri=="bb"){
    n=ceiling(n[1])
	if(n<=0)
	   stop('n is too small.')
  }
  if(distri=="bnb"){
    r=ceiling(r[1])
	if(r<=0)
	   stop('r is too small.')
  }

  ans=stats::rbinom(N,size=1,prob=phi) ##simulate Z_1,Z_2,...,Z_N from Bernoulli(phi).
  m=sum(ans==0)  ##if Z_i=1, let X_i=0, if Z_i==0, generate new samples from original distribution.
  temp1=NULL
  if(m>0){
    temp1=switch(distri,nb=stats::rnbinom(m,size=r,prob=p),bb=extraDistr::rbbinom(m,n,alpha1,alpha2),bnb=extraDistr::rbnbinom(m,r,alpha1,alpha2),stats::rpois(m,lambda))
  }
  temp2=ans
  temp2[ans==1]=0
  temp2[ans==0]=temp1
  return(temp2)
}


