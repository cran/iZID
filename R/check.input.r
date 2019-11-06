##-----check the validity of the inputs-----------

check.input<-function(x,nsim,distri,lowerbound,upperbound){
if(min(x)<0) 
   stop('x should be nonnegative.')
if(length(x)<=30) 
   warning('Small sample size may lead to biased or inaccurate results.')
if(length(unique(x))==1)
   stop('There must be more than one unique values in x.')
if(nsim<=30) 
   stop('nsim is too small: please input nsim larger than 30.')
if(!(distri%in%c('poisson','nb','bb','bnb','zip','zinb','zibb','zibnb','ph','nbh','bbh','bnbh'))) 
   stop('please input a distribution name among poisson,nb,bb,bnb,zip,zinb,zibb,zibnb,ph,nbh,bbh,bnbh.')
if((lowerbound<0)|(lowerbound>0.1))
  stop('lowerbound is negative or larger than 0.1')
if(upperbound<1)
  stop('upperbound is too small.')
}


