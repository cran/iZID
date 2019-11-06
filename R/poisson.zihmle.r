#'Maximum likelihood estimate for hurdle or zero-inflated Poisson distribution.
#'
#'@rdname bb.zihmle
#' @export
poisson.zihmle<-function(x,type=c('zi','h'),lowerbound=1e-2,upperbound=1e4){
  N=length(x)  ##sample size
  y=x[x>0]     ##positive samples
  m=length(y)   #number of positive samples
  sumy=sum(y)  ##sum of the positive samples
  f1<-function(lambda){
    temp1=m*lambda+sumy*exp(-lambda)-sumy
	return(temp1)
  }
  lambda=rootSolve::uniroot.all(f1,c(lowerbound,upperbound))
  if(length(lambda)==0)
    lambda=NA    #to avoid the lambda from uniroot.all is numeric zero.
  p0=exp(-lambda) ##the probability that x_i equals 0
  if((!is.na(p0))&(((m/N)<=(1-p0))|(type=='h'))){
    phi=1-m/N/(1-p0)  ##phi for zero-inflated models.
    loglik=(N-m)*log(1-m/N)+m*log(m/N/(1-p0))+log(lambda)*sumy-m*lambda-sum(lgamma(y+1))   ##log likelihood
	if(type=='h')
	  phi=1-m/N
  }else{
    lambda=sumy/N
    phi=0
	loglik=-N*lambda+log(lambda)*sumy-sum(lgamma(y+1))   ##log likelihood
	warning('cannot obtain mle with the current model type, the output estimate is derived from general Poisson distribution.')
  }
  mle=matrix(c(lambda,phi,loglik),nrow=1)
  colnames(mle)=c('lambda','phi','loglik')
  return(mle)
}

