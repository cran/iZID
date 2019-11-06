
#'Maximum likelihood estimate for hurdle or zero-inflated beta negative binomial distribution.
#'
#'@inheritParams bnb.mle
#'@rdname bb.zihmle
#' @export
bnb.zihmle<-function(x,r,alpha1,alpha2,type=c('zi','h'),lowerbound=1e-2,upperbound=1e4){
  N=length(x)
  t=x[x>0]
  m=length(t)
  neg.log.lik<-function(y){
    r1=y[1]
    a1=y[2]
    b1=y[3]
    logA=lgamma(a1+r1+b1)+lgamma(a1)
    logB=lgamma(a1+r1)+lgamma(a1+b1)
    ans=m*log(1-exp(logB-logA))+m*lgamma(a1)-m*lgamma(a1+r1)-sum(lgamma(r1+t))-sum(lgamma(b1+t))-m*lgamma(a1+b1)+sum(lgamma(t+1))
	    +m*lgamma(r1)+sum(lgamma(a1+r1+b1+t))+m*lgamma(b1)  ##ans=log(multiply_{1}^{m}{P(z_i=k,k>0)/(1-P(Z=0))}). uses the trick that gamma(n)=(n-1)! for integers
    return(ans)  ##pay:m*lgamma(a+r+b) will be eliminated by the term in logA
  }
  gp<-function(y){
    r1=y[1]
    a1=y[2]
    b1=y[3]
    logA=lgamma(a1+r1+b1)+lgamma(a1)
    logB=lgamma(a1+r1)+lgamma(a1+b1)
    dr=-m*exp(logB-logA)*(digamma(a1+r1)-digamma(a1+r1+b1))/(1-exp(logB-logA))-sum(digamma(r1+t))-m*digamma(a1+r1)+m*digamma(r1)
	   +sum(digamma(a1+r1+b1+t))
    da=-m*exp(logB-logA)*(digamma(a1+r1)+digamma(a1+b1)-digamma(a1+r1+b1)-digamma(a1))/(1-exp(logB-logA))+m*digamma(a1)-m*digamma(a1+r1)
	   -m*digamma(a1+b1)+sum(digamma(a1+r1+b1+t))
    db=-m*exp(logB-logA)*(digamma(a1+b1)-digamma(a1+r1+b1))/(1-exp(logB-logA))-sum(digamma(t+b1))-m*digamma(a1+b1)+sum(digamma(a1+r1+t+b1))
	   +m*digamma(b1)
    return(c(dr, da, db))
  }

  estimate1=stats::optim(par=c(r,alpha1,alpha2),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(lowerbound,lowerbound,lowerbound),
            upper = c(upperbound,upperbound,upperbound))
  ans=estimate1$par
  fvalue=estimate1$value
  p0=beta(ans[1]+ans[2],ans[3])/beta(ans[2],ans[3])  ##the probability of z=0
  if((!is.na(p0))&(((m/N)<=(1-p0))|(type=='h'))){
    phi=1-m/N/(1-p0)
	lik=-fvalue+(N-m)*log(1-m/N)+m*log(m/N)  #log likelihood
	if(type=='h')
	  phi=1-m/N
	mle=matrix(c(ans[1],ans[2],ans[3],phi,lik),nrow=1)
	colnames(mle)=c('r','alpha1','alpha2','phi','loglik')
	return(mle)
  }

##if (m/N)<=(1-p0) is not satisfied or type is not hurdle
  warning('cannot obtain mle with the current model type, the output estimate is derived from general beta negative binomial distribution.')
  ans=bnb.mle(x,r,alpha1,alpha2,lowerbound,upperbound)
  mle=matrix(c(ans[1],ans[2],ans[3],0,ans[4]),nrow=1)
  colnames(mle)=c('r','alpha1','alpha2','phi','loglik')
  return(mle)
}
