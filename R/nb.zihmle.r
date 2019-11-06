#'Maximum likelihood estimate for hurdle or zero-inflated negative binomial distribution.
#'
#'@inheritParams nb.mle
#'@rdname bb.zihmle
#' @export
nb.zihmle<-function(x,r,p,type=c('zi','h'),lowerbound=1e-2,upperbound=1e4){ ##to use foreach function, needs to take nb.mle as input.
  N=length(x)
  t=x[x>0]
  tsum=sum(t)
  m=length(t)
  neg.log.lik<-function(y){
    r1=y[1]
    p1=y[2]
    ans=-sum(lgamma(t+r1))+sum(lgamma(t+1))+m*lgamma(r1)-tsum*log(1-p1)-m*r1*log(p1)+m*log(1-p1^r1)
    return(ans)
  }
  gp<-function(y){
    r1=y[1]
    p1=y[2]
    g1=-sum(digamma(t+r1))+m*digamma(r1)-m*log(p1)-m*log(p1)/(1-p1^r1)*((p1)^r1)#=dL/dr
    g2=tsum/(1-p1)-m*r1/p1-m*r1/(1-p1^r1)*((p1)^(r1-1)) #=dL/dp
    return(c(g1,g2))
  }

  estimate=stats::optim(par=c(r,p),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(lowerbound,lowerbound), upper = c(upperbound,1-lowerbound))
  ans=estimate$par  ##c(\hat{r},\hat{p})
  fvalue=estimate$value  ##the estimate negative loglikelihood of truncated distribution.
  p0=(ans[2])^(ans[1])   ##the prob that z=0
  if((!is.na(p0))&(((m/N)<=(1-p0))|(type=='h'))){
    phi=1-m/N/(1-p0)
    lik=-fvalue+(N-m)*log(1-m/N)+m*log(m/N)  #log likelihood
	if(type=='h')
	  phi=1-m/N
	mle=matrix(c(ans[1],ans[2],phi,lik),nrow=1)
	colnames(mle)=c('r','p','phi','loglik')
	return(mle)
  }

  ##if m/N>(1-p^r), do (2), call for mle of general nb models.
  warning('cannot obtain mle with the current model type, the output estimate is derived from general negative binomial distribution.')
  ans=nb.mle(x,r,p,lowerbound,upperbound)
  mle=matrix(c(ans[1],ans[2],0,ans[3]),nrow=1)
  colnames(mle)=c('r','p','phi','loglik')
  return(mle)
}
