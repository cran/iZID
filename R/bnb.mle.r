
#'Maximum likelihood estimate for beta negative binomial distribution.
#'
#' @param r An initial value of the number of success before which m failures are observed, where m is the element of x.
#' Must be a positive number, but not required to be an integer.
#'
#' @export
#'
#' @rdname bb.mle
bnb.mle<-function(x,r,alpha1,alpha2,lowerbound=1e-2,upperbound=1e4){
 N=length(x)
 neg.log.lik<-function(y){
    r1=y[1]
    a1=y[2]
    b1=y[3]
    ans=-N*lgamma(a1+r1)-N*lgamma(a1+b1)+N*lgamma(r1)+N*lgamma(a1)+N*lgamma(b1)-sum(lgamma(x+r1))-sum(lgamma(x+b1))+sum(lgamma(x+1))
	    +sum(lgamma(a1+r1+b1+x))
	return(ans)
  }
 gp<-function(y){
    r2=y[1]
    a2=y[2]
    b2=y[3]
    dr=-N*digamma(a2+r2)+N*digamma(r2)-sum(digamma(r2+x))+sum(digamma(a2+r2+b2+x))
    da=-N*digamma(a2+r2)-N*digamma(a2+b2)+N*digamma(a2)+sum(digamma(a2+r2+b2+x))
    db=-N*digamma(a2+b2)+N*digamma(b2)-sum(digamma(b2+x))+sum(digamma(a2+r2+b2+x))
    return(c(dr,da,db))
  }
 estimate=stats::optim(par=c(r,alpha1,alpha2),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(lowerbound,lowerbound,lowerbound),
          upper = c(upperbound,upperbound,upperbound))
 mle=matrix(c(estimate$par[1],estimate$par[2],estimate$par[3],-estimate$value),nrow=1)
 colnames(mle)=c('r','alpha1','alpha2','loglik')
 return(mle)
}
