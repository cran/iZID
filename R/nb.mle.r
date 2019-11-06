#'Maximum likelihood estimate for negative binomial distribution.
#'
#' @param p An initial value of the probability of success, should be a positive value within (0,1).
#'
#' @export
#'
#' @rdname bb.mle
nb.mle<-function(x,r,p,lowerbound=1e-2,upperbound=1e4){
 N=length(x)
 sx=sum(x)
 neg.log.lik<-function(y) { ##negative loglikelihood function
    r1=y[1]
    p1=y[2]
    ans=-N*r1*log(p1)+N*lgamma(r1)-sx*log(1-p1)-sum(lgamma(x+r1))+sum(lgamma(x+1))
    return(ans)
  }
 gp<-function(y) {
    r2=y[1]
    p2=y[2]
    dr=-N*log(p2)+N*digamma(r2)-sum(digamma(x+r2))
    dp=sx/(1-p2)-N*r2/p2
    return(c(dr,dp))
  }
 estimate=stats::optim(par=c(r,p),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(lowerbound,lowerbound), upper = c(upperbound,1-lowerbound))
 mle=matrix(c(estimate$par[1],estimate$par[2],-estimate$value),nrow=1)
 colnames(mle)=c('r','p','loglik')
 return(mle)
}
