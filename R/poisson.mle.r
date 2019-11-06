#'Maximum likelihood estimate for Poisson distribution.
#'
#' @export
#'
#' @rdname bb.mle
poisson.mle<-function(x){
 N=length(x)
 sx=sum(x)
 lambda=sx/N
 loglik=sx*log(lambda)-N*lambda-sum(lgamma(x+1))  ##the log likelihood with mle plugged in. Note gamma(k)=(k-1)!
 mle=matrix(c(lambda,loglik),nrow=1)
 colnames(mle)=c('lambda','loglik')
 return(mle)
}
