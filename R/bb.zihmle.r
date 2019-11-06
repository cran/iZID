
#' Maximum likelihood estimate for zero-inflated or hurdle beta binomial distributions.
#'
#' Calculate maximum likelihood estimate and the corresponding log likelihood value for zero-inflated or hurdle
#'  beta binomial, beta negative binomial, negative binomial and Poisson distributions.
#'
#'  By setting \code{type='zi'}, \code{bb.zihmle}, \code{bnb.zihmle}, \code{nb.zihmle} and \code{poisson.zihmle} calculate the
#'  maximum likelihood estimate of zero-inflated beta binomial, beta negative binomial, negative binomial and Poisson
#'  distributions, respectively.
#'
#'  By setting \code{type='h'}, \code{bb.zihmle}, \code{bnb.zihmle}, \code{nb.zihmle} and \code{poisson.zihmle} calculate the
#'  maximum likelihood estimate of hurdle beta binomial, beta negative binomial, negative binomial and Poisson
#'  distributions, respectively.
#'
#' Please NOTE that the arguments in the four functions are NOT CHECKED AT ALL! The user must be aware of their inputs to avoid
#' getting suspicious results.
#'
#' For zero-inflated models, zeros occurred by either sampling process or specific structure of data with the structural parameter
#' \eqn{0<\phi<1}. The density function for a zero-inflated model is
#'
#' \eqn{P_{zi}(X=k)=\phi 1_{k=0}+(1-\phi)P(X=k)},
#'
#' where \eqn{P(X=k)} is the probability under standard distributions.
#'
#' Aldirawi et al. (2019) proposed an estimating procedure for zero-inflated models by optimizing over a reparametrization of
#' the likelihood function where \eqn{\phi} and the rest unknown parameters are separable. When \eqn{X} comes from a
#' zero-inflated distribution, the maximum likelihood estimate of parameters except for \eqn{\phi} are obtained by minimizing
#' the truncated version of negative log likelihood function. However, in the zero-deflated case, \eqn{\phi=0} and the sample
#' estimate of other parameters are identical to those for its corresponding standard distributions. Meanwhile, an warning
#' message is shown on the screen such that 'cannot obtain mle with the current model type, the output estimate is derived from
#' general ... distribution'.
#'
#' For hurdle models, all zeros occurred purely by the structure of data with the structural parameter \eqn{0<\phi<1}.
#' The density function for a hurdle model is
#'
#' \eqn{P_{h}(X=k)=\phi 1_{k=0}+(1-\phi)P_{tr}(X=k)},
#'
#' where \eqn{P_{tr}(X=k)} is the truncated probability under standard distributions, where \eqn{P_{tr}(X=0)=0} and
#' \eqn{P_{tr}(X=k)=P(X=k)/(1-P(X=0))}. Since \eqn{\phi} and other unknown parameters are separable in the joint likelihood
#'  function, \eqn{\phi} can be estimated by a value with respect to the number of positive samples. The sample estimate of
#'  other parameters can be obtained by the same procedure for zero-inflated model.
#'
#'  A warning message may also occur when the algorithm of \code{optim} does not converge and the resulting estimates are not
#'  valid. In this case, the results from the corresponding general distribution are output instead.
#'
#' @inheritParams bb.mle
#' @param type The type of distribution used to calculate the sample estimate, where 'zi' and 'h' stand for zero-inflated
#'  and hurdle distributions respectively.
#'
#' @return A row vector containing the maximum likelihood estimate of unknown parameters and the corresponding value of log likelihood.
#'
#' With \code{bb.zihmle}, the following values are returned:
#'\itemize{
#'   \item{n: the maximum likelihood estimate of n.}
#'   \item{alpha1: the maximum likelihood estimate of alpha1.}
#'   \item{alpha2: the maximum likelihood estimate of alpha2.}
#'   \item{phi: the maximum likelihood estimate of \eqn{\phi}.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \code{bnb.zihmle}, the following values are returned:
#'\itemize{
#'   \item{r: the maximum likelihood estimate of r.}
#'   \item{alpha1: the maximum likelihood estimate of alpha1.}
#'   \item{alpha2: the maximum likelihood estimate of alpha2.}
#'   \item{phi: the maximum likelihood estimate of \eqn{\phi}.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \code{nb.zihmle}, the following values are returned:
#'\itemize{
#'   \item{r: the maximum likelihood estimate of r.}
#'   \item{p: the maximum likelihood estimate of p.}
#'   \item{phi: the maximum likelihood estimate of \eqn{\phi}.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \code{poisson.zihmle}, the following values are returned:
#'\itemize{
#'   \item{lambda: the maximum likelihood estimate of lambda.}
#'   \item{phi: the maximum likelihood estimate of \eqn{\phi}.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimate plugged-in.}
#'}
#'
#'@section Reference:
#'\itemize{
#'  \item{H. Aldirawi, J. Yang, A. A. Metwally (2019). Identifying Appropriate Probabilistic Models for Sparse Discrete Omics Data,
#'  accepted for publication in 2019 IEEE EMBS International Conference on Biomedical & Health Informatics (BHI).}
#'  \item{H. Aldirawi, J. Yang (2019). Model Selection and Regression Analysis for Zero-altered or Zero-inflated Data, Statistical
#'  Laboratory Technical Report, no.2019-01, University of Illinois at Chicago.}
#'}
#'
#' @export
#'
#' @examples
#'t1=sample.h(N=2000,phi=0.2,distri='Poisson',lambda=5)  ##hurdle poisson random values
#' t2=sample.h(N=2000,phi=0.2,distri='nb',r=10,p=0.6)   ##hurdle negative binomial
#' t3=sample.zi(N=2000,phi=0.2,distri='bb',alpha1=8,alpha2=9,n=10)   ##zero-inflated beta binomial
#' ##zero-inflated beta negative binomial.
#' t4=sample.zi(N=2000,phi=0.2,distri='bnb',r=10,alpha1=8,alpha2=9)
#' bb.zihmle(t3,3,1,1,type='h')
#' bnb.zihmle(t4, 3.3, 1, 1,type='h')
#' nb.zihmle(t2, 7, 0.5,type='zi')
#' poisson.zihmle(t1,type='zi')
bb.zihmle<-function(x,n,alpha1,alpha2,type=c('zi','h'),lowerbound=1e-2,upperbound=1e4){
  N=length(x)  #the sample size of x
  t=x[x>0]  #the positive samples
  m=length(t)  #the number of positive samples
  neg.log.lik<-function(y) {#negative loglikelihood function of beta binomial model. only use the nonzero samples to estimate.
    n1=y[1]
    a1=y[2]
    b1=y[3]
    logA=lgamma(a1+n1+b1)+lgamma(b1) #lgamma return the natural logarithm of the absolute value of the gamma function for all real x except zero and negative integers (when NaN is returned)..
    logB=lgamma(a1+b1)+lgamma(n1+b1)
    ans=m*log(1-exp(logB-logA))+m*logA-m*lgamma(n1+1)-sum(lgamma(t+a1))-sum(lgamma(n1-t+b1))-m*lgamma(a1+b1)+sum(lgamma(t+1))+sum(lgamma(n1-t+1))
	    +m*lgamma(a1)
    return(ans)  ##ans=log(multiply_{1}^{m}{P(z_i=k,k>0)/(1-P(Z=0))}). uses the trick that C_{n}^{k}=n!/k!/(n-k)! and gamma(n)=(n-1)! for integers
  }
  gp<-function(y) {##gradient function of neg.log.lik
    n1=y[1];
    a1=y[2];
    b1=y[3];
    logA=lgamma(a1+n1+b1)+lgamma(b1)  #digamma returns the first second derivative of the logarithm of the gamma function.
    logB=lgamma(a1+b1)+lgamma(n1+b1)
    dn=-m*exp(logB-logA)*(digamma(n1+b1)-digamma(a1+n1+b1))/(1-exp(logB-logA))-m*digamma(n1+1)-sum(digamma(n1+b1-t))+sum(digamma(n1-t+1))
	   +m*digamma(a1+n1+b1)
    da=-m*exp(logB-logA)*(digamma(a1+b1)-digamma(a1+n1+b1))/(1-exp(logB-logA))-sum(digamma(t+a1))-m*digamma(a1+b1)+m*digamma(a1+n1+b1)+m*digamma(a1)
    db=-m*exp(logB-logA)*(digamma(a1+b1)+digamma(n1+b1)-digamma(a1+n1+b1)-digamma(b1))/(1-exp(logB-logA))+m*digamma(b1)-sum(digamma(n1-t+b1))
	   -m*digamma(a1+b1)+m*digamma(a1+n1+b1)
    return(c(dn,da,db))
  }
  estimate=stats::optim(par=c(n,alpha1,alpha2),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(max(x)-lowerbound,lowerbound,lowerbound),
           upper = c(upperbound,upperbound,upperbound))
  ans=estimate$par
  fvalue=estimate$value   ##the minimal negative loglikelihood value.
  p0=beta(ans[2],ans[1]+ans[3])/beta(ans[2],ans[3])  ##the probability that z=0
  if((!is.na(p0))&(((m/N)<=(1-p0))|(type=='h'))){##if the proportion of 0 in x is greater than P(Z_i=0) in the beta binomial model.
    phi=1-m/N/(1-p0) #phi=1-(m/N)/(1-P(Z=0)) <1
    lik=-fvalue+(N-m)*log(1-m/N)+m*log(m/N)  # likelihood
	if(type=='h')
	   phi=1-m/N
	mle=matrix(c(ans[1],ans[2],ans[3],phi,lik),nrow=1)
    colnames(mle)=c('n','alpha1','alpha2','phi','loglik')
    return(mle)
  }

 ##if (m/N)<=(1-p0) is not satisfied or type is not hurdle
 warning('cannot obtain mle with the current model type, the output estimate is derived from general beta binomial distribution.')
 ans=bb.mle(x,n,alpha1,alpha2,lowerbound,upperbound)
 mle=matrix(c(ans[1],ans[2],ans[3],0,ans[4]),nrow=1)
 colnames(mle)=c('n','alpha1','alpha2','phi','loglik')
 return(mle)
}
