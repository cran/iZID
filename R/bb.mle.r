
#' Maximum likelihood estimate for beta binomial distributions
#'
#' calculate maximum likelihood estimate and the corresponding log likelihood value for beta binomial,
#' beta negative binomial, negative binomial and Poisson distributions.
#'
#' \code{bb.mle}, \code{bnb.mle}, \code{nb.mle} and \code{poisson.mle} calculate the maximum likelihood estimate of beta binomial,
#' beta negative binomial, negative binomial and Poisson distributions, respectively.
#'
#' Please NOTE that the arguments in the four functions are NOT CHECKED AT ALL! The user must be aware of their inputs to avoid
#' getting suspicious results.
#'
#' Suppose that \eqn{X} is a random count variable that only takes non-negative values. If \eqn{p} has a prior distribution
#' \eqn{beta(alpha1,alpha2)} and \eqn{X} follows a binomial distribution \eqn{b(n,p)}, then \eqn{X} follows the beta binomial
#' distribution with
#'
#' \eqn{P(X=k)=C(n,k)Beta(k+alpha1,n-k+alpha2)/Beta(alpha1,alpha2)},
#'
#' where \eqn{C(,)} is the combination function, \eqn{Beta(,)} is the beta function and \eqn{beta(,)} stands for the beta distribution.
#'
#' If \eqn{X} stands for the number of failures observed before the \eqn{r}th success, the probability of \eqn{X} taking the
#' value \eqn{k} under the negative binomial distribution equals
#'
#' \eqn{P(X=k)=C(k+r-1,k)p^r(1-p)^k},
#'
#' As in beta binomial distribution, assume the prior distribution of \eqn{p} is \eqn{beta(alpha1,alpha2)}. \eqn{X} follows
#' a beta negative binomial distribution if \eqn{X} follows a negative binomial distribution with parameters \eqn{r} and \eqn{p}.
#' The probability density function of a beta negative binomial distribution is defined as:
#'
#'\eqn{P(X=k)=\Gamma(r+k)Beta(r+alpha1,k+alpha2)/Beta(alpha1,alpha2)/\Gamma(r)/k!},
#'
#'where \eqn{\Gamma} represents the Gamma function.
#'
#'With the only parameter \eqn{lambda}, the probability density function of a Poisson distribution is
#'
#' \eqn{P(X=k)=lambda^k exp(-lambda)/k!}
#'
#' The maximum likelihood estimate of all four distributions can be derived by minimizing the corresponding negative log likelihood
#' function. It is easy to deduce the sample estimate of \eqn{lambda} which is equal to the sample mean. However, it is not so
#' straightforward to solve the optimization problems of the other three distributions. Thus, we adopt the optimization
#' algorithm "L-BFGS-B" by calling R basic
#' function \code{optim}. Lower and upper bounds on the unknown parameters are required for the algorithm "L-BFGS-B", which are
#' determined by the arguments \code{lowerbound} and \code{upperbound}. But note that for the estimate of \eqn{p}, the upper bound
#' for searching is essentially \code{1-lowerbound}.
#'
#' @param x A vector of count data. Should be non-negative integers.
#' @param n An initial value of the number of trials. Must be a positive number, but not required to be an integer.
#' @param alpha1 An initial value for the first shape parameter of beta distribution. Should be a positive number.
#' @param alpha2 An initial value for the second shape parameter of beta distribution. Should be a positive number.
#' @param lowerbound A lower searching bound used in the optimization of likelihood function. Should be a small positive number.
#'        The default is 1e-2.
#' @param upperbound An upper searching bound used in the optimization of likelihood function. Should be a large positive number.
#'        The default is 1e4.
#'
#' @return A row vector containing the maximum likelihood estimate of unknown parameters and the corresponding value of log likelihood.
#'
#' With \eqn{bb.mle}, the following values are returned:
#'\itemize{
#'   \item{n: the maximum likelihood estimate of n.}
#'   \item{alpha1: the maximum likelihood estimate of alpha1.}
#'   \item{alpha2: the maximum likelihood estimate of alpha2.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \eqn{bnb.mle}, the following values are returned:
#'\itemize{
#'   \item{r: the maximum likelihood estimate of r.}
#'   \item{alpha1: the maximum likelihood estimate of alpha1.}
#'   \item{alpha2: the maximum likelihood estimate of alpha2.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \eqn{nb.mle}, the following values are returned:
#'\itemize{
#'   \item{r: the maximum likelihood estimate of r.}
#'   \item{p: the maximum likelihood estimate of p.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimates plugged-in.}
#'}
#' With \eqn{poisson.mle}, the following values are returned:
#'\itemize{
#'   \item{lambda: the maximum likelihood estimate of lambda.}
#'   \item{loglik: the value of log likelihood with maximum likelihood estimate plugged-in.}
#'}
#'
#'@section Reference:
#'\itemize{
#'  \item{H. Aldirawi, J. Yang, A. A. Metwally (2019). Identifying Appropriate Probabilistic Models for Sparse Discrete Omics Data,
#'  accepted for publication in 2019 IEEE EMBS International Conference on Biomedical & Health Informatics (BHI).}
#'}
#'
#' @export
#'
#' @examples
#' x=extraDistr::rbbinom(2000,12,2,4)
#' bb.mle(x,3,1,1)
#' x=extraDistr::rbnbinom(2000,8,3,5)
#' bnb.mle(x, 3.3, 1, 1)
#' x=stats::rnbinom(2000,size=5,prob=0.3)
#' nb.mle(x, 7, 0.5)
#' x=stats::rpois(2000,7)
#' poisson.mle(x)
bb.mle<-function(x,n,alpha1,alpha2,lowerbound=1e-2,upperbound=1e4){
 N=length(x)
 neg.log.lik<-function(y){
    n1=y[1]
    a1=y[2]
    b1=y[3]
    ans=-N*lgamma(n1+1)-N*lgamma(a1+b1)+N*lgamma(a1)+N*lgamma(b1)+N*lgamma(a1+n1+b1)-sum(lgamma(x+a1))-sum(lgamma(n1-x+b1))+sum(lgamma(x+1))
	    +sum(lgamma(n1-x+1))
	return(ans)
  }
  gp<-function(y){
    n2=y[1]
    a2=y[2]
    b2=y[3]
    dn=-N*digamma(n2+1)+N*digamma(a2+n2+b2)-sum(digamma(n2-x+b2))+sum(digamma(n2-x+1))
    da=-N*digamma(a2+b2)+N*digamma(a2)+N*digamma(a2+n2+b2)-sum(digamma(x+a2))
    db=-N*digamma(a2+b2)+N*digamma(b2)+N*digamma(a2+n2+b2)-sum(digamma(n2-x+b2))
    return(c(dn,da,db))
  }
 estimate=stats::optim(par=c(n,alpha1,alpha2),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(max(x)-lowerbound,lowerbound,lowerbound),
          upper = c(upperbound,upperbound,upperbound))
 mle=matrix(c(estimate$par[1],estimate$par[2],estimate$par[3],-estimate$value),nrow=1)
 colnames(mle)=c('n','alpha1','alpha2','loglik')
 return(mle)
}

