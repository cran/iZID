##, d1 is the null hypothesis, d2 is the alternative.
#' likelihood ratio test for two models
#'
#'Conduct likelihood ratio test for comparing two different models.
#'
#'If the \code{pvalue} of \code{d1} and \code{d2} are greater than the user-specified significance level, which indicates that
#'the original data \code{x} may come from the two distributions in \code{d1} and \code{d2}, a likelihood ratio test is
#'desired to choose a more 'possible' distribution based on the current data. NOTE that the \code{x} in \code{d1} and \code{d2}
#'must be IDENTICAL! Besides, NOTE that the \code{distri} in \code{d1} and \code{d2} must be DIFFERENT!
#'
#'The \code{distri} inherited from \code{d1} is the null distribution and that from \code{d2} is used as the alternative
#'distribution. Following Aldirawi et al. (2019), \code{nsim} bootstrapped or simulated samples will be generated according to
#'\code{bootstrap} of \code{d1}, based on which \code{nsim} maximum likelihood estimates of the parameters in null distribution
#'will be calculated. Remember that we have obtained \code{nsim} such maximum likelihood estimates while calling function
#'\code{dis.kstest}. Thus, the algorithm just adopts the \code{mle_new} from \code{d1} to avoid repetitive work. Using the
#'\code{nsim} maximum likelihood estimates to generate \code{nsim} new samples and calculate \code{nsim} corresponding new
#'likelihood ratio test statistics. The output p-value is the proportion of new samples that have statistics greater than the
#'test statistic of the original data \code{x}.
#'
#'As in \code{\link{dis.kstest}}, the computation is parallelized with the help of packages \code{foreach} and
#' \code{doParallel}.
#'
#'With the output p-value smaller than the user-specified significance level, the \code{distri} of \code{d2} is more appropriate
#'for modelling \code{x}. Otherwise, There is no significant difference between \code{distri} of \code{d1} and \code{distri} of \code{d2},
#'given the current data.
#'
#' @param d1 An object of class 'dis.kstest'.
#' @param d2 An object of class 'dis.kstest'.
#' @param parallel Whether to use multiple threads to parallelize computation. Default is FALSE. Please aware that it may take
#' longer time to execute the program with \code{parallel=FALSE}.
#'
#' @return The p-value of the likelihood ratio test.
#' @export
#'
#'@section Reference:
#'\itemize{
#'  \item{H. Aldirawi, J. Yang, A. A. Metwally (2019). Identifying Appropriate Probabilistic Models for Sparse Discrete Omics Data,
#'  accepted for publication in 2019 IEEE EMBS International Conference on Biomedical & Health Informatics (BHI).}
#'  \item{T. Wolodzko (2019). extraDistr: Additional Univariate and Multivariate Distributions, R package version 1.8.11,
#'   https://CRAN.R-project.org/package=extraDistr.}
#'  \item{R. Calaway, Microsoft Corporation, S. Weston, D. Tenenbaum (2017). doParallel: Foreach Parallel Adaptor for the 'parallel'
#'  Package, R package version 1.0.11, https://CRAN.R-project.org/package=doParallel.}
#'  \item{R. Calaway, Microsoft, S. Weston (2017). foreach: Provides Foreach Looping Construct for R, R package version 1.4.4,
#'   https://CRAN.R-project.org/package=foreach.}
#'}
#'
#'@seealso \code{\link{dis.kstest}}
#'
#' @examples
#' set.seed(2001)
#' temp1=sample.zi(N=300,phi=0.3,distri='poisson',lambda=5)
#'d1=dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='zip')
#' d2=dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='zinb')
#' model.lrt(d1,d2)
#'@importFrom foreach %do%

model.lrt<-function(d1,d2,parallel=FALSE){
  if(!(methods::is(d1,'dis.kstest')&methods::is(d2,'dis.kstest')))
    stop('d1, d2 must be objects from class dis.kstest.')
  if(sum((d1$x-d2$x)==0)<length(d1$x))  ##check if x of d1 is identical to x of d2
    stop('d1$x must be identical to d2$x.')

  lik1_ori=d1$mle_ori  ##the mle of original data under H0.
  length1=length(lik1_ori)
  lik2_ori=d2$mle_ori  ##the mle of original data under H1.
  length2=length(lik2_ori)
  t_ori=lik2_ori[length2]-lik1_ori[length1]##the likehood of original data with mle plugged in under H1 - that under H0. it is the test statistic of lrt.

  #initialize parameters
  mle_new=d1$mle_new  ##the mles using in step 4 of Algorithm 2.
  distri1=d1$distri
  N=d1$N  ##adopt N of d1
  lambda1=d1$lambda
  r1=d1$r
  p1=d1$p
  alpha11=d1$alpha1
  alpha12=d1$alpha2
  n1=d1$n
  phi1=d1$phi
  nsim=d1$nsim  #adopt nsim of d1
  lowerbound1=d1$lowerbound
  upperbound1=d1$upperbound
  distri2=d2$distri
  lambda2=d2$lambda
  r2=d2$r
  p2=d2$p
  alpha21=d2$alpha1
  alpha22=d2$alpha2
  n2=d2$n
  phi2=d2$phi
  lowerbound2=d2$lowerbound
  upperbound2=d2$upperbound

  #generate random deviates
  f1<-function(distri,para){
    temp1=switch(distri,poisson=stats::rpois(N,lambda=para[1]),
	             nb=stats::rnbinom(N,size=ceiling(para[1]),prob=para[2]),
	             bb=extraDistr::rbbinom(N,size=ceiling(para[1]),alpha=para[2],beta=para[3]),
	             bnb=extraDistr::rbnbinom(N,size=ceiling(para[1]),alpha=para[2],beta=para[3]),
	             zip=sample.zi(N,para[2],lambda=para[1]),
				 zinb=sample.zi(N,para[3],distri='nb',r=para[1],p=para[2]),
	             zibb=sample.zi(N,para[4],distri='bb',alpha1=para[2],alpha2=para[3],n=para[1]),
	             zibnb=sample.zi(N,para[4],distri='bnb',r=para[1],alpha1=para[2],alpha2=para[3]),
	             ph=sample.h(N,para[2],lambda=para[1]),
				 nbh=sample.h(N,para[3],distri='nb',r=para[1],p=para[2]),
	             bbh=sample.h(N,para[4],distri='bb',alpha1=para[2],alpha2=para[3],n=para[1]),
	             bnbh=sample.h(N,para[4],distri='bnb',r=para[1],alpha1=para[2],alpha2=para[3]))

	return(temp1)
  }

  #calculate mle
  f2<-function(x,distri,lambda,r,p,alpha1,alpha2,n,phi,lowerbound,upperbound){
  	temp2=switch(distri,poisson=poisson.mle(x),nb=nb.mle(x,r,p,lowerbound,upperbound),bb=bb.mle(x,n,alpha1,alpha2,lowerbound,upperbound),
	bnb=bnb.mle(x,r,alpha1,alpha2,lowerbound,upperbound),zip=poisson.zihmle(x,type='zi',lowerbound,upperbound),
	zinb=nb.zihmle(x,r,p,type='zi',lowerbound,upperbound),zibb=bb.zihmle(x,n,alpha1,alpha2,type='zi',lowerbound,upperbound),
	zibnb=bnb.zihmle(x,r,alpha1,alpha2,type='zi',lowerbound,upperbound),ph=poisson.zihmle(x,type='h',lowerbound,upperbound),
	nbh=nb.zihmle(x,r,p,type='h',lowerbound,upperbound),bbh=bb.zihmle(x,n,alpha1,alpha2,type='h',lowerbound,upperbound),
	bnbh=bnb.zihmle(x,r,alpha1,alpha2,type='h',lowerbound,upperbound))

	return(temp2)
  }

  if(parallel){
    cl_cores=parallel::detectCores()
    cl=parallel::makeCluster(cl_cores-2)
  }else{
    cl=parallel::makeCluster(1)
  }
  doParallel::registerDoParallel(cl)
  j=0

  t_new=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
        {
		  x=f1(distri1,mle_new[,j])  #simulate random deviates under d1
		  new1=f2(x,distri1,lambda1,r1,p1,alpha11,alpha12,n1,phi1,lowerbound1,upperbound1)  #calculate mle of d1
          new2=f2(x,distri2,lambda2,r2,p2,alpha21,alpha22,n2,phi2,lowerbound2,upperbound2)
		  dif=new2[length2]-new1[length1]   ##new test statistics in step 5 of Algorithm 2.
	    }

  parallel::stopCluster(cl)
  pvalue=sum(t_new>t_ori)/nsim
  return(pvalue)
}
