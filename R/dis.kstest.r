
#' The Monte Carlo estimate for the p-value of a discrete KS Test
#'
#'Computes the Monte Carlo estimate for the p-value of a discrete one-sample Kolmogorov-Smirnov (KS) Test for Poisson, negative
#'binomial, beta binomial, beta negative binomial distributions and their zero-inflated as well as hurdle versions.
#'
#'For arguments \code{nsim}, \code{bootstrap}, \code{distri}, if the length is larger than 1, only the first element will be used.
#'For other arguments except for \code{x}, the first valid value will be used if the input is not \code{NULL}, otherwise some
#' naive sample estimates will be fed into the algorithm. Note that only the initial values that are occurred in the null
#' distribution \code{distri} are needed. For example, with \code{distri=poisson}, user may provide a value for \code{lambda} but
#' not for \code{r} or \code{p}, though it won't disturb the algorithm.
#'
#' With an output p-value less than some user-specified significance level, \code{x} is very likely from a distribution other
#' than the \code{distri}, given the current data. If p-values of more than one distributions are greater than the pre-specified
#' significance level, user may consider a following likelihood ratio test to select a 'better' distribution.
#'
#' The methodology of computing Monte Carlo p-value is taken from Aldirawi et al. (2019). With \code{bootstrap=TRUE}, \code{nsim}
#' bootstrapped samples will be generated by resampling \code{x} without replacement. Otherwise, \code{nsim} samples are
#' simulated from the null distribution with the maximum likelihood estimate of original data \code{x}. Then compute the maximum
#' likelihood estimates of \code{nsim} bootstrapped or simulated samples, based on which \code{nsim} new samples are generated
#' under the null distribution. \code{nsim} KS statistics are calculated for the \code{nsim} new samples, then the Monte Carlo
#' p-value is resulted from comparing the \code{nsim} KS statistics and the statistic of original data \code{x}.
#'
#' During the process of computing maximum likelihood estimates, the negative log likelihood function is minimized via basic R
#'  function \code{optim} with the searching interval decided by \code{lowerbound} and \code{upperbound}, except that the optimization
#'   of \code{p} takes \code{1-lowerbound} as the upper searching bound.
#'
#' To accelerate the whole process, the algorithm uses the parallel strategy via the packages \code{foreach} and
#'   \code{doParallel}.
#'
#'
#' @param x A vector of count data. Should be non-negative integers. If elements of x are not integers, they will be
#' automatically rounded up to the smallest integers that are no less than themselves.
#' @param nsim The number of bootstrapped samples or simulated samples generated to compute p-value. If it is not an integer,
#' nsim will be automatically rounded up to the smallest integer that is no less than nsim. Should be greater than 30. Default is
#'  100.
#' @param bootstrap Whether to generate bootstrapped samples or not. See Details. 'TRUE' or any numeric non-zero value indicates
#' the generation of bootstrapped samples. The default is 'TRUE'.
#' @param distri The distribution used as the null hypothesis. Can be one of \{'poisson','nb','bb',
#' 'bnb','zip','zinb','zibb', zibnb','ph','nbh','bbh','bnbh'\}, which corresponds to Poisson, negative binomial, beta binomial
#' and beta negative binomial distributions and their zero-inflated as well as hurdle versions, respectively. Default is 'Poisson'.
#' @param r An initial value of the number of success before which m failures are observed, where m is the element of x.
#' Must be a positive number, but not required to be an integer.
#' @param p An initial value of the probability of success, should be a positive value within (0,1).
#' @param alpha1 An initial value for the first shape parameter of beta distribution. Should be a positive number.
#' @param alpha2 An initial value for the second shape parameter of beta distribution. Should be a positive number.
#' @param n An initial value of the number of trials. Must be a positive number, but not required to be an integer.
#' @param lowerbound A lower searching bound used in the optimization of likelihood function. Should be a small positive number.
#'        The default is 1e-2.
#' @param upperbound An upper searching bound used in the optimization of likelihood function. Should be a large positive number.
#'        The default is 1e4.
#' @param parallel whether to use multiple threads to parallelize computation. Default is FALSE. Please aware that it may take
#' longer time to execute the program with \code{parallel=FALSE}.
#'
#' @return An object of class 'dis.kstest' including the following elements:
#'\itemize{
#'   \item{x: \code{x} used in computation.}
#'   \item{nsim: nsim used in computation.}
#'   \item{bootstrap: bootstrap used in computation.}
#'   \item{distri: distri used in computation..}
#'   \item{lowerbound: lowerbound used in computation.}
#'   \item{upperbound: upperboound used in computation.}
#'   \item{mle_new: A matrix of the maximum likelihood estimates of unknown parameters under the null distribution, using
#'   \eqn{nsim} bootstrapped or simulated samples.}
#'   \item{mle_ori: A row vector of the maximum likelihood estimates of unknown parameters under the null distribution, using the
#'   original data \code{x}.}
#'   \item{pvalue: Monte Carlo p-value of the one-sample KS test.}
#'   \item{N: length of \code{x}.}
#'   \item{r: initial value of r used in computation.}
#'   \item{p: initial value of p used in computation.}
#'   \item{alpha1: initial value of alpha1 used in computation.}
#'   \item{alpha2: initial value of alpha2 used in computation.}
#'   \item{n: initial value of n used in computation.}
#'}
#'
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
#'@seealso \code{\link{model.lrt}}
#'
#' @examples
#'
#' set.seed(2001)
#' temp1=sample.zi(N=300,phi=0.3,distri='poisson',lambda=5)
#' dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='Poisson')$pvalue
#' dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='nb')$pvalue
#' dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='zip')$pvalue
#' dis.kstest(temp1,nsim=100,bootstrap=TRUE,distri='zinb')$pvalue
#'
#'@importFrom foreach %do%
#'
dis.kstest<-function(x,nsim=100,bootstrap=TRUE,distri='Poisson',r=NULL,p=NULL,alpha1=NULL,alpha2=NULL,n=NULL,lowerbound=1e-2,upperbound=1e4,
 parallel=FALSE){
 x=ceiling(as.vector(x))  ##to make x integers
 nsim=ceiling(nsim[1])
 bootstrap=bootstrap[1]
 distri=tolower(distri)[1]  ##to make all letters in distri be of the lower case.
 check.input(x,nsim,distri,lowerbound,upperbound)  ##check the validity of x,nsim,distri
 init=init.para(x,r,p,alpha1,alpha2,n)  #obtain the initial values of parameters
 r=init$r
 p=init$p
 alpha1=init$alpha1
 alpha2=init$alpha2
 n=init$n
 N=length(x)

 ##initialization for parallel computation
 if(parallel){
    cl_cores=parallel::detectCores()
    cl=parallel::makeCluster(cl_cores-2)
 }else{
    cl=parallel::makeCluster(1)
 }
 doParallel::registerDoParallel(cl)
 j=0

##--------for possion distribution------
 if(distri=='poisson'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=poisson.mle(x)  ##mle using the original data, \hat{\theta} in Algorithm 1.
    probs_ori=stats::ppois(0:max(x),lambda=mle_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
	z=stats::knots(step_ori)  ##Line 38-40 are taken from R code of function disc_ks_test in gsgeneral package
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           poisson.mle(sample(x, size=N, replace=T))
	   mle_new=t(mle_new)
    }else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,lambda=mle_new[1,j],distri='poisson',type='general')

	temp=list(r=NULL,p=NULL,alpha1=NULL,alpha2=NULL,n=NULL)
 }

##--------for nb distribution------
 if(distri=='nb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=nb.mle(x,r,p,lowerbound,upperbound)  ##mle using the original data
    probs_ori=stats::pnbinom(0:max(x),size=ceiling(mle_ori[1]),prob=mle_ori[2])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)  ##Line 60-62 are taken from R code of function disc_ks_test in gsgeneral package
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           nb.mle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],lowerbound,upperbound)  ##use the mle of original data as initial values
	   mle_new=t(mle_new)
    }else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,r=mle_new[1,j],p=mle_new[2,j],distri='nb',type='general')

	temp=list(r=r,p=p,alpha1=NULL,alpha2=NULL,n=NULL)
 }

##--------for bb distribution------
 if(distri=='bb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bb.mle(x,n,alpha1,alpha2,lowerbound,upperbound)  ##mle using the original data
    probs_ori=extraDistr::pbbinom(0:max(x), size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)  ##Line 83-85 are taken from R code of function disc_ks_test in gsgeneral package
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	   bb.mle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],lowerbound,upperbound)  ##use the mle of original data as initial values
	   mle_new=t(mle_new)
    }else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,n=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bb',type='general')

	temp=list(r=NULL,p=NULL,alpha1=alpha1,alpha2=alpha2,n=n)
 }

##--------for bnb distribution------
 if(distri=='bnb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bnb.mle(x,r,alpha1,alpha2,lowerbound,upperbound)  ##mle using the original data
    probs_ori=extraDistr::pbnbinom(0:max(x),size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)  ##Line 106-108 are taken from R code of function disc_ks_test in gsgeneral package
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	  	  bnb.mle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],lowerbound,upperbound)  ##use the mle of original data as initial values
       mle_new=t(mle_new)
    }else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,r=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bnb',type='general')

    temp=list(r=r,p=NULL,alpha1=alpha1,alpha2=alpha2,n=NULL)
 }

 ##--------for zip distribution------
 if(distri=='zip'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=poisson.zihmle(x,type='zi',lowerbound,upperbound)  ##mle using the original data
    probs_ori=mle_ori[2]+(1-mle_ori[2])*stats::ppois(0:max(x),lambda=mle_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages=c('iZID','rootSolve')) %do%
	           poisson.zihmle(sample(x, size=N, replace=T),type='zi',lowerbound,upperbound)  ##use the mle of original data as initial values
	   mle_new=t(mle_new)
    }else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,lambda=mle_new[1,j],distri='poisson',type='zi',phi=mle_new[2,j])

    temp=list(r=NULL,p=NULL,alpha1=NULL,alpha2=NULL,n=NULL)
 }

  ##--------for zinb distribution------
 if(distri=='zinb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=nb.zihmle(x,r,p,type='zi',lowerbound,upperbound)  ##mle using the original data
    probs_ori=mle_ori[3]+(1-mle_ori[3])*stats::pnbinom(0:max(x),size=ceiling(mle_ori[1]),prob=mle_ori[2])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           nb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],type='zi',lowerbound,upperbound)  ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,r=mle_new[1,j],p=mle_new[2,j],distri='nb',type='zi',phi=mle_new[3,j])

    temp=list(r=r,p=p,alpha1=NULL,alpha2=NULL,n=NULL)
 }

   ##--------for zibb distribution------
 if(distri=='zibb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bb.zihmle(x,n,alpha1,alpha2,type='zi',lowerbound,upperbound)  ##mle using the original data
    probs_ori=mle_ori[4]+(1-mle_ori[4])*extraDistr::pbbinom(0:max(x), size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           bb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],type='zi',lowerbound,upperbound) ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,n=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bb',type='zi',phi=mle_new[4,j])

	temp=list(r=NULL,p=NULL,alpha1=alpha1,alpha2=alpha2,n=n)
 }

    ##--------for zibnb distribution------
 if(distri=='zibnb'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bnb.zihmle(x,r,alpha1,alpha2,type='zi',lowerbound,upperbound)  ##mle using the original data
    probs_ori=mle_ori[4]+(1-mle_ori[4])*extraDistr::pbnbinom(0:max(x),size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           bnb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],type='zi',lowerbound,upperbound) ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,r=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bnb',type='zi',phi=mle_new[4,j])

    temp=list(r=r,p=NULL,alpha1=alpha1,alpha2=alpha2,n=NULL)
 }

    ##--------for ph distribution------
 if(distri=='ph'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=poisson.zihmle(x,type='h',lowerbound,upperbound)  ##mle using the original data
    probs_ori=stats::ppois(0:max(x),lambda=mle_ori[1])
	probs_ori=mle_ori[2]+(1-mle_ori[2])*(probs_ori-probs_ori[1])/(1-probs_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages=c('iZID','rootSolve')) %do%
	           poisson.zihmle(sample(x, size=N, replace=T),type='h',lowerbound,upperbound) ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,lambda=mle_new[1,j],distri='poisson',type='h',phi=mle_new[2,j])

    temp=list(r=NULL,p=NULL,alpha1=NULL,alpha2=NULL,n=NULL)
 }

    ##--------for nbh distribution------
 if(distri=='nbh'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=nb.zihmle(x,r,p,type='h',lowerbound,upperbound)  ##mle using the original data
    probs_ori=stats::pnbinom(0:max(x),size=ceiling(mle_ori[1]),prob=mle_ori[2])
	probs_ori=mle_ori[3]+(1-mle_ori[3])*(probs_ori-probs_ori[1])/(1-probs_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           nb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],type='h',lowerbound,upperbound)##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages='iZID') %do%
	   general.ks(N,r=mle_new[1,j],p=mle_new[2,j],distri='nb',type='h',phi=mle_new[3,j])

    temp=list(r=r,p=p,alpha1=NULL,alpha2=NULL,n=NULL)
 }

     ##--------for bbh distribution------
 if(distri=='bbh'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bb.zihmle(x,n,alpha1,alpha2,type='h',lowerbound,upperbound)  ##mle using the original data
    probs_ori=extraDistr::pbbinom(0:max(x), size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
	probs_ori=mle_ori[4]+(1-mle_ori[4])*(probs_ori-probs_ori[1])/(1-probs_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           bb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],type='h',lowerbound,upperbound) ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,n=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bb',type='h',phi=mle_new[4,j])

    temp=list(r=NULL,p=NULL,alpha1=alpha1,alpha2=alpha2,n=n)
 }

      ##--------for bnbh distribution------
 if(distri=='bnbh'){
    ##------calculate KS statistic based on the original data and its mle of theta-----
    mle_ori=bnb.zihmle(x,r,alpha1,alpha2,type='h',lowerbound,upperbound)  ##mle using the original data
    probs_ori=extraDistr::pbnbinom(0:max(x),size=ceiling(mle_ori[1]),alpha=mle_ori[2],beta=mle_ori[3])
	probs_ori=mle_ori[4]+(1-mle_ori[4])*(probs_ori-probs_ori[1])/(1-probs_ori[1])
    step_ori=stats::stepfun(0:max(x),c(0,probs_ori))  ##step function for the original samples
    z=stats::knots(step_ori)
    dev=c(0, stats::ecdf(x)(z)-step_ori(z))
    Dn_ori=max(abs(dev))  ##KS statistic for the original samples, D_n in Algorithm 1.

	##--------obtain the bootstrapped or simulated mle of nsim iterations
    if(bootstrap){#if bootstrap=TRUE, obtain bootstrapped mles
	   mle_new=foreach::foreach(j=1:nsim,.combine=rbind,.packages='iZID') %do%
	           bnb.zihmle(sample(x, size=N, replace=T),mle_ori[1],mle_ori[2],mle_ori[3],type='h',lowerbound,upperbound) ##use the mle of original data as initial values
       mle_new=t(mle_new)
	}else{ ##when simulate new samples, uses the mle of original samples
       mle_new=matrix(rep(mle_ori,nsim),ncol=nsim)
	}

	##--------obtain the bootstrapped or simulated KS statistics---------
	D2=foreach::foreach(j=1:nsim,.combine=c,.packages=c('iZID','extraDistr')) %do%
	   general.ks(N,r=mle_new[1,j],alpha1=mle_new[2,j],alpha2=mle_new[3,j],distri='bnb',type='h',phi=mle_new[4,j])

    temp=list(r=r,p=NULL,alpha1=alpha1,alpha2=alpha2,n=NULL)
 }

 ##--------calculate the bootstrapped or simulated p value
 parallel::stopCluster(cl)
 pvalue=sum(D2>Dn_ori)/nsim
 ans=list(x=x,nsim=nsim,bootstrap=bootstrap,distri=distri,lowerbound=lowerbound,upperbound=upperbound,mle_new=mle_new,mle_ori=mle_ori,pvalue=pvalue,
 N=N)
 ans=c(ans,temp)
 class(ans)='dis.kstest'
 return(ans)
}

##examples
#x=sample.zi(N=2000,phi=0.2,distri='Poisson',lambda=5)+0.3
#dis.kstest(x,nsim=100,bootstrap=TRUE,distri='Poisson',lambda=3)