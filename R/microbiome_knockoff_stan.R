#' Produce knockoff copy for microbiome data
#'
#' @param microbiomeCount A matrix if microbiome sequence counts, rows correspond to observations, columns correspond to taxa.
#' @param burnin The number of MCMC burnin for stan
#' @return knockoffCopy - Generate knockoff copy for microbiome data
#'
#' @importFrom rstan sampling
#' @importFrom rstan extract
#' @importFrom stats rmultinom

#' @export microbiome_knockoff_stan

microbiome_knockoff_stan<-function(microbiomeCount,burnin=1000){
  if(max(apply(microbiomeCount,1,sum))<=1){
    stop("Please provide the matrix of microbiome sequence counts not matrix of abundance.")
  }
  N <- nrow(microbiomeCount)
  lengthZ<-ncol(microbiomeCount)-1
  mu0<-rep(0,lengthZ)
  tau0<-diag(rep(1/lengthZ,lengthZ))
  Z<-mu0
  Sigma<-diag(lengthZ)
  initsfun <-function() list(Z=Z, Sigma=Sigma)
  data <- list(N=N, Xvec=t(microbiomeCount), mu0=mu0, tau0=tau0,lengthZ=lengthZ)
  fit_rstan <- rstan::sampling(stanmodels$microbiomeKO,
    data = data,
    init=initsfun,
    chains=1,
    warmup=burnin,
    iter=burnin+N,
    pars=c("Z")
  )

  ZMCMC<-extract(fit_rstan)$Z
  probMCMC<-exp(cbind(rep(0,nrow(ZMCMC)),ZMCMC))
  probMCMC<-apply(probMCMC,1,function(x) x/sum(x)) # probMCMC num of taxa * num of obs
  knockoffCopy<-probMCMC
  seqDep<-apply(microbiomeCount,1,sum)
  for(i in 1:ncol(probMCMC)){
    knockoffCopy[,i]<-rmultinom(1,seqDep[i],probMCMC[,i])
  }
  knockoffCopy
}

