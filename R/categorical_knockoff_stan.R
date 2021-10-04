#' Use Bayesian method to generate knockoff copy for data with categorical covariates.
#'
#' @param dataX Predictor dataframe, rows correspond to observations, columns correspond to covariates.
#' @param burnin The number of MCMC burnin for stan
#' @return A dataframe with generated knockoff copy
#'
#' @importFrom rstan sampling
#' @importFrom rstan extract
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack rwish
#' @importFrom stats rmultinom
#'
#' @export categorical_knockoff_stan

categorical_knockoff_stan<-function(dataX,burnin=1000){
  if(is.null(colnames(dataX))){
    stop("The matrix of predictors must have column names.")
  }
  colnamesX<-colnames(dataX)
  N <- nrow(dataX)
  factorIndi<- sapply(seq(1, ncol(dataX)), function(x){nlevels(dataX[,x])})
  factorIndi[factorIndi==0]<-1
  contNum<-sum(factorIndi==1)
  cateNum<-sum(factorIndi>1)

  if(cateNum==0){
    stop("This function requires the data frame provided must have a categorical covariate.")
  }

  cateVar<-colnames(dataX)[factorIndi>1]
  contVar<-colnames(dataX)[factorIndi==1]

  if(contNum>0){
    contX <- scale(dataX[,contVar])
    contXScale<-attr(contX, "scaled:scale")
    contXCenter<-attr(contX, "scaled:center")
    contX<-as.matrix(contX)
    contNumReal<-contNum
  }else{
    contNumReal<-0
    contNum<-10
    contX<-matrix(0,nrow=nrow(dataX),ncol=contNum)
  }

  cateX0<-dataX[,cateVar]
  cateX<-matrix(unlist(lapply(1:cateNum, function(x){as.numeric(cateX0[,x])})),ncol=cateNum)

  cateDummyNum<-sum(factorIndi[factorIndi>1])
  lengthZ<-contNum+cateDummyNum-cateNum

  cumCateDummyNum<-c(0,cumsum(factorIndi[factorIndi>1]))
  cumCateDummyNumInter<-c(0,cumsum(factorIndi[factorIndi>1]-1))
  mu0<-rep(0,lengthZ)
  tau0<-diag(rep(0.1/lengthZ,lengthZ))
  Z<-mu0
  Sigma<-MCMCpack::rwish(lengthZ,tau0)
  initsfun <-function() list(Z=Z, tau=1,Sigma=Sigma)

  data <- list(N=N, contX=contX, cateX= cateX, contNum=contNum, contNumReal=contNumReal,
               cateNum=cateNum, cateDummyNum=cateDummyNum,
               cumCateDummyNum=cumCateDummyNum, cumCateDummyNumInter=cumCateDummyNumInter,
               mu0=mu0, tau0=tau0,lengthZ=lengthZ)
  fit_rstan <- rstan::sampling(stanmodels$categoricalKO,
    data = data,
    init=initsfun,
    chains=1,
    warmup=burnin,
    iter=N+burnin
  )

  ZMCMC<-extract(fit_rstan)$Z
  tauMCMC<-extract(fit_rstan)$tau
  if(contNumReal>0){
    contMCMC<- ZMCMC[,1:contNum]
    contSample<-matrix(unlist(lapply(1:N,function(x) MASS::mvrnorm(1,contMCMC[x,],diag(rep(tauMCMC[x],contNum))))),
                       ncol=contNum,byrow=TRUE)
    contSample<-contSample*matrix(rep(contXScale,each=N),nrow=N)
    contSample<-contSample+matrix(rep(contXCenter,each=N),nrow=N)
  }
  cateMCMC<- ZMCMC[,(contNumReal+1):lengthZ]
  cateSample<-matrix(NA,nrow=N,ncol=cateNum)
  for(i in 1:cateNum){
    for(n in 1:N){
      temp<-c(0,cateMCMC[n,(1+cumCateDummyNumInter[i]):cumCateDummyNumInter[i+1]])
      temp<-exp(temp)/sum(exp(temp))
      cateSample[n,i]<- which(rmultinom(1, 1, temp)>0)
    }
    cateSample[,i]<-as.factor(levels(cateX0[,i])[cateSample[,i]])
  }
  if(contNumReal>0){
    knockoffX<-cbind(contSample,cateSample)
  }else{
    knockoffX<-cateSample
  }

  colnames(knockoffX)<-c(contVar,cateVar)
  knockoffX<-knockoffX[,colnamesX]
  knockoffX
}
