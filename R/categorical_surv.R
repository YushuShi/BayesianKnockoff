#' Select features for data with categorical covariates and survival outcome using knockoff filter with group lasso.
#'
#' @param dataX The original data frame, where rows correspond to observations, while columns correspond to predictors. Categorical covariates need to be factors.
#' @param knockoffCopy The Knockoff copy containing categorical covariates. Can be produced by the 'categorical_knockoff_stan' of this package.
#' @param y The 'Surv' type right-censored survival outcome
#' @param fdr The target false discovery rate
#' @param offset Offset for the group knockoff
#'
#' @return The name of the predictors selected under the target FDR
#' @importFrom knockoff knockoff.threshold
#' @importFrom grpreg grpsurv
#' @importFrom survival Surv
#' @importFrom methods hasArg
#' @importFrom fastDummies dummy_cols
#' @export categorical_surv
#'
#' @examples
#' library(BayesianKnockoff)
#' library(survival)
#' set.seed(1)
#' data0<-data.frame(matrix(sample(1:3,1000,replace = TRUE),ncol=10))
#' data1<-apply(data0,2,as.character)
#' data1 <- as.data.frame(unclass(data1),stringsAsFactors = TRUE)
#' data2<-data.frame(matrix(rnorm(1000),ncol=10))
#' colnames(data2)<-paste0("X",11:20)
#' dataX<-cbind(data1,data2)
#' knockoffCopy<-categorical_knockoff_stan(dataX)
#'
#' time0<-(data0$X5-2)*10-(data0$X9-2)*10+data2$X11*10+data2$X12*10
#' time0<-time0-min(time0)+1
#' time1<-rexp(nrow(dataX),0.01)
#' cens<-time0<time1
#' time<-ifelse(cens,time0,time1)
#' tempsurv<-categorical_surv(dataX,knockoffCopy,Surv(time,cens))


categorical_surv<- function(dataX,knockoffCopy,y,fdr=0.2,offset=0){
  N <- nrow(dataX)
  factorIndi<- sapply(seq(1, ncol(dataX)), function(x){nlevels(dataX[,x])})
  cateVar<-colnames(dataX)[factorIndi>1]
  contVar<-colnames(dataX)[factorIndi==0]
  dataC<-rbind(dataX,knockoffCopy)
  dataC<-fastDummies::dummy_cols(dataC, select_columns=cateVar,
                                 remove_selected_columns = TRUE, remove_most_frequent_dummy = TRUE)
  combX<-cbind(dataC[1:N,],dataC[(1+N):(2*N),])
  group<-c(1:length(contVar),rep(1:length(cateVar),times=(factorIndi[factorIndi>1]-1))+length(contVar))
  groupknockoffCopy<-paste0("knockoffCopy",group)

  groupComb<-c(group,groupknockoffCopy)
  grpregFit <- grpreg::grpsurv(as.matrix(combX), y, group=groupComb, penalty="grLasso")
  betas <- grpregFit$beta
  Lambdas <- matrix(NA, nrow=ncol(dataX), ncol=2)
  lambda<-as.numeric(colnames(betas))
  for (i in 1:ncol(dataX)) {
    Lambdas[i, ] <- c(max(lambda[apply(matrix(betas[which(group==i),],ncol=ncol(betas)),2,function(x) sum(x != 0))>0]),
                      max(lambda[apply(matrix(betas[which(group==i)+length(group),],ncol=ncol(betas)),2,function(x) sum(x != 0))>0]))
  }
  W <- apply(Lambdas, 1, function(x) max(x) * sign(x[1]-x[2]))
  W[is.nan(W)]<-0
  threshold.group<-knockoff::knockoff.threshold(W, fdr=fdr, offset=offset)
  selected<-c(contVar,cateVar)[which(W>=threshold.group)]
  selected
}
