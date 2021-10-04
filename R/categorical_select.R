#' Select features for data with categorical covariates using knockoff filter with group lasso.
#'
#' @param dataX The original data frame, where rows correspond to observations, while columns correspond to predictors. Categorical covariates need to be factors.
#' @param knockoffCopy The Knockoff copy containing categorical covariates. Can be produced by the 'categorical_knockoff_stan' of this package.
#' @param y The outcome vector. Can be continous, count or binary.
#' @param fdr The target false discovery rate
#' @param offset Offset for the group knockoff
#' @param family The family of the response variable, can be 'binomial', 'poisson' or 'gaussian'.
#' @param nlambda Number of penalty parameter used for calculating the knockoff statistics
#' @param lambda A user supplied sequence of lambda values. Typically, this is left unspecified, and the function automatically computes a grid of lambda values that ranges uniformly on the log scale over the relevant range of lambda values.
#'
#' @return The name of the predictors selected under the target FDR
#' @importFrom knockoff knockoff.threshold
#' @importFrom grpreg grpreg
#' @importFrom methods hasArg
#' @importFrom fastDummies dummy_cols
#' @export categorical_select
#'
#' @examples
#' library(BayesianKnockoff)
#' library(grpreg)
#' set.seed(1)
#' data0<-data.frame(matrix(sample(1:3,1000,replace = TRUE),ncol=10))
#' data1<-apply(data0,2,as.character)
#' data1 <- as.data.frame(unclass(data1),stringsAsFactors = TRUE)
#' data2<-data.frame(matrix(rnorm(1000),ncol=10))
#' colnames(data2)<-paste0("X",11:20)
#' dataX<-cbind(data1,data2)
#' knockoffCopy<-categorical_knockoff_stan(dataX)
#'
#' y<-(data0$X5-2)*10-(data0$X9-2)*10+data2$X11*10+data2$X12*10
#' temp<-categorical_select(dataX,knockoffCopy,y)
#'
#' y2<-as.numeric(cut(y,breaks = quantile(y),include.lowest = TRUE))-1
#' temp2<-categorical_select(dataX,knockoffCopy,y2,family = "poisson")
#'
#' y3<-as.numeric(y>median(y))
#' temp3<-categorical_select(dataX,knockoffCopy,y3,family = "binomial")

#' @references Dai, R., Barber, R.. The knockoff filter for FDR control in group-sparse and multitask regression. (2016) \url{	arXiv:1602.03589}


categorical_select<- function(dataX,knockoffCopy,y,fdr=0.2,offset=0,family="gaussian",nlambda=500,lambda=NULL){
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

    combX <- scale(combX)
    combX[is.nan(combX)]<-0

    if (!methods::hasArg(lambda) ) { # if we do not have lambda vector provided
      if(identical(family, "gaussian") ) {
        if(!is.numeric(y)) {
          stop('Input y must be numeric.')
        }
        lambda_max = max(abs(t(combX) %*% y)) / N
        lambda_min = lambda_max / 2e3 #keep it same with the original
        k = (0:(nlambda-1)) / nlambda
        lambda = lambda_max * (lambda_min/lambda_max)^k
      }
      else {
        lambda = NULL
      }
    }
    groupComb<-c(group,groupknockoffCopy)
    if(identical(family, "gaussian") ) {
      grpregFit <- grpreg::grpreg(combX, y, group=groupComb, penalty="grLasso", family=family,lambda=lambda)
    }else{
      grpregFit <- grpreg::grpreg(combX, y, group=groupComb, penalty="grLasso", family=family, nlambda=nlambda)
    }

    betas <- grpregFit$beta[-1, ]
    Lambdas <- matrix(NA, nrow=ncol(dataX), ncol=2)
    for (i in 1:ncol(dataX)) {
      Lambdas[i, ] <- c(max(lambda[apply(matrix(betas[which(group==i),],ncol=nlambda),2,function(x) sum(x != 0))>0]),
      max(lambda[apply(matrix(betas[which(group==i)+length(group),],ncol=nlambda),2,function(x) sum(x != 0))>0]))
    }
    W <- apply(Lambdas, 1, function(x) max(x) * sign(x[1]-x[2]))
    W[is.nan(W)]<-0
    threshold.group<-knockoff::knockoff.threshold(W, fdr=fdr, offset=offset)
    selected<-c(contVar,cateVar)[which(W>=threshold.group)]
    selected
}
