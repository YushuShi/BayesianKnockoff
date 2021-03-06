% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorical_select.R
\name{categorical_select}
\alias{categorical_select}
\title{Select features for data with categorical covariates using knockoff filter with group lasso.}
\usage{
categorical_select(
  dataX,
  knockoffCopy,
  y,
  fdr = 0.2,
  offset = 0,
  family = "gaussian",
  nlambda = 500,
  lambda = NULL
)
}
\arguments{
\item{dataX}{The original data frame, where rows correspond to observations, while columns correspond to predictors. Categorical covariates need to be factors.}

\item{knockoffCopy}{The Knockoff copy containing categorical covariates. Can be produced by the 'categorical_knockoff_stan' of this package.}

\item{y}{The outcome vector. Can be continous, count or binary.}

\item{fdr}{The target false discovery rate}

\item{offset}{Offset for the group knockoff}

\item{family}{The family of the response variable, can be 'binomial', 'poisson' or 'gaussian'.}

\item{nlambda}{Number of penalty parameter used for calculating the knockoff statistics}

\item{lambda}{A user supplied sequence of lambda values. Typically, this is left unspecified, and the function automatically computes a grid of lambda values that ranges uniformly on the log scale over the relevant range of lambda values.}
}
\value{
The name of the predictors selected under the target FDR
}
\description{
Select features for data with categorical covariates using knockoff filter with group lasso.
}
\examples{
library(BayesianKnockoff)
library(grpreg)
set.seed(1)
data0<-data.frame(matrix(sample(1:3,1000,replace = TRUE),ncol=10))
data1<-apply(data0,2,as.character)
data1 <- as.data.frame(unclass(data1),stringsAsFactors = TRUE)
data2<-data.frame(matrix(rnorm(1000),ncol=10))
colnames(data2)<-paste0("X",11:20)
dataX<-cbind(data1,data2)
knockoffCopy<-categorical_knockoff_stan(dataX)

y<-(data0$X5-2)*10-(data0$X9-2)*10+data2$X11*10+data2$X12*10
temp<-categorical_select(dataX,knockoffCopy,y)

y2<-as.numeric(cut(y,breaks = quantile(y),include.lowest = TRUE))-1
temp2<-categorical_select(dataX,knockoffCopy,y2,family = "poisson")

y3<-as.numeric(y>median(y))
temp3<-categorical_select(dataX,knockoffCopy,y3,family = "binomial")
}
\references{
Dai, R., Barber, R.. The knockoff filter for FDR control in group-sparse and multitask regression. (2016) \url{	arXiv:1602.03589}
}
