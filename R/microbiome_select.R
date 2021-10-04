#' Run feature selection for microbiome data
#'
#' @param microbiomeCount A matrix of the original microbiome sequence counts, rows correspond to observations, columns correspond to taxa. Users must provide column names.
#' @param knockoffCopy A matrix of knockoff microbiome sequence counts.
#' @param y A continuous outcome
#' @param Tree A matrix reflecting the phylogenetic relationships between taxa. The rows correspond to the tip node of the phylogenetic tree, while the columns correspond to the tip and internal nodes of the phylogenetic tree. For each row(tip node), the entries that correspond  to the phylogenentic tree node (including the tip node)that locate on the path from the root of the tree to the tip node are 1s, while other entries are 0s.
#' @param fdr Target fdr
#' @param offset The offset for knockoff+
#' @return A list with selected taxa or tree node (if phylogenetic tree is provided)
#' \itemize{
#'  \item selected - The taxa selected
#'  \item nodeSelected - If phylogenetic tree is provided, the node selected will also be provided
#' }

#' @importFrom reticulate import
#' @importFrom knockoff knockoff.threshold
#' @importFrom stats sd
#' @export microbiome_select
#'
#' @references Bien, J., Yan, X., Simpson, L. et al. Tree-aggregated predictive modeling of microbiome data. Sci Rep 11, 14505 (2021). \url{https://doi.org/10.1038/s41598-021-93645-3}
#' @examples
#' \dontrun{
#' library(BayesianKnockoff)
#' data("HIVdata")
#' knockoffCopy<-microbiome_knockoff_stan(HIVdata$microbiomeCount) #may take a few hours
#' temp1<-microbiome_select(HIVdata$microbiomeCount,knockoffCopy,HIVdata$y)
#' temp2<-microbiome_select(HIVdata$microbiomeCount,knockoffCopy,HIVdata$y,HIVdata$Tree)
#' }

microbiome_select<-function(microbiomeCount,knockoffCopy,y,Tree=NULL,fdr=0.2,offset=0){
  if(is.null(colnames(microbiomeCount))){
    stop("The matrix of microbiome sequence counts must have column names.")
  }
  if(max(apply(microbiomeCount,1,sum))<=1){
    stop("Please provide the matrix of microbiome sequence counts not matrix of abundance.")
  }

  zko<-t(knockoffCopy)
  zko[zko == 0] <- 0.5
  zko <- zko/rowSums(zko)
  zko <- log(zko)

  microbiomeCount[microbiomeCount == 0] <- 0.5
  microbiomeCount <- microbiomeCount/rowSums(microbiomeCount)
  ztr <- log(microbiomeCount)

  classo <- reticulate::import("classo", delay_load = TRUE)
  #reticulate::py_install("c-lasso", method = "auto", conda = "auto", pip = TRUE)

  Z<-cbind(ztr,zko)
  y<-scale(y)
  if(is.null(Tree)){
    Tree<-diag(ncol(microbiomeCount))
  }

  ZeroPad<-matrix(0,nrow=nrow(Tree),ncol=ncol(Tree))
  Adouble<-rbind(cbind(Tree,ZeroPad),cbind(ZeroPad,Tree))
  Z <- as.matrix(Z %*% Adouble)
  Zclasso<-scale(Z, scale=apply(Z, 2, sd)*sqrt(nrow(Z)-1))
  fac <- 1/attr(Zclasso, "scaled:scale")
  C<-rbind(c(apply(Tree,2,sum)*fac[1:ncol(Tree)],rep(0,ncol(Tree))),
           c(rep(0,ncol(Tree)),apply(Tree,2,sum)*fac[(1+ncol(Tree)):ncol(Zclasso)]))
  nlambda<-2000
  lam.max <- max(abs(crossprod(Zclasso, y)))
  fraclist <- lam.max*exp(seq(0, log(1/nlambda), length=nlambda))
  # set up CLASSO problem:

  prob <- classo$classo_problem(X = Zclasso,
                                C = C,
                                y = array(y))
  prob$formulation$classification <- FALSE
  prob$formulation$concomitant <- FALSE
  prob$formulation$huber <- FALSE
  prob$model_selection$PATH <- TRUE
  prob$model_selection$CV <- FALSE
  prob$model_selection$LAMfixed <- FALSE
  prob$model_selection$StabSel <- FALSE
  prob$model_selection$PATHparameters$lambdas <- fraclist
  prob$solve()
  gamma<- t(prob$solution$PATH$BETAS)
  beta <- Adouble %*% gamma
  lambdas <- matrix(NA, nrow=nrow(gamma)/2, ncol=2)
  for (i in 1:nrow(lambdas)) {
    lambdas[i, ] <- c(max(fraclist[which(gamma[i,]!=0)]),
                      max(fraclist[which(gamma[i+nrow(lambdas),]!=0)]))
  }
  W <- apply(lambdas, 1, function(x) max(x) * sign(x[1]-x[2]))
  W[is.nan(W)]<-0
  threshold<-knockoff::knockoff.threshold(W, fdr=fdr, offset=offset)
  gammaSelect<-rep(0,nrow(lambdas))
  gammaSelect[which(W>threshold)]<-1
  betaSelect<-Tree %*%gammaSelect
  tipselected<-colnames(microbiomeCount)[which(betaSelect>0)]
  print(tipselected)
  if(is.null(colnames(Tree))&(!identical(Tree,diag(ncol(microbiomeCount))))){
    colnames(Tree)<-paste("node",1:ncol(Tree))
  }
  nodeSelected<-colnames(Tree)[which(W>threshold)]
  if(identical(as.matrix(Tree),diag(ncol(microbiomeCount)))){
    microbiomeKOResult<-list(selected=tipselected)
  }else{
    microbiomeKOResult<-list(selected=tipselected,nodeSelected=nodeSelected)
  }
  microbiomeKOResult
}

