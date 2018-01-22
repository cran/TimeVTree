\name{bootstrap}
\alias{bootstrap}
\alias{boot.coxfit}
\alias{free}

\title{
Bootstrap to Correct for Over-optimism due to Adaptive Splitting
}
\description{
	This function is used to obtain the bias-corrected cost.
	One may select the final subtree with the lowest bootstrap estimated cost, with or without the additional AIC/BIC as in Xu and Adak (2002).
}
\usage{
bootstrap(B = 20, nodetree, subtrees, survtime, survstatus, 
          x, D = 4, minfail = 30, alphac = 2)
}

\arguments{
  \item{B}{
Number of bootstrap samples. Default is 20.
}
  \item{nodetree}{
Full grown tree with original data. Output from \code{\link{output.coxphout}}
}
  \item{subtrees}{
Pruned subtrees with original data. Output from \code{\link{prune}}
}
  \item{survtime}{
	survival time/follow up time of subjects
}
  \item{survstatus}{
	survival status of subjects.
}
  \item{x}{
	a data frame of covariates. In case of single covariate, use [,,drop =FALSE] to keep the data frame structure
}
  \item{D}{
	maximum depth the tree will grow. Default depth is 4.
}
  \item{minfail}{
	minimum number of unique event required in each block. Default is 10
}
  \item{alphac}{
	Predetermined penalty parameter
}
}
\details{
The implemented cost here is the negative log partial likelihood. Each 
bootstrap sample is used to grow a full tree and then pruned to obtain the set of subtrees. The bias is estimated by the average of the differences between the cost of a bootstrapped subtree itself and the cost of sending the original data down the bootstrapped subtree. The bias-corrected cost is then obtained by subtracting this bias from the original cost.  
Predetermined penalty parameter can be used to account for the dimension of covariates, via Akaike information criteria (AIC), Schwarz Bayesian information criteria (BIC), or the 0.95 quantile of the chi-square distribution.
}
\value{
\item{bcoef }{coefficient values from each bootstrap sample}
\item{btree }{Tree related information from each bootstrapped sample. Types of information are the same as the ones from \code{output.coxphout}}
\item{bomega }{Bias at each subtree for each bootstrapped data, the average of which gives the overall bootstrap estimated bias}
\item{bootcost}{cost based on the bootstrapped data}
\item{ori.boot}{negative log partial likelihood of the original data fitted to the model given by bootstrapped data}
}
\references{
Xu, R. and Adak, S. (2002), Survival Analysis with Time-Varying Regression Effects Using a Tree-Based Approach. Biometrics, 58: 305-315.
}

\examples{
\dontrun{
data('alcohol')
require(survival)

coxtree <- coxph.tree(alcohol[,'time'], alcohol[,'event'], 
                      x = alcohol[,'alc', drop = FALSE], D = 4)
nodetree <- output.coxphout(coxtree)

subtrees <- prune(nodetree)

#This function requires output from output.coxphout, prune, and the original data set.

store.mult.cont <- bootstrap(B=20, nodetree, subtrees, alcohol[,'time'],
                                alcohol[,'event'], x = alcohol[,'alc', drop = FALSE], 
                                D=4,minfail=20, alphac=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
