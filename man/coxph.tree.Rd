\name{coxph.tree}
\alias{coxph.tree}
\alias{infmul}
\alias{rawscore}
\alias{xname}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Grow the Tree Using the Score Statistic
}
\description{
This funtion finds the optimal cutpoints for the time-varying regression effects and grows the 'full tree' using the score statistic.
}
\usage{
coxph.tree(survtime, survstatus, x, D = 3, method = "breslow", minfail = 10, 
            iter.max = 20, eps = 1e-06, type = 'mod')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survtime}{
	survival time/ follow up time of subjects
}
  \item{survstatus}{
	survival status of subjects. 0 for censored and 1 for event
}
  \item{x}{
	a data frame of covariates. In case of single covariate, use \code{[,,drop =F]} to keep the data frame structure
}
  \item{D}{
	maximum depth the tree will grow. Default depth is 3.
}
  \item{method}{
	argument for coxph function. Default is 'breslow'. See  \code{\link[survival]{coxph}} for more details.
}
  \item{minfail}{
	minimum number of unique events required in each block. Default is 10	
}
  \item{iter.max}{
	the maximum number of iteration in coxph; default is 20. See  \code{\link[survival]{coxph}} for more details.
}
  \item{eps}{
	argument for coxph function; default is 0.000001. See  \code{\link[survival]{coxph}} for more details.
}
  \item{type}{
	method to calculate the score statistic. Two options are available: 'mod' for the modified score statistic and 'ori' for the original score statistic.
	Default value is 'mod.' Modified score statistic is used in the bootstrap part
}
}
\details{
\code{coxph.tree} takes in survival time, survival status, and covariates to grow the full tree. 
It follows one of the stopping rules: 1) when the pre-specified depth is reached, or 2) the number of events in a node is less than a prespecified number, or 3) the maximized score statistic is less than a default value (0.0001).

 Currently data need to be arranged in descending order of time and with no missing. 
}
\value{
\code{coxph.tree} returns an object of class 'coxphtree.'

The function \code{\link{output.coxphout}} is used to obtain and print a summary of the result.

An object of class 'coxphtree' is a list containing the following components:
\item{D}{Depth value specified in the argument}
\item{coef}{coefficient values of predictors. First number represents depth and second number represents block number}
\item{lkl}{Likelihood ratio value of each node}
\item{breakpt}{Starting point of each node. Starting point of node at Depth= 0 to maximum Depth = D+1 is shown.}
\item{ntree}{Number of cases in each node}
\item{nevent}{Number of events in each node}
}
\references{
Xu, R. and Adak, S. (2002), Survival Analysis with Time-Varying Regression Effects Using a Tree-Based Approach. Biometrics, 58: 305-315.
}

\examples{
##Call in alcohol data set
data('alcohol')
require(survival)

coxtree <- coxph.tree(alcohol[,'time'], alcohol[,'event'], 
                      x = alcohol[,'alc', drop = FALSE], D = 4)

nodetree <- output.coxphout(coxtree)

subtrees <- prune(nodetree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

