\name{final.tree}
\alias{final.tree}

\title{
Finding the Final Tree After Bootstrap
}
\description{
\code{final.tree} uses bias-corrected costs obtained from \code{bootstrap} function and the predetermined penalty parameter to find the optimal tree from the set of subtrees.
}
\usage{
final.tree(nodetree=nodetree, subtrees=subtrees, omega, alphac=2)
}

\arguments{
  \item{nodetree}{
Fully grown tree from the original data. Output from \code{\link{output.coxphout}}
}
  \item{subtrees}{
Pruned subtrees from the original data. Output from \code{\link{prune}}
}
  \item{omega}{
	Bias (i.e. third index of the output) from \code{\link{bootstrap}}. Look at the value section of \code{\link{bootstrap}} for more information.
}
  \item{alphac}{
	Predetermined penalty parameter
}
}
\details{
	\code{final.tree} is part of the \code{\link{bootstrap}} function but can be used  to try different penalty parameters without re-running \code{bootstrap}. 
}
\value{
\item{subtree }{output from \code{prune} with an additional column 'cost' that contains bootstrap estimate of each subtree}
\item{final }{A tree with lowest cost value after applying predetermined penalty}
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

store.mult.cont <- bootstrap(B=20, nodetree, subtrees, alcohol[,'time'],
                                alcohol[,'event'], x = alcohol[,'alc', drop = FALSE], 
                                D=4,minfail=20, alphac=2)
                                
Balph <- 0.5 * 2 * log(nrow(alcohol))                                
final.tree <- final.tree(nodetree, subtrees, store.mult.cont[[3]], alphac= Balph)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
