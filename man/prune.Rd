\name{prune}
\alias{prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Prune Using the Score Statistic
}
\description{
This function merges over-segmented intervals to create optimally pruned subtrees.
}
\usage{
prune(fulltree)
}
\arguments{
  \item{fulltree}{
  output from \code{output.coxphout}
}
}
\details{
\code{prune} uses the CART algorithm and -log (partial likelihood) as cost to find the optimally pruned subtrees.
}
\value{
\code{prune} returns a matrix with the following columns, where each row is an optimally pruned subtree:
\item{K}{subtrees number 1, 2, etc. Tree #1 is the full tree}
\item{N[1]}{Number of terminal nodes}
\item{alpha}{penalty parameter corresponding to the subtree}
\item{S[1]}{-log(partial likelihood) of the subtree}
\item{pruneoff}{Node that was removed from the previous larger subtree to obtain the current subtree}
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

