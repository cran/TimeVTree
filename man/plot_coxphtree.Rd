\name{plot_coxphtree}
\alias{plot_coxphtree}
\title{
Plotting of Full Tree and Subtrees
}
\description{
This functin uses the full tree and subtrees (optional) to create visual outputs of the tree(s) and segments.
}
\usage{
plot_coxphtree(fulltree, subtrees = NULL, mm = 3, start = 0, pdf = FALSE, file.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fulltree}{
output of \code{output.coxphout}.
}
  \item{subtrees}{
(Optional) output of \code{prune}.
}
  \item{mm}{
Number of subtrees plot to be placed in one page. Default is 3
}
  \item{start}{
Sets starting point for segments. Useful if the minimum event time is far away from 0.
}
  \item{pdf}{
Do you want to export the plots in pdf format? Default is FALSE.
When set as FALSE, all plots need to be cleared before running this function to avoid 'Plot rendering error.'
}
  \item{file.name}{
Name for the pdf file output. 
}
}
\details{
\code{plot_coxphtree} takes an output from \code{output.coxphout} and creates treeplot and barplot showing blocks at each depth.
If an output from \code{prune} is also included in the argument, the function creates treeplot and barplot for each subtree.
In the barplot, end nodes are in dark blue color.
}

\references{
Xu, R. and Adak, S. (2002), Survival Analysis with Time-Varying Regression Effects Using a Tree-Based Approach. Biometrics, 58: 305-315.
}

\examples{
#This function requires output from output.coxphout and prune(optional)
data('alcohol')
require(survival)

coxtree <- coxph.tree(alcohol[,'time'], alcohol[,'event'], 
                      x = alcohol[,'alc', drop = FALSE], D = 4)
nodetree <- output.coxphout(coxtree)

subtrees <- prune(nodetree)

plot_coxphtree(nodetree, subtrees, start = 70, pdf = FALSE)
}

