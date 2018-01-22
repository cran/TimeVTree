\name{output.coxphout}
\alias{output.coxphout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of coxph.tree output
}
\description{
This funtion organizes \code{coxph.tree} output into a format that can be used as an input for \code{prune}, \code{plot_coxphtree}, and \code{mat.tvbeta}.
}
\usage{
output.coxphout(coxout) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coxout}{
	output from \code{coxph.tree}
}
}

\value{
\code{output.coxphout} returns a table with following columns.

\item{Depth}{Depth value specified in the argument}
\item{Block}{Time intervals present at each depth}
\item{Node}{Unique number assigned to each block}
\item{Left}{Node of a block that was divided into the left side in the next depth}
\item{Right}{Node of a block that was divided into the right side in the next depth}
\item{Score}{Modified score statistic of each node}
\item{lkl}{Likelihood ratio value of each node}
\item{Start}{Starting time of the node}
\item{End}{Ending time of the node}
\item{# of Cases}{Number of observations in each node}
\item{# of Events}{Number of events in each node}
}
\references{
Xu, R. and Adak, S. (2002), Survival Analysis with Time-Varying Regression Effects Using a Tree-Based Approach. Biometrics, 58: 305-315.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
