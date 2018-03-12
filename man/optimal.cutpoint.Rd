\name{optimal.cutpoint}
\alias{optimal.cutpoint}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Find the First Cutpoint and its P Value
}
\description{
This function finds the first optimal cutpoint for the time-varying regression effects based on the maximized score statistics
and calculates p-value based on a formula from Davies (1987) and O'Quigley and Pessione (1991). This is for depth 1 only.
}
\usage{
optimal.cutpoint(survtime, survstatus, x, method = "breslow", acpf = 10, 
            iter.max = 20, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survtime}{
	survival time/ follow up time of subjects
}
  \item{survstatus}{
	survival status of subjects. 0 for censored and 1 for an event
}
  \item{x}{
	a data frame of covariates. In case of a single covariate, use \code{[,,drop =F]} to keep the data frame structure
}
  \item{method}{
	argument for coxph function. Default is 'breslow'. See  \code{\link[survival]{coxph}} for more details.
}
  \item{acpf}{
	The search for the optimal cutpoint starts from the ((acpf/2)+1)th event until the (k - (acpf/2))th event, where k is the total number of events. Default is 10.	
}
  \item{iter.max}{
	the maximum number of iteration in coxph; default is 20. See  \code{\link[survival]{coxph}} for more details.
}
  \item{eps}{
	argument for coxph function; default is 0.000001. See  \code{\link[survival]{coxph}} for more details.
}
}
\details{
\code{optimal.cutpoint} takes in survival time, survival status, and covariates to find the first optimal cutpoint. 

 Currently, data need to be arranged in descending order of time and with no missing. 
}
\value{
\code{optimal.cutpoint} returns the following information:

\item{breakpt}{optimal cutpoint}
\item{scoretest}{Maximum score associated with the optimal cut point}
\item{summary}{3 output from \code{coxph} fitted with 1) entire data, 2) data before the optimal cutpoint, and 3) data after the optimal cutpoint.}
\item{pvalue}{p-value to test the existance of a change point against none}
}
\references{
Davies, R. (1987). Hypothesis Testing when a Nuisance Parameter is Present Only Under the Alternatives. Biometrika,  74(1), 33-43.

O'Quigley, J., and Pessione, F. (1991). The Problem of a Covariate-Time Qualitative Interaction in a Survival Study. Biometrics, 47(1), 101-115.
}

\examples{
##Call in alcohol data set
data('alcohol')
require(survival)

coxtree <- optimal.cutpoint(alcohol[,'time'], alcohol[,'event'], 
                      x = alcohol[,'alc', drop = FALSE])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

