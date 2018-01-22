\name{mat.tvbeta}
\alias{mat.tvbeta}
\title{
Beta coefficient estimate at each time point
}
\description{
Function that ouputs beta coefficient estimate of each covariate at each observation time point for a given tree, which can be used to plot the time-varying coefficients.
}
\usage{
mat.tvbeta(indx, fulltree, subtrees = NULL, survtime, survstatus, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indx}{
Index number of a subtree that needs to be analyzed
}
  \item{fulltree}{
output of \code{output.coxphout}.
}
  \item{subtrees}{
(Optional) output of \code{prune}. 
}
  \item{survtime}{
	survival time/ follow up time of subjects
}
  \item{survstatus}{
	survival status of subjects. 0 for alive and 1 for dead
}
  \item{x}{
	a data frame of covariates. In case of single covariate, use [,,drop =F] to keep the data frame structure
}
}

\value{
For each predictor, \code{mat.tvbeta} gives the coefficient values at each observation time for a given subtree.
The function outputs a matrix that can be used to plot the time-varying coefficient estimates over time.
The number of rows in the matrix is the # of observations and the number of columns is the product of the # of covariates and the # of specified subtrees.
}
\references{
Xu, R. and Adak, S. (2002), Survival Analysis with Time-Varying Regression Effects Using a Tree-Based Approach. Biometrics, 58: 305-315.
}

\examples{
#This function requires output from output.coxphout, prune, and the original data set.
data('alcohol')
require(survival)

coxtree <- coxph.tree(alcohol[,'time'], alcohol[,'event'], 
                      x = alcohol[,'alc', drop = FALSE], D = 4)
nodetree <- output.coxphout(coxtree)

subtrees <- prune(nodetree)

#creating matrix of beta coefficients at each event time point for all subtrees
k <- nrow(subtrees)
for (l in 1:k) {
    print(paste("Tree #",l))
    coeftmp <- mat.tvbeta(l,nodetree,subtrees,alcohol[,'time'], alcohol[,'event'],
                           x = data.frame(model.matrix(~alc, data=alcohol)[,-c(1), drop = FALSE]))
    if (l == 1) coef <- coeftmp
    if (l > 1) coef <- cbind(coef,coeftmp) 
  }


##Creating plot of all subtrees for each predictor:

p <- ncol(coef)/k #Number of variables
x = data.frame(model.matrix(~alc, data=alcohol)[,-c(1), drop = FALSE])
xnames <- xname(x)
xnames <- c('Alcohol 1', 'Alcohol 4')
#Subsetting data
coefnew <- data.frame(coef)
survtime <- alcohol[,'time']
#Setting desired depth (All the subtrees)
kk <- nrow(subtrees)
for (j in 1:p) {
  matplot(survtime,coefnew[,seq(from=j,to=kk*p,by=p)],type="l",lty=1:kk,col= (1:kk)+1
          ,xlab="Survival Time",ylab=" ")
  title(main=paste('all:', xnames[j]))
  legend('bottomleft', legend = paste('tree number', 1:kk), lty=1:kk,col= (1:kk)+1)
  }
  
##Creating a plot showing changes in coefficient of two predictors in full tree
#creating matrix of beta coefficients at each event time point for full tree
coeftmp <- mat.tvbeta(1,nodetree,subtrees,alcohol[,'time'], alcohol[,'event'],
                       x = data.frame(model.matrix(~alc, data=alcohol)[,-c(1), drop = FALSE]))
coefnew <- coeftmp
matplot(survtime,coefnew,type="l",lty=1:2,col= (1:2)+1,xlab="Survival Time",ylab=" ")
legend('bottomleft', legend = c("Alcohol 1", "Alcohol 4"), lty=1:2,col= (1:2)+1)  
}

