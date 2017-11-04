#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
library(mvtnorm)
require(mvtnorm)
agents <- as.numeric(args[1])
numtypesa <- as.numeric(args[2])
numtypesb <- as.numeric(args[3])
discrete <- matrix(nrow=numtypesa, ncol=numtypesb)
mean <- rep(0, 2)
corrr <- diag(2)
corrr[1,1] <- 1
corrr[1,2] <- as.numeric(args[4])
corrr[2,1] <- as.numeric(args[4])
corrr[2,2] <- 1
prob <- pmvnorm(lower=-Inf, upper=Inf, mean, sigma=corrr)
edges = 2
rowlen <- (2 * edges)/numtypesa
collen <- (2 * edges)/numtypesb
lowvec <- c(-Inf, -Inf)
#print(pmvnorm(lower=-Inf, upper=c(Inf,0), mean, sigma=corrr))
for (i in 1:numtypesa)
{
	for (j in 1:numtypesb)
	{
		ay <- 2 - ((i - 1) * rowlen)
		ax <- -2 + ((j - 1) * collen)
		apoint <- c(ax, ay)
		bpoint <- c(ax + collen, ay)
		cpoint <- c(ax, ay - rowlen)
		dpoint <- c(ax + collen, ay - rowlen)
		acdf <- pmvnorm(lower=lowvec, upper=apoint, mean, corr=NULL, sigma=corrr)
		bcdf <- pmvnorm(lower=lowvec, upper=bpoint, mean, corr=NULL, sigma=corrr)
		ccdf <- pmvnorm(lower=lowvec, upper=cpoint, mean, corr=NULL, sigma=corrr)
		dcdf <- pmvnorm(lower=lowvec, upper=dpoint, mean, corr=NULL, sigma=corrr)
		discrete[numtypesa - i + 1, j] <- bcdf - acdf - dcdf + ccdf
	}
}
write.table(discrete, file=args[5], row.names=FALSE, col.names=FALSE)
#print(discrete)
