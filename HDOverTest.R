library(MASS)
library(glmnet)
library(matrixStats) 
library(dplyr)

HDOverTest_IV <- function(D, X, Z, Y) {
  n <- length(Y)
  px <- ncol(X); pz <- ncol(Z); p <- px + pz
  
  # gamma.hat = (gamma.x.hat, gamma.z.hat)
  V <- cbind(X, Z)
  gamma.lasso <- cv.glmnet(V, D, nfolds=10, intercept=FALSE)
  D.hat <- predict(gamma.lasso, V, s="lambda.min")
  
  # beta.hat, alpha.x.star.hat
  X.new <- cbind(D.hat, X)
  alpha.star.lasso <- cv.glmnet(X.new, Y, nfolds=10, intercept=FALSE)
  alpha.star.coef <- coef(alpha.star.lasso, s="lambda.min")[-1]
  beta.hat <- alpha.star.coef[1]
  
  # alpha.x.hat
  alpha.lasso <- cv.glmnet(X, Y-D*beta.hat, nfolds=10, intercept=FALSE)
  alpha.x.hat <- coef(alpha.lasso, s="lambda.min")[-1]

  # construction of test statistics
  Z.mat <- Z %*% t(Z)
  tr.SigZ2.hat <- (sum(Z.mat^2) - sum(diag(Z.mat^2)))/(n*(n-1))
  err.mat <- outer(as.numeric(Y-D*beta.hat-X%*%alpha.x.hat), as.numeric(Y-D*beta.hat-X%*%alpha.x.hat), "*")
  Tn <- (sum(err.mat*Z.mat) - sum(diag(err.mat*Z.mat)))/n
  sigma2.hat <- mean((Y-D*beta.hat-X%*%alpha.x.hat)^2)
  
  Sn <- Tn/(sigma2.hat*sqrt(2*tr.SigZ2.hat))
  pval <- 1 - pnorm(Sn)
  
  return(pval)
}

## Examples
library(MASS)
library(glmnet)
library(matrixStats) 
library(dplyr)

# generate data
gene_data <- function(n=100, px=100, pz=50, s.alphax=10, s.gammax=10, s.alphaz=10, s.gammaz=10,
                      varsigma=5, IV.strength=1, a0=0) {
  
  # the dimension of V=(X, Z)
  p <- px + pz
  
  # the coefficients
  alpha.x <- c(0.5^(0:(s.alphax-1)), rep(0, px-s.alphax))
  gamma.x <- c(0.6^(0:(s.gammax-1)), rep(0, px-s.gammax))
  
  # the validty of IV
  alpha.z <- c(rep(varsigma/sqrt(s.alphaz), s.alphaz), rep(0, pz-s.alphaz))
  
  # the strength of IV
  if(IV.strength == 1) {
    gamma.z <- c(rep(1, s.gammaz), rep(0, pz-s.gammaz))
  }
  if(IV.strength == 2) {
    gamma.z <- c(seq(0.1, 1, length.out=s.gammaz), rep(0, pz-s.gammaz))
  }
  
  beta.ture <- 0.5
  
  # the covariance matrix structure
  Sig <- toeplitz(0.5^seq(0, p-1))
  
  # covariates and instruments
  V <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=Sig)
  if (px == 0) {
    X <- NULL; Z <- V
  } else {
    X <- V[,1:px]; Z <- V[,-(1:px)]
  }
  
  # error terms: the correlation between epsilon and delta is not zero
  err.epsilon <- a0*Z[,1]*rnorm(n, 0, 1) + sqrt(1-a0^2)*rnorm(n, 0, 1)
  err.delta <- 0.5*err.epsilon + sqrt(1-0.5^2)*rnorm(n, 0, 1)
  
  # endogenous variable and response
  D <- X%*%gamma.x + Z%*%gamma.z + err.delta
  Y <- D*beta.ture + X%*%alpha.x + Z%*%alpha.z + err.epsilon
  
  return(list(D=D, X=X, Z=Z, Y=Y))
}

set.seed(2023)

# H0 + dense
out.data.null <- gene_data(n=200, px=600, pz=600, s.alphax=10, s.gammax=10, s.alphaz=300, 
                           s.gammaz=10, varsigma=0, IV.strength=1, a0=0)
pval.null <- HDOverTest_IV(D=out.data.null$D, X=out.data.null$X, Z=out.data.null$Z, Y=out.data.null$Y)
cat("(n=200, px=600, pz=600, a0=0, gammaz1, s.alphaz=0.5pz, varsigma=0) for H0+dense,", "p-value is", pval.null)

# H1 + dense
out.data.alter <- gene_data(n=200, px=600, pz=600, s.alphax=10, s.gammax=10, s.alphaz=300, 
                           s.gammaz=10, varsigma=0.5, IV.strength=1, a0=0)
pval.alter <- HDOverTest_IV(D=out.data.alter$D, X=out.data.alter$X, Z=out.data.alter$Z, Y=out.data.alter$Y)
cat("(n=200, px=600, pz=600, a0=0, gammaz1, s.alphaz=0.5pz, varsigma=0.5) for H1+dense,", "p-value is", pval.alter)



