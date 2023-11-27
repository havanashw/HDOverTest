library(MASS)
library(glmnet)

HDOverTest <- function(D, X, Z, Y) {
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

