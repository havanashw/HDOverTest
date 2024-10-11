# HDOverTest
R code for "Testing overidentifying restrictions on high-dimensional instruments and covariates".

## Usage

```{R}
HDOverTest(D, X, Z, Y)
```

## Required Packages
- `glmnet`

## Inputs
- `D`: Endogenous treatment variable with length n.
- `X`: Exogenous covariates, a matrix of n*px, where px is the dimension of `X`.
- `Z`: Instrumental variables, a matrix of n*pz, where pz is the dimension of `Z`.
- `Y`: Response variable with length n.

## Examples
```{R}
library(MASS)
library(glmnet)

source("HDOverTest.R")

## generate data (See the Simulations section of the paper for detailed descriptions)
gene_data <- function(n, px, pz, s.alphax, s.gammax, s.alphaz,
                      s.gammaz, varsigma, IV.strength, a0) {
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

  # treatment effect
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

## implement
set.seed(2023)
# H0 + dense
out.data.null <- gene_data(n=200, px=600, pz=600, s.alphax=10, s.gammax=10, s.alphaz=300, 
                           s.gammaz=10, varsigma=0, IV.strength=1, a0=0)
pval.null <- HDOverTest(D=out.data.null$D, X=out.data.null$X, Z=out.data.null$Z, Y=out.data.null$Y)
cat("(n=200, px=600, pz=600, a0=0, gammaz1, s.alphaz=0.5pz, varsigma=0) for H0+dense,", "p-value is", pval.null)

# H1 + dense
out.data.alter <- gene_data(n=200, px=600, pz=600, s.alphax=10, s.gammax=10, s.alphaz=300, 
                           s.gammaz=10, varsigma=0.5, IV.strength=1, a0=0)
pval.alter <- HDOverTest(D=out.data.alter$D, X=out.data.alter$X, Z=out.data.alter$Z, Y=out.data.alter$Y)
cat("(n=200, px=600, pz=600, a0=0, gammaz1, s.alphaz=0.5pz, varsigma=0.5) for H1+dense,", "p-value is", pval.alter)
```
