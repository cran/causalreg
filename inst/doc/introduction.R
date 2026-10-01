## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----poisson-glm--------------------------------------------------------------
library(causalreg)

n <- 1000
set.seed(123)
X1 <- rnorm(n)
Y  <- rpois(n, exp(X1))
X2 <- log(Y + 1) + rnorm(n, 0, 0.3)
data <- data.frame(X1, X2, Y)

## ----poisson-all--------------------------------------------------------------
result <- cglm(Y ~ X1 + X2, "poisson", data, pval = "chi-square", search = "all")
result$model.opt

## ----poisson-all-details------------------------------------------------------
# All models considered
unlist(result$models)

# Their p-values (acceptance means no evidence to reject Pearson risk = 1)
result$pv

# Their BIC values
result$bic

## ----poisson-step-------------------------------------------------------------
result_step <- cglm(Y ~ X1 + X2, "poisson", data, pval = "chi-square", search = "stepwise")
result_step$model.opt

# Models visited during the search
unlist(result_step$models)

## ----binomial-glm-------------------------------------------------------------
n <- 2000
set.seed(123)
X1 <- rnorm(n)
Y  <- rbinom(n, 1, exp(X1) / (1 + exp(X1)))
flip <- rbinom(n, 1, 0.1)
X2 <- (1 - flip) * Y + rnorm(n, 0, 0.3)
data <- data.frame(X1, X2, Y)

set.seed(1)
result <- cglm(Y ~ X1 + X2, "binomial", data, pval = "bootstrap", search = "all")
result$model.opt

## ----poisson-gam--------------------------------------------------------------
n <- 1000
set.seed(123)
X1 <- rnorm(n)
Y  <- rpois(n, exp(sin(X1)))
X2 <- log(Y + 1) + rnorm(n, 0, 0.5)
data <- data.frame(X1, X2, Y)

result <- cgam(Y ~ s(X1) + s(X2), "poisson", data, pval = "chi-square", search = "all")
result$model.opt

## ----five-cov, eval = FALSE---------------------------------------------------
# set.seed(12)
# n <- 3000
# X1 <- rnorm(n)
# X2 <- rnorm(n, X1, 0.5)
# X3 <- rnorm(n, 0, 1)
# X4 <- rnorm(n, X2, 0.5)
# Y  <- rbinom(n, 1, exp(0.8 * X2 - 0.9 * X3) / (1 + exp(0.8 * X2 - 0.9 * X3)))
# flip <- rbinom(n, 1, 0.1)
# X5 <- (1 - flip) * Y + flip * (1 - Y) + rnorm(n, 0, 0.3)
# dat <- data.frame(X1, X2, X3, X4, X5, Y)
# 
# # Exhaustive search (evaluates all 2^5 - 1 = 31 subsets)
# set.seed(1)
# mod_all <- cglm(Y ~ X1 + X2 + X3 + X4 + X5, "binomial", dat,
#                 pval = "bootstrap", search = "all")
# mod_all$model.opt
# #> [1] "Y ~ X2 + X3"
# 
# # Stepwise search (much faster)
# set.seed(1)
# mod_step <- cglm(Y ~ X1 + X2 + X3 + X4 + X5, "binomial", dat,
#                  pval = "bootstrap", search = "stepwise")
# mod_step$model.opt
# #> [1] "Y ~ X2 + X3"

