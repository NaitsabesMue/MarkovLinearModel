################################
#
# This file contains an implementation of the estimators introduced in
# "Estimations of means and variances in a Markov linear model" 
# by A. Gutierrez and S. Mueller
# arxiv:
#
# We follow the notations of the above paper, where more information and details
# can be found.
#
# Note that each observation consists of a vector=path and a value b.
## We use tibble as a data structure.
################################


#### Packages

library("ggplot2")
library("dplyr")
library("ggpubr")
library("tibble")
library("list")


############
# The function sample_path samples a sample path according to the transition 
# matrices Q
##########

sample_path <- function(number_col, number_lines, Q)
{
  result <- numeric(number_col)
  result[1] <-
    sample(1:number_lines[1], 1, replace = FALSE, prob = Q[[1]])
  for (i in 2:number_col)
    result[i] <-
    sample(1:number_lines[i], 1, replace = FALSE, prob = Q[[i]][result[i - 1],])
  return(result)
}


############
# The function sample_path_quality samples a sample path according to the transition 
# matrices Q together with a realization of quality b; the matrix S is Gaussian with mean 
# ES and variance V 
##########

sample_path_quality <- function(number_col, number_lines, Q, ES, V)
{
  path <- numeric(number_col)
  path[1] <-
    sample(1:number_lines[1], 1, replace = FALSE, prob = Q[[1]])
  for (i in 2:number_col)
    path[i] <-
    sample(1:number_lines[i], 1, replace = FALSE, prob = Q[[i]][path[i - 1], ])
  b <- 0
  for (i in 1:number_col)
    b <- b + rnorm(1, ES[path[i], i], sqrt(V[path[i], i]))
  df <-  as.data.frame(matrix(c(path, b), nrow = 1))
  names(df)[number_col + 1] <- "b"
  return(df)
}

############
# The function T_U calculates the estimator T_U^(n)
##########

T_U <- function(observation, number_col, number_lines)
{
  T <- matrix(0, nrow = max(number_lines), ncol = number_col)
  D <- matrix(0, nrow = max(number_lines), ncol = number_col)
  n <- nrow(observation)
  for (k in 1:n) {
    path <- unlist(observation[k, 1:number_col])
    b <- observation$b[k]
    for (j in 1:number_col) {
      D[path[j], j] <- D[path[j], j] + 1
      T[path[j], j] <- T[path[j], j] + b
    }
  }
  return(T / D)
}

############
# The function Q_estimator calculates an estimator for the transition kernels Q
##########

Q_estimator <- function(observation, number_col, number_lines) {
  n <- nrow(observation)
  Q <- list()
  estimate <- numeric(number_lines[1])
  for (k in 1:n) {
    estimate[as.numeric(observation[k, 1])] <-
      estimate[as.numeric(observation[k, 1])] + 1
  }
  Q[[1]] <- estimate / n
  for (j in 2:(number_col)) {
    estimate <- matrix(0, nrow = number_lines[j - 1], ncol = number_lines[j])
    for (k in 1:n) {
      index1 <- as.numeric(observation[k, j - 1])
      index2 <- as.numeric(observation[k, j])
      estimate[index1, index2] <- estimate[index1, index2] + 1
    }
    Q[[j]] <- estimate / rowSums(estimate)
  }
  return(Q)
}

############
# The function P_path calculates the probability of a given path according to Q
##########

P_path <- function(path, Q){
  n <- length(path)
  result <- Q[[1]][path[1]]
  for (i in 2:n){
    result <- result * Q[[i]][path[i-1], path[i]]
  }
  return(result)
}


############
# The function T_QU calculates the unbiased estimator \widehat{T}^{(n)}_{Q,U}
##########

T_QU <- function(observation, number_col, number_lines, Q)
{
  T <- matrix(0, nrow = max(number_lines), ncol = number_col)
  D <- matrix(0, nrow = max(number_lines), ncol = number_col)
  n <- nrow(observation)
  for (k in 1:n) {
    path <- unlist(observation[k, 1:number_col])
    b <- observation$b[k]
    for (j in 1:number_col) {
      const <- prod(number_lines[-j])
      D[path[j], j] <- D[path[j], j] + 1
      T[path[j], j] <- T[path[j], j] + b / P_path(path, Q)  / const
    }
  }
  return(T / D * D / n)
}

############
# The function S_U calculates the estimator \widehat{\Sigma}^{(n)}
##########

S_U <- function(observation, number_col, number_lines)
{
  T <- matrix(0, nrow=max(number_lines), ncol=number_col)
  S <- matrix(0, nrow=max(number_lines), ncol=number_col)
  D <- matrix(0, nrow=max(number_lines), ncol=number_col)
  n <- nrow(observation)
  for (k in 1:n){
    path <- unlist(observation[k, 1: number_col])
    b <- observation$b[k]
    for (j in 1:number_col){
      D[path[j],j] <- D[path[j],j]+1
      T[path[j],j] <- T[path[j],j]+b
      S[path[j],j] <- S[path[j],j]+b^2
    }
  }
  return(S/D - (T/D)^2)    
} 

############
# The function S_QU calculates the unbiased estimator \widehat{\Sigma}^{(n)}_{Q,U}
##########

S_QU <- function(observation, number_col, number_lines, Q)
{
  T <- matrix(0, nrow = max(number_lines), ncol = number_col)
  S <- matrix(0, nrow = max(number_lines), ncol = number_col)
  n <- nrow(observation)
  for (k in 1:n) {
    path <- unlist(observation[k, 1:number_col])
    b <- observation$b[k]
    for (j in 1:number_col) {
      const <- prod(number_lines[-j])
      T[path[j], j] <- T[path[j], j] + b / P_path(path, Q)  / const
      S[path[j], j] <- S[path[j], j] + b ^ 2 / P_path(path, Q)  / const
    }
  }
  return(S / n - (T / n) ^ 2)
}


