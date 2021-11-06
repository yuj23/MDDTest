####################################################
MDD test
####################################################

#set work directory
setwd("")

#load packages
library("Rcpp")
Rcpp::sourceCpp("MV_Statistic.cpp")
library(movMF)

##  MV statistics 
MV <- function(X, groups){
  n <- nrow(X)
  R <- length(unique(groups))
  group <- as.integer(factor(groups))
  distance <- as.matrix(dist(X))
  result <- MV_Statistic(distance, n, group, R)
  return(result)
}

## permutation test to get p-value
MV_test <- function(X, groups, k=200){
  n <- nrow(X)
  R <- length(unique(groups))
  group <- as.integer(factor(groups))
  distance <- as.matrix(dist(X))
  result <- sapply(2:k, function(x){perm <- sample(1:n); MV_Statistic(distance[perm, perm], n, group, R)})
  result <- c(result, MV_Statistic(distance, n, group, R))
  p_value <- sum(result>=result[k])/k
  return(p_value)
}

# independence test  2-sample
simu = 200
inde_test = function(i,R = 2,sample_size=120){
       num1 <- as.integer(sample_size*0.5)
       num2 <- sample_size-num1
       #X <- rmovMF(sample_size,c(0,0,0))
       X <- matrix(rnorm(sample_size*3,mean=0),nrow=sample_size)
       groups <- c(sample(c(1,1),num1,replace=T),sample(c(2,2),num2,replace=T))
       perm <- sample(1:sample_size)
       X <- X[perm,]
       groups <- groups[perm]
       p_value = MV_test(X, groups)
       return(p_value)}
p_values_inde = sapply(1:simu,inde_test)
type_I_inde = round(sum(p_values_inde<=0.05)/simu,3)
type_I_inde

## dependence test
simu = 200
depen_test = function(i,R = 2,sample_size=120){
  num1 = 0.3*sample_size
  num2 = 0.7*sample_size
  X1 <-rmovMF(num1,c(0,0,0))
  X2 <-rmovMF(num2,c(1,1,1))
  X <- rbind(X1,X2)
  groups <- c(sample(c(1,1),num1,replace=T),sample(c(2,2),num2,replace=T))
  perm <- sample(1:sample_size)
  X <- X[perm,]
  groups <- groups[perm]
  p_value = MV_test(X, groups)
  return(p_value)}
p_values_de = sapply(1:simu,depen_test)
type_I_de = round(sum(p_values_de<=0.05)/simu,4)
type_I_de

##set datasets
# test 1
X1 <-matrix(rnorm(300,mean=1),nrow=150)
X2 <-matrix(rnorm(100,mean=5),nrow=50)
X_ori <- matrix(rnorm(600,mean=0),nrow=200)

# test 2
X1 <-rmovMF(50, c(0,0,0))
X2 <-rmovMF(150,c(1,1,1))
X_ori <- rmovMF(200,c(0,0,0))
