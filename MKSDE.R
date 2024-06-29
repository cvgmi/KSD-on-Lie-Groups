library(SMFilter)
library(fastmatrix)

skew_symmetrization <- function(X) {
  # X: n x n 
  return(0.5 * (X - t(X)))
}

symmetrization <- function(X) {
  # X: n x n 
  return(0.5 * (X + t(X)))
}

kernel <- function(X, Y, tau = 0.1) {
  # matrix kernel
  # X: n x n SO(n)
  # Y: n x n SO(n)
  XTY <- t(X) %*% Y
  trace_XTY <- sum(diag(XTY))
  return(exp(tau * trace_XTY))
}

calculate_Bij <- function(X,Y,tau = 0.1){
  XTY <- t(X) %*% Y
  k <- kernel(X,Y)
  return((X-Y) %*% skew_symmetrization(XTY) * tau * k)  
}

calculate_Aij <- function(X,Y,tau = 0.1, snn){
  n = dim(X)[1]
  I = diag(n)
  XYT <- (X) %*% t(Y)
  IKXYT <- kronecker.prod(I,XYT)
  XKY <- kronecker.prod(t(X),Y)
  k <- kernel(X,Y)
  return((IKXYT - XKY %*% snn) * k)  
}

von_mises_kernel <- function(F, X, Y, tau = 0.1) {
  # X: n x n SO(n)
  # Y: n x n SO(n)
  n <- dim(X)[1]
  k <- kernel(X, Y, tau)
  trace_XTY <- sum(diag(t(X) %*% Y))
  c <- (n - 1) * tau * trace_XTY * k
  XY <- t(X) %*% (F + tau * Y)
  YX <- t(Y) %*% (F + tau * X)
  trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
  return(c + trace * k)
}

# Define the vectors
vector1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
vector4 <- c(0, 0, 0, 1, 0, 0, 0, 0, 0)
vector7 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
vector2 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
vector5 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
vector8 <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
vector3 <- c(0, 0, 1, 0, 0, 0, 0, 0, 0)
vector6 <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
vector9 <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)

# Concatenate the vectors using rbind
S_NN <- rbind(vector1, vector4, vector7, vector2, vector5, vector8, vector3, vector6, vector9)

print(S_NN)



library(Directional)

set.seed(123)
Ns = seq(100, 200, by=50)
d_mle = rep(0, length(Ns))
d_ksd = rep(0, length(Ns))
K0s = rep(0, length(Ns))
K_ksds = rep(0, length(Ns))
K_mles = rep(0, length(Ns))


#F <- matrix( c(25, 0, 0, 0, 5, 0, 0, 0, 1), ncol = 3)

F <- matrix( c(8.5, 1.1, 4.1, 7.8, 3.9, 6, 4.3, 6.4, 4.8), ncol = 3)   ### An arbitrary F matrix

F_svd <- svd(F)

tau <- 1
n <- 3
I <- diag(n)
for (kk in 1:length(Ns)) {
  n_samples <- Ns[kk]
  X <- rmatrixfisher(n_samples, F)
  X_mean <-rowMeans(X, dims = 2)
  X_mean_svd <- svd(X_mean)
  

  F_mle_D = X_mean_svd$d*3 #approximation
  

  F_mle <- X_mean_svd$u %*% diag(F_mle_D) %*% t(X_mean_svd$v)
  d_mle[kk] <- FDist2(F, F_mle)
  
  tra <- matrix(0, n_samples, n_samples)
  ker <- matrix(0, n_samples, n_samples)
  A <- 0
  B <- 0
  start = Sys.time()
  for (i in 1:n_samples){
    for (j in 1:n_samples){
      Xi <- X[1:3,1:3,i]
      Xj <- X[1:3,1:3,j]
      XiTXj <- t(Xi) %*% Xj
      XiXjT <- (Xi) %*% t(Xj)
      # trace_XiTXj <- sum(diag(XiTXj))
      # k <- exp(tau * trace_XiTXj)
      tra[i,j] <- sum(diag(XiTXj))
      ker[i,j] <- exp(tau * tra[i,j])
      Bij <- ((Xi-Xj) %*% skew_symmetrization(XiTXj) * tau * ker[i,j]) 
      IKXiXjT <- kronecker.prod(I,XiXjT)
      XiKXj <- kronecker.prod(t(Xi),Xj)
      Aij <- ((IKXiXjT - XiKXj %*% S_NN) * ker[i,j]) 
      
      B <- B + Bij
      A <- A + Aij
    }
  }
  B <- B / (n_samples*n_samples*2)
  A <- A / (n_samples*n_samples*2)
  b <- fastmatrix::vec(B)
  AI  <- solve(A)
  F_ksd_vector <- AI %*% b
  F_ksd <- matrix(F_ksd_vector, ncol = 3)
  d_ksd[kk] <- FDist2(F, F_ksd)
  
  
  K0 <- 0
  K_ksd <- 0
  K_mle <- 0
  for (i in 1:n_samples){
    for (j in 1:n_samples){
      Xi <- X[1:3,1:3,i]
      Xj <- X[1:3,1:3,j]
      c <- (n - 1) * tau * tra[i,j] * ker[i,j]
      
      
      XY <- t(Xi) %*% (F + tau * Xj)
      YX <- t(Xj) %*% (F + tau * Xi)
      trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
      K0 <- K0 + (c + trace * ker[i,j])
      
      XY <- t(Xi) %*% (F_mle + tau * Xj)
      YX <- t(Xj) %*% (F_mle + tau * Xi)
      trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
      K_mle <- K_mle + (c + trace * ker[i,j])
      
      XY <- t(Xi) %*% (F_ksd + tau * Xj)
      YX <- t(Xj) %*% (F_ksd + tau * Xi)
      trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
      K_ksd <- K_ksd + (c + trace * ker[i,j])
    }
  }
  K0s[kk] = K0 / (n_samples*n_samples)
  K_ksds[kk] = K_ksd / (n_samples*n_samples)
  K_mles[kk] = K_mle / (n_samples*n_samples)
  end = Sys.time()
  end-start
}

plot(Ns, d_mle, type='l', col='red', ylim = c(min(d_mle, d_ksd), max(d_mle, d_ksd)))
lines(Ns, d_ksd, col='green')

plot(Ns, K0s, type='l', ylim = c(min(K0s, K_ksds, K_mles), max(K0s, K_ksds, K_mles)))
lines(Ns, K_ksds, col='green')
lines(Ns, K_mles, col='red')

df <- cbind(Ns,K0s,K_ksds,K_mles,d_ksd,d_mle)
#write.csv(df,file='2551I.csv')
