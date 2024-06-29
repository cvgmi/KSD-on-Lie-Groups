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
library(rotations)

set.seed(123)

F <- matrix( c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#F <- matrix( c(8.5, 1.1, 4.1, 7.8, 3.9, 6, 4.3, 6.4, 4.8), ncol = 3)   ### An arbitrary F matrix

F_svd <- svd(F)

tau <- 1
n <- 3
I <- diag(n)

n_samples <- 100


# X <- rmatrixfisher(n_samples, F)
# X_mean <-rowMeans(X, dims = 2)
# X_mean_svd <- svd(X_mean)

cayley_kappa <- 2.0

X <- ruars(n = n_samples, rangle = rcayley, kappa = cayley_kappa, space = 'SO3')
X <- t(X)
X <- array(X, dim=c(3, 3, n_samples))



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
d_ksd <- FDist2(F, F_ksd)
  
  
K0 <- 0
K_ksd <- 0
gram_K_ksd <- matrix(0,n_samples,n_samples)
for (i in 1:n_samples){
  for (j in 1:n_samples){
    Xi <- X[1:3,1:3,i]
    Xj <- X[1:3,1:3,j]
    c <- (n - 1) * tau * tra[i,j] * ker[i,j] / 2
      
    XY <- t(Xi) %*% (F + tau * Xj)
    YX <- t(Xj) %*% (F + tau * Xi)
    trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
    K0 <- K0 + (c + trace * ker[i,j])
      
    XY <- t(Xi) %*% (F_ksd + tau * Xj)
    YX <- t(Xj) %*% (F_ksd + tau * Xi)
    trace <- sum(diag(t(skew_symmetrization(XY)) %*% skew_symmetrization(YX)))
    gram_K_ksd_ij <- c + trace * ker[i,j]
    gram_K_ksd[i,j] <- gram_K_ksd_ij
    K_ksd <- K_ksd + gram_K_ksd_ij
    }
}

K0 = K0 / (n_samples*n_samples)
K_ksd = K_ksd / (n_samples*n_samples)

end = Sys.time()
end-start

# gram_K_ksd_svd <- svd(gram_K_ksd)
lambda <- eigen(gram_K_ksd)$values  / (n_samples)

n_bootstrap <- 1000
beta <- 0.1
Z <- matrix(rnorm(n_samples*n_bootstrap, sd=1), nrow = n_samples, ncol = n_bootstrap)

W <- lambda %*% (Z**2)

gamma <- quantile(W, probs = 1-beta)
statistics <- n_samples * K_ksd 
gamma
statistics

quantiles <- c(0*100)

for (i in 1:100){
  quantiles[i] <- quantile(W, probs = i*0.01)
}

results <-c(n_samples,K_ksd,cayley_kappa,statistics,quantiles)
results

file_name <- sprintf("%f_cayley.csv",cayley_kappa)

# write.csv(results,file=file_name)

