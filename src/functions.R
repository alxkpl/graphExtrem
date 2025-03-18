#---------------------- Functions for travariate document ----------------------

#' 3x3 variogram generator for a trivariate Husler-Reiss model
#'
#' @param k A positive number to set the range of the coefficients.
#' @return The coefficients of a valid Husler-Reiss 3x3 variogram where the corresponding
#' graphical model doesn't have edge between nodes 1 and 2.
#' The vector correspond to : 
#'      - a : Gamma_12
#'      - b : Gamma_13
#'      - c : Gamma_23
#' @examples
#' random_Gamma12(5)
random_Gamma12 <- function(k){
  b <- k * runif(1)
  c <- k * runif(1)
  a <- b + c 
  Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
  
  while(!(checkGamma(Gamma, returnBoolean = TRUE, alert = FALSE))){
    b <- k * runif(1)
    c <- k * runif(1)
    a <- b + c
    Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
  }
  return(matrix(c(a, b, c), nc = 3))
}

#' Transform a set of variogram's parameter to a matrix.
#'
#' @param Gamma_params A numerical vector: parameters of the variogram
#'
#' @return The corresponding symetric matrix.
#'
#' @examples
to_matrix <- function(Gamma_params){
  # Parameters of the variogram
  a <- Gamma_params[1]
  b <- Gamma_params[2]
  c <- Gamma_params[3]
  
 return(
   matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
   )
}

#' 2-variate extremal coefficient of a bivariate Husler-Reiss distribution.
#'
#' @param gamma A real number: the value of the Husler-Reiss parameter.
#'
#' @return The value of the 2-variate extremal coefficient for a Husler-Reiss distribution
#' with parameter gamma. This value is \Lambda(1, 1) = 2 phi(sqrt(gamma)/2), where 
#' phi is the distribution function a standard gaussian.
#'
#' @examples
#' theta(1)
theta <- function(gamma){
  return(
    2 * pnorm(sqrt(gamma) / 2)
    )
  }


#' 3-variate extremal coefficient of a trivariate Husler-Reiss distribution.
#'
#' @param Gamma_parms A vector of size 3 symbolizing a conditionnal negtaive matrix :
#' the Husler-Reiss parameters.
#'
#' @returnThe value of the 3-variate extremal coefficient for a Husler-Reiss distribution
#' with variogram represented by the vector Gamma_params.
#'
#' @examples
#' Gamma_params <- c(5, 1, 4)
#' lambda_2(Gamma_params)
lambda_2 <- function(Gamma_params){
  # Parameters of the variogram
  a <- Gamma_params[1]
  b <- Gamma_params[2]
  c <- Gamma_params[3]
  
  # Computation of the correlation for the Gaussian distribution function
  rho_1 <- (a + b - c) / (2 * sqrt(a * b))
  rho_2 <- (a + c - b) / (2 * sqrt(a * c))
  rho_3 <- (b + c - a) / (2 * sqrt(b * c))
  
  
  if(sum(round(abs(c(rho_1, rho_2, rho_3)), 8) >= 1)){
    return(NA)
  }else{
  return(
    pnorm2d(x = sqrt(a) / 2, y = sqrt(b) / 2, rho = rho_1)[1] +
      pnorm2d(x = sqrt(a) / 2, y = sqrt(c) / 2, rho = rho_2)[1] +
      pnorm2d(x = sqrt(b) / 2, y = sqrt(c) / 2, rho = rho_3)[1] 
    )
  }
}


#' Trivariate extremal coefficient in a Husler-Reiss model.
#'
#' @param Gamma_params The coefficient of the HR variogram of a size-3 (sub-)vector . 
#'
#' @return The value of the trivariate coefficient in the corresponding HR 
#' random vector.
#'
#' @examples
#' chi_trivariate_HR(random_Gamma12())
chi_trivariate_HR <- function(Gamma_params){
  return(
    3 - theta(Gamma_params[1]) - theta(Gamma_params[2]) - theta(Gamma_params[3]) + lambda_2(Gamma_params)
    )
}


#' Title
#'
#' @param Gamma_params 
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
cond_trivariate_HR <- function(Gamma_params, i){
    return(
      chi_trivariate_HR(Gamma_params) / (2 - theta(Gamma_params[i]))
      )
}


#' Give a matrix of null variogram, parameterised by the diagonal values
#'
#' @param diagonal Numeric vector : diagonal's value of a matrix
#'
#' @returns A matrix which belong to the kernel of the linear application gamma, 
#' computing the variogram for a given matrix. It is characterized by the relation :
#'                          a_ij = (a_ii + a_jj) / 2
#' @examples
#' ker_gamma(1:4)
ker_gamma <- function(diagonal){
  
  n <- length(diagonal)                 # dimension of the matrix
  
  one <- rep(1, n)                      # vector with only ones
  
  return(0.5 * (one%*%t(diagonal) + diagonal%*%t(one)))
}


#' Semi-definite checker
#'
#' @param matrix A symmetric matrix.
#'
#' @returns TRUE if the matrix is semi-definite positive, FALSE otherwise. The 
#' verification is done by checking the sign of the eigen-values of the matrix.
#'
#' @examples
#' 
semi_def <- function(matrix){
  return(
    !sum(eigen(matrix)$values < -1e-10)        # for numerical errors
    )
}

#----------------------- Gradient descent for clustering -----------------------
##========= Functions for the building of some particular matrices =============

#' Computation of the matrix of clusters
#'
#' @param clusters a list of vector : each vector gives the element of a cluster.
#'
#' @returns The d x K matrix U such that u_jk = 1 if the variable j belongs to 
#' cluster C_k and 0 otherwise.
#'
#' @examples
#' clusters <- list(c(1,2,3),c(4,5))
#' U_matrix(clusters)
U_matrix <- function(clusters){
  K <- length(clusters)
  d <- length(unlist(clusters))
  U <- matrix(rep(0, d * K), nc = K)
  for(k in 1:K){
    for(j in 1:d){
      if(j %in% clusters[[k]]){
        U[j, k] <- 1
      }
    }
  }
  return(U)
}

#' From R matrix compute theta matrix
#'
#' @param R a K x K matrix : the matrix of the clusters coefficients
#' @param clusters a list of vector : each vector gives the element of a cluster.
#'
#' @returns Return the theta matrix from the value of the R matrix of clusters
#' coefficient using :
#' 
#'                            theta = U R U^t + A 
#'                            
#' where U is the cluster matrix and a diagonal matrix such that the rows of theta sum 
#' to zero : 
#'                            a_ii = - sum_l p_l r_kl
#' for i in cluster C_k.
#'
#' @examples
#' 
#' R <- matrix(c(1, 2,
#'               2, 5), nr = 2)
#' clusters <- list(c(1,2,3),c(4,5))
#' build_theta(R, clusters)
#' 
build_theta <- function(R, clusters){
  U <- U_matrix(clusters)
  URUt <- U %*% R %*% t(U)
  a <-  - rowSums(URUt)
  return(
    URUt + diag(a)
  )
}


#' Compute the trace cluster matrix vector.
#'
#' @param gamma a dxd matrix : an estimation of the variogram gamma.
#' @param clusters a list of vector : each vector gives the element of a cluster.
#'
#' @returns Returns a vector of size K of entries tr(Gamma_C_l).
#' @export
#'
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' gamma <- matrix(1:25, nc = 5)
#' trace_vector(gamma, clusters)
trace_vector <- function(gamma, clusters){
  K <- length(clusters)
  T.vector <- rep(0, K)
  i <- 0
  for(index in clusters){
    i <- i + 1
    T.vector[i] <-  sum(diag(gamma[index, index]))
  }
  return(T.vector)
}

#' Variogram transformation application gamma
#'
#' @param sigma A d x d numeric matrix.
#'
#' @returns For a symmetric positive matrix sigma (covariance matrix), return the 
#' corresponding variogram matrix. Can be used for other but with no interpretation.
#'
#' @examples
#' s_sigma <- matrix(rnorm(16, 2), nc = 4)
#' gamma_function(s_sigma %*% t(s_sigma))
gamma_function <- function(sigma){
  indic <- rep(1, nrow(sigma))
  return(
    tcrossprod(diag(sigma), indic) + tcrossprod(indic, diag(sigma)) - 2 * sigma
  )
}


#' Crout factorization algorithm.
#'
#' @param A a d x d symmetric positive matrix.
#' @param tol a positive value : tolerance for the zero diagonal
#'
#' @returns The LU Crout decomposition of a matrix A. For A a symmetric positive 
#' matrix, there exists a LU decomposition such that : 
#'                                   A = L' L 
#'
#' @examples
#'A <- matrix(c(1,2,3,
#'              2,5,6,
#'              3,6,9), nc = 3)
#' L <- crout_factorisation(A)
#' L %*% t(L)
crout_factorisation <- function(A, tol = 1e-12){
  A <- unname(as.matrix(A))
  n <- nrow(A)
  if(n != ncol(A)){
    stop("no square matrix.")
  }
  if(!semi_def(A)){
    stop("no positive semi definite matrix.")
  }
  d <- rep(0, n)
  L <- diag(rep(1, n))
  d[1] <- A[1, 1]
  for(i in 2:n){
    for(j in 1:(i - 1)){
      L[i, j] <- (1 / d[j]) * (A[i, j] - sum(L[i, ] * L[j, ] * d))
    }
    d[i] <- A[i, i] - sum(d*(L[i, ])**2)
  }
  
  return(L %*% diag(sqrt(d * (d > tol))))
}

#' Moore-Penrose pseudo inverse
#'
#' @param A a d x d symmetric positive semi-definite matrix.
#'
#' @returns Computes the Moore-Penrose inverse of a matrix. The calculation is 
#' done thanks to an article and if  : 
#'                                A = L L^t 
#' (with L having no zero-columns) then we have : 
#'                        A^+ = L (L^t L)^-1 (L^t L)^-1 L^t
#'
#' @examples
#'A <- matrix(c(1,2,3,
#'              2,5,6,
#'              3,6,9), nc = 3)
#' psolve(A)
psolve <- function(A, tol = 1e-12){

  S <- crout_factorisation(A, tol = tol)
  
  L <- S[, which(diag(S) != 0)]         # to get no null columns
  
  return(
    L %*% solve(t(L) %*% L) %*% solve(t(L) %*% L) %*% t(L)
  )
}

#' Computation of the clustered weight matrix
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters : the weight matrix for each pair of clusters :
#'                          W_kl = sum_(i in C_k) sum_(j in C_l) w_ij
#'
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' W <- matrix(1:25, nc = 5)
#' W_c <- weight_clustered(W)
#' W_c(clusters)
weight_clustered <- function(weights){
  function(clusters){
    U <- U_matrix(clusters) 
    return(
      t(U) %*% weights %*% U
    )
  }
}


#' Generalised determinant
#'
#' @param A a d x d real valued matrix.
#' @param tol a positive value.
#'
#' @returns Compute the generalised determinant of a matrix A. We recall that the 
#' generalised determinant is an extension of the determinant for singular matrix. 
#' It corresponds to the product of all the non zero eigen values.
#'
#' @examples
#' A <- matrix(c(1,2,3,
#'               2,5,6,
#'               3,6,9), nc = 3)
#' gen_det(A)
gen_det <- function(A, tol = 1e-10){
  res <- eigen(A, only.values = TRUE)$values
  return(
    prod(res[res > tol])
  )
}
##======= Computation of the gradient of the initial log-likelihood ============

#' Negative log-likelihood computation.
#'
#' @param gamma A d x d matrix: the empirical variogram matrix.
#'
#' @returns For a fixed variogram gamma, compute for a set of clusters and 
#' corresponding R matrix, the value of the associated negative likelihood defined
#' by : 
#'                  nllh = - log(|theta|_+) - 1/2 trace(gamma * theta)
#' where |.|_+ is the generalised determinant.
#' 
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' nllh <- neg_likelihood(gamma) 
#' nllh(R, clusters)
neg_likelihood <- function(gamma){
  function(R, clusters){
    # Building of the theta matrix from R
    theta <- build_theta(R, clusters)
    
    # Computation of the log-determinant part
    log_det <- log(gen_det(theta))
    
    # Computation of the trace part
    tr <- sum(diag(gamma %*% theta))
    
    return(
      - log_det - .5 * tr
    )
  }
}


#' Gradient of the negative log likelihood without penalty
#'
#' @param gamma a d x d variogram matrix.
#'
#' @returns A function of the R matrix and clusters and compute the gradient 
#' matrix of the negative log likelihood for a fixed variogram gamma. The gradient 
#' matrix can be computed by : 
#' 
#'                        dnllh = dlog + dtrace
#' where :
#'      - dlog = t(U) g(Theta_p) U - 0.5 diag(t(U) g(Theta_p) U)  
#'      - dtrace = - t(U) gamma U  + 0.5 diag(t(U) gamma U)  
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' gradient <- nloglike_grad_np(gamma)
#' gradient(R, clusters)
nloglike_grad_np <- function(gamma){
  function(R, clusters){
    # Get tu U matrix of clusters indicators
    U <- U_matrix(clusters) 
    p <- colSums(U)
    
    # Computation of gamma(Theta_p) with Theta_p the Penrose inverse of Theta
    G_theta_p <- build_theta(R, clusters) |> 
                    psolve() |> 
                    gamma_function()
    
    # Gradient (factorisation of the formula)
    dlog <- t(U) %*% (G_theta_p - gamma) %*% U 
    diag <- - .5 * diag(diag(dlog * (p == 1)))
    
    return(dlog + diag)
    
  }
}


##================= Computation of the gradient of the penalty =================

#' Cluster distance squared function
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of a cluster.
#'
#' @returns A function of cluster number : compute the square distance between two 
#' cluster for the distance defined in section 4.2 in cluster document.
#'
#' @examples
#' 
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' D2 <- D_tilde2_r(R, clusters)
#' D2(1, 2)
D_tilde2_r <- function(R, clusters){
  function(k, l){
    if(k == l){
      return(0)         # 0 distance for the same cluster even with the symmetry of the matrix
    }else{
      # Parameters of clusters
      K <- length(clusters)               # Number of clusters
      p <- sapply(clusters, length)       # Vector of cluster's size
      return(
        # for fixed k, l the square difference is multiplied by 1-p_k and 1-p_l
        sum((p - ((1:K) %in% c(k, l))) * (R[k, ] - R[l, ])**2)
        )
    }
  }
}

#' Penalty function.
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns For fixed weights, returns a function which compute the value of the 
#' penalty for chosen clusters and corresponding R matrix. We recall the penalty 
#' is given by : 
#'                    P(R) = sum_{k<l} W_kl D^2(r_.k, r_.l)
#' with W_kl the clustered weight for clusters k and l. 
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4)) 
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' P <- penalty(W)
#' P(R, clusters)
penalty <- function(weights){
  # Fixing the weights for computing weights clustered
  get_W <- weight_clustered(weights)
  function(R, clusters){
    # Initialization
    K <- length(clusters)                   # Number of clusters
    D2 <- D_tilde2_r(R, clusters)           # Function for distance between clusters
    W <- get_W(clusters)                    # Weights clustered
    D <- matrix(rep(0, K * K), nc = K)      # Distance matrix for clusters
    
    # Computation of the distance matrix
    for(l in 2:K){
      for(k in 1:(l-1)){        # we keep only the lower triangular part
        D[k, l] <- D2(k, l) 
      }
    }
    return(
      sum(D * W)
    )
  }
}

#' Gradient matrix of distance between two columns
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of a cluster.
#'
#' @returns Return a functon of indices (k', l') computing the gradient matrix of 
#' tilde D^2(r_k', r_l'). See section 4.3.3 in cluster document for details.
#'
#' @examples
#' R <- matrix(c(-1,0,-2,
#'               0,-3,-1,
#'               -2,-1,-1), 3)
#' clusters <- list(c(1,3), c(2), 4)
#' grad <- gradient_D2(R, clusters)
#' grad(1, 2)
gradient_D2 <- function(R, clusters){
  # Initialization 
  K <- length(clusters)                   # Number of clusters
  p <- sapply(clusters, length)           # Vector of cluster's size
  function(k, l){
    A <- matrix(rep(0, K * K), nc = K)
    
    # Computation of the non zero row and column
    A[k, ] <- 2 * p * (R[k, ]- R[l, ])
    A[, l] <- 2 * p * (R[l, ]- R[k, ])
    A[k, l] <- 2 * ((p[k] - 1) * (R[k,l] - R[k, k]) + (p[l] - 1) * (R[k,l] - R[l, l]))
    
    # Build the symmetry of the gradient (except for the diagonal)
    grad <- A + t(A)
    
    # Computation of the diagonal
    grad[k, k] <- 2 * (p[k] - 1) * (R[k, k] - R[k, l])
    grad[l, l] <- 2 * (p[l] - 1) * (R[l, l] - R[k, l])
    
    return(
      grad
    )
    
  }
}


#' Computation of the penalty's gradient
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters and corresponding R matrix. Compute the gradient 
#' with fixed weight. The expression of the gradient is just the weighted sum of the 
#' gradient of each tilde D^2 where the weights are the clustered weights.
#' See equations in section 4.3.3 for details.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4)) 
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' dpen <- penalty_grad(W)
#' dpen(R, clusters)
penalty_grad <- function(weights){
  get_W <- weight_clustered(weights)
  function(R, clusters){
    # Initialization
    W <- get_W(clusters)                # Weight clustered
    K <- length(clusters)               # Number of clusters
    
    # Function of gradient of indices
    grad_D2 <- gradient_D2(R, clusters)
    
    res <- matrix(rep(0, K * K), nc = K)
    
    for(k in 1:(K - 1)){
      for(l in (k + 1):K){
        res <- res + W[k, l] * grad_D2(k, l)        # Weighted sum
      }
    }
    
    return(
      res
    )
    
  }
}

##================== Computation of the global neg-loglikelihood ===============

#' Computation of the penalised negative log-likelihood
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#'
#' @returns A function of clusters and the R matrix which compute the penalised 
#' negative log-likelihood of the model.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4)) 
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' f <- neg_likelihood_pen(gamma, W, 0.5)
#' f(R, clusters)
neg_likelihood_pen <- function(gamma, weights, lambda){
  nllh <- neg_likelihood(gamma)
  pen <- penalty(weights)
  function(R, clusters){
    return(
      nllh(R, clusters) + lambda * pen(R, clusters)
    )
  }
}



##=========================== Positive condition on R ==========================

sub_theta <- function(R, clusters){
  p <- sapply(clusters, length)           # Vector of cluster's size
  tilde_R <- - as.numeric(R %*% p)
  
  return(
    R %*% diag(p) + diag(tilde_R)
  )
}

##========================= Gradient descent algorithm =========================

#' Step for the gradient descent
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#' @param size_grid integer : size of the search grid for the optimal step.
#'
#' @returns A function of clusters and R matrix which returns the next step of
#' the optimisation for the gradient descent algorithm.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4)) 
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' f <- step_gradient(gamma, W, 0.5)
#' f(R, clusters)
step_gradient <- function(gamma, weights, lambda, size_grid = 100){
  # Initialization of functions
  dlog <- nloglike_grad_np(gamma)                     # Neg-lklh gradient part
  dpen <- penalty_grad(weights)                       # Penalty gradient part
  nllh <- neg_likelihood_pen(gamma, weights, lambda)  # Penalised negative log-likelihood function
  
  function(R, clusters){
    # Initialization
    p <- sapply(clusters, length)           # Vector of cluster's size
    
    # Gradient matrix computation
    grad <- dlog(R, clusters) + lambda * dpen(R, clusters)      
    
    # Grid line search for optimal gradient step
    # Grid line construction
    s_max <- min(
      abs((R %*% p) / (grad %*% p))         # Maximum step size to get positive matrix
    )
    if(s_max == 0){
      s <- seq(0, 1, 0.01)
    }else{
      s <- seq(0, s_max, s_max / size_grid)
    }
    
    # Searching
    s_opt <- 0
    score <- nllh(R, clusters)

    for(i in 2:(length(s) - 1)){
      # Better step size 
      if(score > nllh(R - s[i] * grad, clusters)){
        # Positive matrix checking
        if(semi_def(sub_theta(R - s[i] * grad, clusters))){ 
          s_opt <- s[i] 
          score <- nllh(R - s_opt * grad, clusters)
        }
      }
    }
    return(list(step = s_opt, gradient = grad))
  }
}

#' Function which merges clusters
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of a cluster.
#' @param eps positive value : minimal tolerance for merging clusters
#' @param cost a function : Cost function of the optimisation
#'
#' @returns Returns, if merging, a list of the new clusters and the corresponding R matrix, 
#' where the coefficient of the new clustered is computing by averaging the coefficient of 
#' the two previous clusters.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4)) 
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' cost <- neg_likelihood_pen(gamma, weights, 100000)
#' merge_clusters(R, clusters, cost = cost)
merge_clusters <- function(R, clusters, eps=1e-1, cost){
  # Initialization
  D <- D_tilde2_r(R, clusters)                    # Function of clusters distance
  K <- length(clusters)                           # Actual number of clusters
  
  # Computation of the distance matrix 
  distance <- matrix(rep(Inf, K * K), nc = K)
  
  for(k in 1:(K - 1)){
    for(l in (k + 1):K){
      distance[k, l] <- D(k, l)
    }
  }
  
  # Search of the two potential clusters to merge
  index <- as.numeric(which(distance == min(distance), arr.ind = T))
  k <- index[1] 
  l <- index[2]
  
  # Checking uselessness of merging
  if(distance[k, l] > eps){
    return(
      list(
        R = R, 
        clusters = clusters
      )
    )
  }
  
  p <- sapply(clusters, length)           # Vector of cluster's size
  
  # Case when merging give only one cluster
  if(nrow(R)==2){
    return(
      list(
        R = (p[1] * R[1, 1] + p[2] * R[1, 2]) / (p[1] + p[2]),
        clusters = c(clusters[[1]], clusters[[2]])
        )
    )
  }
  
  # Merging and computation of merged coefficient in the case where K>2
  new_clusters <- clusters[-l]
  
  new_clusters[[k]] <- c(clusters[[k]], clusters[[l]])      # New cluster
  
  # Coefficient calculation
  R_new <- R[-l, -l]
  
  R_new[k, -k] <- ((p[k] * R[k, ] + p[l] * R[l, ]) / (p[k] + p[l]))[-c(k, l)]
  R_new[-k, k] <- ((p[k] * R[k, ] + p[l] * R[l, ]) / (p[k] + p[l]))[-c(k, l)]
  R_new[k, k] <- R[k, l]
  
  # Final checking : decreasing of the negative log-likelihood
  if(cost(R, clusters) > cost(R_new, new_clusters)){
    return(
      list(
        R = R_new, 
        clusters = new_clusters
      )
    )
  }else{
    return(
      list(
        R = R, 
        clusters = clusters
      )
    )
  }
}

#' Gradient descent algorithm for Husler-Reiss graphical models clustering
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#' @param ... 
#'
#' @returns A function which compute the maximum likelihood estimator using 
#' cluster-path gradient descent and returns the estimation of the clusters and 
#' the corresponding R matrix.
#'
#' @examples
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- generate_random_Gamma(d = 4)
#' R <- matrix(c(1,0,0,-1,
#'               0,1,1,-2,
#'               0,1,1,-1,
#'               -1,-2,-1,1), nc = 4)                
#' Cluster_HR <- get_cluster(gamma, W, 100)
#' Cluster_HR(R)
get_cluster <- function(gamma, weights, lambda, ...){
  L <- neg_likelihood_pen(gamma, weights, lambda)
  step <- step_gradient(gamma, weights, lambda,...)
  function(R.init, it_max = 1000, eps_g = 1e-3){
    # Initialization
    d <- nrow(gamma)
    R <- R.init
    clusters <- as.list(1:d)
    gradstep <- list(gradient = eps_g + 1)
    cpt <- 1
    while((cpt < it_max)&(length(R) != 1)&(sum(gradstep$gradient**2)>eps_g)){
      # Gradient step
      gradstep <- step(R, clusters)
      
      R <- R - gradstep$step * gradstep$gradient
      # Try for merging
      res.merge <- merge_clusters(R, clusters, cost = L, ...)
      
      if(length(res.merge$R) != length(R)){
        R <- res.merge$R
        clusters <- res.merge$clusters
      }
      cpt <- cpt + 1 
    } 
    
    if(length(R) == 1){
      return(
        list(
          R = R,
          clusters = clusters,
          nllh = -(d - 1) * (d - 2) * R
        )
      )
    }
    return(
      list(
        R = R,
        clusters = clusters, 
        nllh = L(R, clusters)
      )
    )
  }
}



#' Search of optimal lambda for the penalty in optimization
#'
#' @param data n x d matrix : the data.
#' @param chi a positive number : tuned parameter for the exponential weights
#' @param l_grid a numerix vector of psoitive number : the grid line for lambda.
#' @param include_zero Boolean : if FALSE (default) avoid computation for  
#' non-penalization setting (i.e. lambda = 0).
#'
#' @returns Returns the optimal optimization from get_clusters() with the best
#'  lambda in the grid line. 
#'
#' @examples
best_clusters <- function(data, chi, l_grid, include_zero = FALSE){
  # Initialization 
  Gamma_est <- emp_vario(data)
  d <- ncol(data)
  R.init <- Gamma2Theta(Gamma_est)
  
  # Exponential weights construction 
  D <- D_tilde2_r(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)
  
  for(k in 1:(d - 1)){
    for(l in (k + 1):d){
      W[k, l] <- exp(-chi * D(k, l))
    }
  }
  W <- W + t(W)
  
  if(include_zero){
    lambda <- c(0, l_grid[l_grid != 0])
  }else{
    lambda <- l_grid[l_grid != 0]
  }
  
  # for one grid lambda
  if(length(lambda) == 1){
    Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = lambda)
    res_base <- Cluster_HR(R.init, it_max = 200)
    return(
      list(
        R = res_base$R,
        clusters = res_base$clusters,
        nllh = res_base$nllh,
        lambda = lambda
      )
    )
  }
  
  l_opt <- lambda[1]
  
  # First estimation 
  Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = l_opt)
  res_base <- Cluster_HR(R.init, it_max = 200)

  
  # Search of an optimal penalty
  for(i in 2:length(lambda)){
    Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = lambda[i])
    res <- Cluster_HR(R.init, it_max = 200)
    if(res$nllh < res_base$nllh){
      res_base <- res
      l_opt <- lambda[i]
    }
  }
  
  return(
    list(
      R = res_base$R,
      clusters = res_base$clusters,
      nllh = res_base$nllh,
      lambda_optim = l_opt
    )
  )
  
}

##============================= Simulation study ===============================
#' Give the lambda for a list of optimization results from replication.
#'
#' @param list List of results from best_clusters with include_zero = FALSE.
#'
#' @returns A tibble giving the optimal lambda (non null) for each replication.
#'
#' @examples
extract_lambda <- function(list){
  n <- length(list)
  lambda <- rep(NA, n)
  for(i in 1:n){
    lambda[i] <- list[[i]]$lambda_optim
  }
  return(
    tibble(simulation = 1:n, lambda_opt = lambda)
  )
}

#' Extract negative loglikelihood for best penalised results and with no penalty.
#'
#' @param list_pen List of results from best_clusters with include_zero = FALSE.
#' @param list_nopen List o results from best_clusters with one size grid equal to 0.
#'
#' @returns A tibble of the negative loglikelihood in both situation, for each 
#' replication.
#'
#' @examples
extract_nllh <- function(list_pen, list_nopen){
  n <- length(list_pen)
  nllh_pen <- rep(NA, n)
  nllh_nopen <- rep(NA, n)
  for(i in 1:n){
    nllh_pen[i] <- list_pen[[i]]$nllh
    nllh_nopen[i] <- list_nopen[[i]]$nllh
  }
  return(
    tibble(simulation = 1:n, nllh_pen = nllh_pen, nllh_nopen = nllh_nopen)
  )
}

#' Computation of the Adjusted Rand Index for two clusters
#'
#' @param cluster1 a list of vector : the first cluster.
#' @param cluster2 a list of vector : the second cluster.
#'
#' @returns The ARI for the two clusters.
#'
#' @examples
ARI <- function(cluster1, cluster2){
  U1 <- U_matrix(cluster1)
  U2 <- U_matrix(cluster2)
  # Contingency table 
  N <- t(U1) %*% U2
  n <- sum(N)
  
  # total rows and columns
  a <- rowSums(N)
  b <- colSums(N)
  
  num <- sum(N * (N - 1) / 2) - (sum(a * (a - 1) / 2) * sum(b * (b - 1) / 2)) / (n * (n - 1) / 2)
  denom <- 0.5 * (sum(a * (a - 1) / 2) + sum(b * (b - 1) / 2)) - (sum(a * (a - 1) / 2) * sum(b * (b - 1) / 2)) / (n * (n - 1) / 2)
  return(
    num / denom
  )
}

#' Computation of the Rand Index for two clusters
#'
#' @param cluster1 a list of vector : the first cluster.
#' @param cluster2 a list of vector : the second cluster.
#'
#' @returns The RI for the two clusters.
#'
#' @examples
RI <- function(cluster1, cluster2){
  U1 <- U_matrix(cluster1)
  U2 <- U_matrix(cluster2)
  n <- sum(U1)
  a <- 0
  b <- 0 
  for(i in 1:(n - 1)){
    for(j in ((i + 1):n)){
      if(sum(U1[i, ] * U1[j ,]) * sum(U2[i, ] * U2[j ,])){
        a <- a + 1 
      }else if((1 - sum(U1[i, ] * U1[j ,])) * (1 - sum(U2[i, ] * U2[j ,]))){
        b <- b + 1
      }
    }
  }
  return(
    2 * (a + b) / n / (n - 1)
  )
}


#' Build the RI coefficients tibble for each replication
#'
#' @param cluster_init a list of vector : the true cluster.
#' @param list List of results from best_clusters with include_zero = FALSE.
#'
#' @returns A tibble with Ri and ARI values (according to the true cluster) for 
#' each replication.
#'
#' @examples
get_rand_index <- function(cluster_init, list){
  n <- length(list)
  ari <- rep(NA, n)
  ri <- rep(NA, n)
  for(i in 1:n){
    ari[i] <- ARI(cluster_init, list[[i]]$clusters)
    ri[i] <- RI(cluster_init, list[[i]]$clusters)
  }
  
  return(
    tibble(simulation = 1:n, RI = ri, ARI = ari)
  )
}

#' Summary table results for replications.
#'
#' @param list_pen List of results from best_clusters with include_zero = FALSE.
#' @param list_nopen List o results from best_clusters with one size grid equal to 0.
#' @param cluster_init a list of vector : the true cluster.
#'
#' @returns A tibble which corresponds to the inner join of all the tibble from 
#' extract_lambda, extract_nllh and get_rand_index functions.
#'
#' @examples
get_info_replicate <- function(list_pen, list_nopen, cluster_init){
  d_lambda <- extract_lambda(list_pen)
  d_nllh <- extract_nllh(list_pen, list_nopen)
  d_RI <- get_rand_index(cluster_init, list_pen)
  
  return(
    d_lambda |> inner_join(d_nllh) |> inner_join(d_RI)
  )
}

#------------------------------- Others functions ------------------------------

## Decomposition of the gradient computation using finite difference

#' Function applied after one step in one coordinate.
#'
#' @param f a function which takes values in the set of real number.
#' @param theta a matrix dxd.
#' @param h a positive number : the step size.
#'
#' @returns Return a function giving the value of the function f when we apply 
#' one step in one coefficient of the symmetric matrix theta.
#'
#' @examples
#' f <- \(.) sum(.)
#' theta <- diag(1, nr = 5, nc = 5)
#' f_h <- coord_function_step(f, theta, 0.1)
#' f_h(1, 2)
coord_function_step <- function(f, theta, h){
  function(i, j){
    d <- nrow(theta)
    if(i == j){
      step <- diag(h * (1:d == i))
      return(
        f(theta + step)
      )
    }
    u_step <- matrix(ifelse((1:(d**2)%%d == i%%d)*(0:(d**2-1)%/%d + 1 == j), h, 0),
                   nc = d)
    
    step <- u_step + t(u_step) - diag(diag(u_step))
    
    return(f(theta + step))
  }
}

#' Computation of the partial derivative of a function.
#'
#' @param f a function which takes values in the set of real number.
#' @param theta a matrix dxd.
#' @param h a positive number : the step size.
#'
#' @returns Return a function giving the partial derivative of the function f 
#' for one coefficient of theta.
#' 
#' @examples
#' f <- \(.) sum(.)
#' theta <- diag(1, nr = 5, nc = 5)
#' df_ij <- coord_diff_finite(f, theta, 0.1)
#' df_ij(1, 2)
coord_diff_finite <- function(f, theta, h){
  function(i, j){
    f_h <- coord_function_step(f, theta, h)
    return((f_h(i, j) - f(theta)) / h)
  }
}

#' Computation of the gradient matrix of the function
#'
#' @param f a function which takes values in the set of real number.
#' @param theta a matrix dxd.
#' @param h a positive number : the step size.
#'
#' @returns Return the full gradient (approximated by finite difference) of the
#' function f at the point theta
#'
#' @examples
#' f <- \(.) sum(.)
#' theta <- diag(1, nr = 5, nc = 5)
#' gradient_diff_finite(f, theta, 0.1)
gradient_diff_finite <- function(f, theta, h){
  df_ij <- coord_diff_finite(f, theta, h)
  
  d <- nrow(theta)
  u_gradient <- matrix(rep(0, d * d), nc = d)
  for(i in 1:d){
    for(j in i:d){
      u_gradient[i, j] <- df_ij(i, j)
    }
  }
  gradient <- u_gradient + t(u_gradient) - diag(diag(u_gradient))
  return(gradient)
}


#' Computation of the first element of the trace's gradient
#'
#' @param gamma a dxd matrix : an estimation of the variogram gamma.
#'
#' @returns A function of clusters which returns the following matrix :
#'                         2(U' Gamma U)
#' which represents the derivative of tr(Gamma U R U') with symmetric R matrix.
#'
#' @examples
#' clusters <- list(c(1,2,3),c(4,5))
#' gamma <- matrix(1:25,nc = 5)
#' dt1 <- dtrace_1(gamma)
#' dt1(clusters)
dtrace_1 <- function(gamma){
  function(clusters){
    U <- U_matrix(clusters)
    return(2 * t(U) %*% gamma %*% U)
  }
}

#' Computation of the second element of the trace's gradient
#'
#' @param gamma a dxd matrix : an estimation of the variogram gamma.
#'
#' @returns A function of clusters which returns the following matrix :
#'                       Tilde Gamma_kl = p_k tr(Gamma_C_l) + p_l tr(Gamma_C_k)
#'                       Tilde Gamma_kk = p_k tr(Gamma_C_k) 
#' where p_k is the cardinal of cluster C_k, and Gamma_l the submatrix of Gamma
#' for the cluster C_l.
#' 
#' You can get a matrix expression of this gradient given by : 
#'                Tilde Gamma = p T^t + T p^t + diag(p*T)
#' where T is the vector of entries tr(Gamma_C_l) and p*T is the vector product 
#' componentwise.
#'
#' This derivative correspond to the linear function of r_kl : 
#'                    sum_k(p_k sum_l(r_kl tr(Gamma_l)))
#'
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' gamma <- matrix(1:25, nc = 5)
#' dt2 <- dtrace_2(gamma)
#' dt2(clusters)
dtrace_2 <-  function(gamma){
  function(clusters){
    p <- sapply(clusters, length)
    T.vector <- trace_vector(gamma, clusters)
    return(
      p %*% t(T.vector) + T.vector %*% t(p) - diag(T.vector * p)
    )
  }
}

#' Compute the full gradient of the trace part
#'
#' @param gamma a dxd matrix : an estimation of the variogram gamma.
#'
#' @returns A function of the clusters which compute the value of the gradient
#' for the trace part of the log-likelihood.
#'
#' @examples
#' clusters <- list(c(1,2,3),c(4,5))
#' gamma <- matrix(1:25,nc = 5)
#' dt <- dtrace(gamma)
#' dt(clusters)
dtrace <- function(gamma){
  dt1 <- dtrace_1(gamma)
  dt2 <- dtrace_2(gamma)
  function(clusters){
    dt1(clusters) - dt2(clusters)
  }
}

## Computation of the gradient of the penalty

#' The composed gradient between f(R) and another function
#'
#' @param clusters a list of vector : each vector gives the element of a cluster.
#' @param A a dxd matrix : the jacobian matrix of a function g.
#' @returns If A = J(g) where J is the jacobian matrix of g, row_f_gradient compute
#' the derivative of (g o f)(R), with f the function defined in section 4.3.1 of 
#' cluster document
#'
#' @examples
#' clusters <- list(c(1,3), 2)
#' A <- matrix(c(0,2,1,2,0,4,1,4,0), nc = 3)
#' row_f_gradient(A, clusters)
row_f_gradient <- function(A, clusters){
  U <- U_matrix(clusters)
  indic <- rep(1, length(clusters))
  T.vector <- trace_vector(A, clusters)
  
  tUAU <- t(U) %*% A %*% U
  trAp <- tcrossprod(T.vector * p, indic)
  return(
    tUAU - (trAp + t(trAp)) + diag((1 + p) * T.vector - .5 * diag(tUAU))
  )
}



#--------------------------- Dead lands of functions ---------------------------

# where the functions are false or do not work correctly.



# #' Intermediate function for gradient penalty computation
# #'
# #' @param R K x K symmetric matrix.
# #' @param clusters a list of vector : each vector gives the element of a cluster.
# #' @param weights a d x d symmetric matrix with a zero diagonal.
# #'
# #' @returns A function which compute a particular scalar product which appears in
# #' the computation of the gradient of the penalty. The formula is given by : 
# #'                        sum_{k'<l} W_k'l (r_kl - r_kk')
# #'                        sum_{k<l'} W_kl' (r_kl - r_ll')
# #'                        sum_{q} W_kq (r_kl - r_kq)
# #' see equations in section 4.3.3 for details.
# #'
# #' @examples
# #' 
# #' 
# #' 
# inter_prod <- function(R, clusters, weights){
#   get_W <- weight_clustered(weights)
#   K <- length(clusters)               # Number of clusters
#   W <- get_W(clusters)                # Weights clustered
#   
#   function(k, l){
#     if(k==l){
#       return(
#         sum(W[k, ] * (R[k, k] - R[, k]))
#       )
#     }
#     if(k < l){
#       indic <- 1 * (1:K < l)
#       
#       return(sum(W[l, ] * (R[k, l] - R[, k]) * indic))
#     }else{
#       indic <- 1 * (1:K > k)
#       
#       return(sum(W[k, ] * (R[k, l] - R[, l]) * indic))
#     }
#   }
# }





# #' Computation of the penalty's gradient
# #'
# #' @param weights a d x d symmetric matrix with a zero diagonal.
# #'
# #' @returns A function of clusters and corresponding R matrix. Compute the gradient 
# #' with fixed weight. The expression of the gradient is given by :
# #'  dpen_kl = 2p_k inter(k,l) + 2p_l inter(l, k) + 2 W_kl (r_kk + r_ll - 2r_kl)
# #'  dpen_kk = 2 (p_k-1) inter(k, k)
# #' see equations in section 4.3.3 for details.
# #'
# #' @examples
# #' R <- matrix(c(0.5, -1,
# #'               -1, -1), nr = 2)
# #' clusters <- list(c(1,3), c(2,4)) 
# #' W <- matrix(c(0, 1, 1, 1,
# #'               1, 0, 1, 1,
# #'               1, 1, 0, 1,
# #'               1, 1, 1, 0), nc = 4)
# #' dpen <- penalty_grad(W)
# #' dpen(R, clusters)
# #' 
# penalty_grad <- function(weights){
#   get_W <- weight_clustered(weights)
#   function(R, clusters){
#     # Initialization
#     K <- length(clusters)               # Number of clusters
#     p <- sapply(clusters, length)       # Vector of cluster's size
#     W <- get_W(clusters)                # Weights clustered
#     
#     dpen <- matrix(rep(0, K * K), nc = K)
#     f <- inter_prod(R, clusters, weights)   # build the intermediate function
# 
#     for(k in 1:(K-1)){
#       for(l in (k+1):K){
#         dpen[k, l] <- 2 * p[k] * f(k, l) + 2 * p[l] * f(l, k)
#       }
#       dpen[k, k] <- 2 * (p[k] - 1) * f(k, k)
#     }
#     dpen[K, K] <-  2 * (p[K] - 1) * f(K, K)
#     # To get symmetry of the gradient matrix
#     dpen <- t(dpen) + dpen - diag(diag(dpen))
#     
#     # No problem for the diagonal because gamma give a zero diagonal matrix
#     return(
#       dpen + 2 * W * gamma_function(R)
#     )
#   }
# }




