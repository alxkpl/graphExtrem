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
## Computation of the gradient of the trace split in tow elements

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




#------------------------------- Others functions ------------------------------

## Decomposition of the gradient computation using finite difference

#' Function applied after one step in one coordinate.
#'
#' @param f a function which takes values in the set of real number.
#' @param theta a matrix dxd.
#' @param h a positive number : the step size.
#'
#' @returns Return a function giving the value of the function f when we apply 
#' one step in one coefficient of the matrix theta.
#'
#' @examples
#' f <- \(.) sum(.)
#' theta <- diag(1, nr = 5, nc = 5)
#' f_h <- coord_function_step(f, theta, 0.1)
#' f_h(c(1,1))
coord_function_step <- function(f, theta, h){
  function(index){
    d <- nrow(theta)
    i <- index[1]
    j <- index[2]
    step <- matrix(ifelse((1:(d**2)%%d == i%%d)*(0:(d**2-1)%/%d + 1 == j), h, 0),
                   nc = d)
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
#' df_ij(c(1,1))
coord_diff_finite <- function(f, theta, h){
  function(index){
    f_h <- coord_function_step(f, theta, h)
    return((f_h(index) - f(theta)) / h)
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
  gradient <- matrix(rep(0, d*d), nc = d)
  for(i in 1:d){
    for(j in 1:d){
      gradient[i, j] <- df_ij(c(i,j))
    }
  }
  return(gradient)
}

