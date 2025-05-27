pacman::p_load(
  dplyr,
  ggplot2,
  ggraph,
  tidygraph,
  igraph,
  gridExtra,
  tidyr,
  evd,
  tibble,
  graphicalExtremes
)
#---------------------- Functions for travariate document ----------------------

#' 3x3 variogram generator for a trivariate Husler-Reiss model
#'
#' @param k A positive number to set the range of the coefficients.
#' @return The coefficients of a valid Husler-Reiss 3x3
#' variogram where the corresponding graphical model doesn't
#' have edge between nodes 1 and 2.
#' The vector correspond to :
#'      - a : Gamma_12
#'      - b : Gamma_13
#'      - c : Gamma_23
#' @examples
#' random_Gamma12(5)
random_Gamma12 <- function(k) {
  b <- k * runif(1)
  c <- k * runif(1)
  a <- b + c
  Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)

  while (!(graphicalExtremes::checkGamma(Gamma, returnBoolean = TRUE, alert = FALSE))) {
    b <- k * runif(1)
    c <- k * runif(1)
    a <- b + c
    Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
  }
  matrix(c(a, b, c), nc = 3)
}

#' Transform a set of variogram's parameter to a matrix.
#'
#' @param Gamma_params A numerical vector: parameters of the variogram
#'
#' @return The corresponding symetric matrix.
#'
#' @examples
to_matrix <- function(Gamma_params) {
  # Parameters of the variogram
  a <- Gamma_params[1]
  b <- Gamma_params[2]
  c <- Gamma_params[3]

  matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
}

#' 2-variate extremal coefficient of a bivariate Husler-Reiss distribution.
#'
#' @param gamma A real number: the value of the Husler-Reiss parameter.
#'
#' @return The value of the 2-variate extremal coefficient for a Husler-Reiss
#' distribution with parameter gamma. This value is
#' \Lambda(1, 1) = 2 phi(sqrt(gamma)/2), where
#' phi is the distribution function a standard gaussian.
#'
#' @examples
#' theta(1)
theta <- function(gamma) {
  2 * pnorm(sqrt(gamma) / 2)
}


#' 3-variate extremal coefficient of a trivariate Husler-Reiss distribution.
#'
#' @param Gamma_parms A vector of size 3 symbolizing a conditionnal negative
#' matrix : the Husler-Reiss parameters.
#'
#' @return The value of the 3-variate extremal coefficient for a Husler-Reiss
#' distribution with variogram represented by the vector Gamma_params.
#'
#' @examples
#' Gamma_params <- c(5, 1, 4)
#' lambda_2(Gamma_params)
lambda_2 <- function(Gamma_params) {
  # Parameters of the variogram
  a <- Gamma_params[1]
  b <- Gamma_params[2]
  c <- Gamma_params[3]

  # Computation of the correlation for the Gaussian distribution function
  rho_1 <- (a + b - c) / (2 * sqrt(a * b))
  rho_2 <- (a + c - b) / (2 * sqrt(a * c))
  rho_3 <- (b + c - a) / (2 * sqrt(b * c))

  if (sum(round(abs(c(rho_1, rho_2, rho_3)), 8) >= 1)) {
    NA
  }else {
    fMultivar::pnorm2d(x = sqrt(a) / 2, y = sqrt(b) / 2, rho = rho_1)[1] +
      fMultivar::pnorm2d(x = sqrt(a) / 2, y = sqrt(c) / 2, rho = rho_2)[1] +
      fMultivar::pnorm2d(x = sqrt(b) / 2, y = sqrt(c) / 2, rho = rho_3)[1]
  }
}


#' Trivariate extremal coefficient in a Husler-Reiss model.
#'
#' @param Gamma_params The coefficient of the HR variogram of a size-3
#' (sub-)vector.
#'
#' @return The value of the trivariate coefficient in the corresponding HR
#' random vector.
#'
#' @examples
#' chi_trivariate_HR(random_Gamma12())
chi_trivariate_HR <- function(Gamma_params) {
  3 - theta(Gamma_params[1]) - theta(Gamma_params[2]) -
    theta(Gamma_params[3]) + lambda_2(Gamma_params)
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
cond_trivariate_HR <- function(Gamma_params, i) {
  chi_trivariate_HR(Gamma_params) / (2 - theta(Gamma_params[i]))
}


#' Give a matrix of null variogram, parameterised by the diagonal values
#'
#' @param diagonal Numeric vector : diagonal's value of a matrix
#'
#' @returns A matrix which belong to the kernel of the linear application
#' gamma, computing the variogram for a given matrix. It is characterized
#' by the relation :
#'                          a_ij = (a_ii + a_jj) / 2
#' @examples
#' ker_gamma(1:4)
ker_gamma <- function(diagonal) {

  n <- length(diagonal)                 # dimension of the matrix

  one <- rep(1, n)                      # vector with only ones

  0.5 * (one %*% t(diagonal) + diagonal %*% t(one))
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
semi_def <- function(matrix) {
  !sum(Re(eigen(matrix)$values) < -1e-10)        # for numerical errors
}

#----------------------- Gradient descent for clustering -----------------------
##========= Functions for the building of some particular matrices =============

#' Computation of the matrix of clusters
#'
#' @param clusters a list of vector : each vector gives the element of a
#' cluster.
#'
#' @returns The d x K matrix U such that u_jk = 1 if the variable j belongs to
#' cluster C_k and 0 otherwise.
#'
#' @examples
#' clusters <- list(c(1,2,3),c(4,5))
#' U_matrix(clusters)
#'
U_matrix <- function(clusters) {
  K <- length(clusters)
  d <- length(unlist(clusters))
  U <- matrix(rep(0, d * K), nc = K)
  for (k in 1:K){
    for (j in 1:d){
      if (j %in% clusters[[k]]) {
        U[j, k] <- 1
      }
    }
  }
  U
}

#' From R matrix compute theta matrix
#'
#' @param R a K x K matrix : the matrix of the clusters coefficients
#' @param clusters a list of vector : each vector gives the element of a
#' cluster.
#'
#' @returns Return the theta matrix from the value of the R matrix of clusters
#' coefficient using :
#'
#'                            theta = U R U^t + A
#'
#' where U is the cluster matrix and a diagonal matrix such that the rows of
#' theta sum to zero :
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
build_theta <- function(R, clusters) {
  U <- U_matrix(clusters)
  URUt <- U %*% R %*% t(U)
  a <-  - rowSums(URUt)

  URUt + diag(a)

}


#' Compute the trace cluster matrix vector.
#'
#' @param gamma a dxd matrix : an estimation of the variogram gamma.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#'
#' @returns Returns a vector of size K of entries tr(Gamma_C_l).
#' @export
#'
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' gamma <- matrix(1:25, nc = 5)
#' trace_vector(gamma, clusters)
trace_vector <- function(gamma, clusters) {
  K <- length(clusters)
  T.vector <- rep(0, K)
  i <- 0
  for (index in clusters){
    i <- i + 1
    T.vector[i] <-  sum(diag(gamma[index, index]))
  }

  T.vector
}

#' Variogram transformation application gamma
#'
#' @param sigma A d x d numeric matrix.
#'
#' @returns For a symmetric positive matrix sigma (covariance matrix),
#' return the corresponding variogram matrix. Can be used for other
#' but with no interpretation.
#'
#' @examples
#' s_sigma <- matrix(rnorm(16, 2), nc = 4)
#' gamma_function(s_sigma %*% t(s_sigma))
gamma_function <- function(sigma) {
  indic <- rep(1, nrow(sigma))

  tcrossprod(diag(sigma), indic) + tcrossprod(indic, diag(sigma)) - 2 * sigma
}


#' Crout factorization algorithm.
#'
#' @param A a d x d symmetric positive matrix.
#' @param tol a positive value : tolerance for the zero diagonal
#'
#' @returns The LU Crout decomposition of a matrix A. For A a symmetric
#' positive matrix, there exists a LU decomposition such that :
#'                                   A = L' L
#'
#' @examples
#'A <- matrix(c(1,2,3,
#'              2,5,6,
#'              3,6,9), nc = 3)
#' L <- crout_factorisation(A)
#' L %*% t(L)
crout_factorisation <- function(A, tol = 1e-12) {
  A <- unname(as.matrix(A))
  n <- nrow(A)
  if (n != ncol(A)) {
    stop("no square matrix.")
  }
  if (!semi_def(A)) {
    stop("no positive semi definite matrix.")
  }
  d <- rep(0, n)
  L <- diag(rep(1, n))
  d[1] <- A[1, 1]
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      L[i, j] <- (1 / d[j]) * (A[i, j] - sum(L[i, ] * L[j, ] * d))
    }
    d[i] <- A[i, i] - sum(d * (L[i, ])**2)
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
psolve <- function(A, tol = 1e-12) {

  S <- crout_factorisation(A, tol = tol)

  L <- S[, which(diag(S) != 0)]         # to get no null columns

  return(
    L %*% solve(t(L) %*% L) %*% solve(t(L) %*% L) %*% t(L)
  )
}

#' Computation of the first weight matrix
#' @param data the data
#'
#' @returns The initial wieght matrix without fused variables.
#'
#' @example
#'
compute_W <- function(data) {
  Gamma_est <- graphicalExtremes::emp_vario(data)
  d <- ncol(data)
  R.init <- graphicalExtremes::Gamma2Theta(Gamma_est)
  # Exponential weights construction
  D <- D_tilde2_r(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)) {
    for (l in (k + 1):d){
      W[k, l] <- exp(-1 * D(k, l))
    }
  }
  W + t(W)
}

#' Transform a chi matrix to the corresponding variogram
#'
#' @param Chi_matrix the matrix with the chi coefficient
#'
#' @return The corresponding variogram matrix $\Gamma$.
#'
#' @example
#'
ChiToGamma <- function(Chi_matrix) {

  (2 * qnorm((2 - Chi_matrix) / 2))**2

}



#' Function to extract the coefficient of the reduced matrix R
#' 
#' @param matrix A block matrix which can be factorizable.
#' @param clusters a list of vectors : the variable index per clusters.
#' 
#' @return A matrix on size the number of clusters (length(clusters)). It
#' corresponds to the reduced matrix R of the original matrix (see clusters document
#' for definition).)
extract_R_matrix <- function(matrix, clusters) {

  K <- length(clusters)
  indx <- 1:K
  R <- matrix(rep(NA, K * K), nc = K)

  for (i in indx) {
    k <- clusters[[i]][1]
    if (length(clusters[[i]]) > 1) {
      l <- clusters[[i]][2]
      R[i, i] <- matrix[k, l]
    }
    for (j in indx[-i]){
      l <- clusters[[j]][1]
      R[i, j] <- matrix[k, l]
    }
  }

  R
}

#' Computation of the clustered weight matrix
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters : the weight matrix for each pair
#' of clusters :
#'                          W_kl = sum_(i in C_k) sum_(j in C_l) w_ij
#'
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' W <- matrix(1:25, nc = 5)
#' W_c <- weight_clustered(W)
#' W_c(clusters)
weight_clustered <- function(weights) {
  function(clusters) {
    U <- U_matrix(clusters)

    t(U) %*% weights %*% U
  }
}


#' Generalised determinant
#'
#' @param A a d x d real valued matrix.
#' @param tol a positive value.
#'
#' @returns Compute the generalised determinant of a matrix A. We recall
#' that the generalised determinant is an extension of the determinant for
#' singular matrix. It corresponds to the product of all the non zero eigen
#' values.
#'
#' @examples
#' A <- matrix(c(1,2,3,
#'               2,5,6,
#'               3,6,9), nc = 3)
#' gen_det(A)
gen_det <- function(A, tol = 1e-10) {
  res <- Re(eigen(A, only.values = TRUE)$values)      # avoid complex number (impossible for symmetric matrices)
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
#' corresponding R matrix, the value of the associated negative likelihood
#' defined by :
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
neg_likelihood <- function(gamma) {
  function(R, clusters) {
    # Building of the theta matrix from R
    theta <- build_theta(R, clusters)

    # Computation of the log-determinant part
    log_det <- log(gen_det(theta))

    # Computation of the trace part
    tr <- sum(diag(gamma %*% theta))


    - log_det - .5 * tr
  }
}


#' Gradient of the negative log likelihood without penalty
#'
#' @param gamma a d x d variogram matrix.
#'
#' @returns A function of the R matrix and clusters and compute the gradient
#' matrix of the negative log likelihood for a fixed variogram gamma. The
#' gradient matrix can be computed by :
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
nloglike_grad_np <- function(gamma) {
  function(R, clusters) {
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

    dlog + diag
  }
}


##================= Computation of the gradient of the penalty =================

#' Cluster distance squared function
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#'
#' @returns A function of cluster number : compute the square distance between
#' two clusters for the distance defined in section 4.2 in cluster document.
#'
#' @examples
#'
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' D2 <- D_tilde2_r(R, clusters)
#' D2(1, 2)
D_tilde2_r <- function(R, clusters) {
  function(k, l) {
    if (k == l) {
      # 0 distance for the same cluster even with the symmetry of the matrix
      0
    }else {
      # Parameters of clusters
      K <- length(clusters)               # Number of clusters
      p <- sapply(clusters, length)       # Vector of cluster's size

      # for fixed k, l the square difference is multiplied by 1-p_k and 1-p_l
      sum((p - ((1:K) %in% c(k, l))) * (R[k, ] - R[l, ])**2)
    }
  }
}

#' Penalty function.
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns For fixed weights, returns a function which compute the value of
#' the penalty for chosen clusters and corresponding R matrix. We recall the
#' penalty is given by :
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
penalty <- function(weights) {
  # Fixing the weights for computing weights clustered
  get_W <- weight_clustered(weights)
  function(R, clusters) {
    # Initialization
    K <- length(clusters)              # Number of clusters
    D2 <- D_tilde2_r(R, clusters)      # Function for distance between clusters
    W <- get_W(clusters)               # Weights clustered
    D <- matrix(rep(0, K * K), nc = K) # Distance matrix for clusters

    # Computation of the distance matrix
    for (l in 2:K) {
      for (k in 1:(l - 1)) {        # we keep only the lower triangular part
        D[k, l] <- D2(k, l)
      }
    }

    sum(D * W)
  }
}

#' Gradient matrix of distance between two columns
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#'
#' @returns Return a functon of indices (k', l') computing the gradient
#' matrix of tilde D^2(r_k', r_l'). See section 4.3.3 in cluster document
#' for details.
#'
#' @examples
#' R <- matrix(c(-1,0,-2,
#'               0,-3,-1,
#'               -2,-1,-1), 3)
#' clusters <- list(c(1,3), c(2), 4)
#' grad <- gradient_D2(R, clusters)
#' grad(1, 2)
gradient_D2 <- function(R, clusters) {
  # Initialization
  K <- length(clusters)                   # Number of clusters
  p <- sapply(clusters, length)           # Vector of cluster's size
  function(k, l) {
    A <- matrix(rep(0, K * K), nc = K)

    # Computation of the non zero row and column
    A[k, ] <- 2 * p * (R[k, ] - R[l, ])
    A[, l] <- 2 * p * (R[l, ] - R[k, ])
    A[k, l] <- 2 * ((p[k] - 1) * (R[k, l] - R[k, k]) +
                      (p[l] - 1) * (R[k, l] - R[l, l]))

    # Build the symmetry of the gradient (except for the diagonal)
    grad <- A + t(A)

    # Computation of the diagonal
    grad[k, k] <- 2 * (p[k] - 1) * (R[k, k] - R[k, l])
    grad[l, l] <- 2 * (p[l] - 1) * (R[l, l] - R[k, l])

    grad

  }
}


#' Computation of the penalty's gradient
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters and corresponding R matrix. Compute
#' the gradient with fixed weight. The expression of the gradient is
#' just the weighted sum of the gradient of each tilde D^2 where the
#' weights are the clustered weights.
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
penalty_grad <- function(weights) {
  get_W <- weight_clustered(weights)
  function(R, clusters) {
    # Initialization
    W <- get_W(clusters)                # Weight clustered
    K <- length(clusters)               # Number of clusters

    # Function of gradient of indices
    grad_D2 <- gradient_D2(R, clusters)

    res <- matrix(rep(0, K * K), nc = K)

    for (k in 1:(K - 1)) {
      for (l in (k + 1):K) {
        res <- res + W[k, l] * grad_D2(k, l)        # Weighted sum
      }
    }

    res
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
neg_likelihood_pen <- function(gamma, weights, lambda) {
  nllh <- neg_likelihood(gamma)
  pen <- penalty(weights)

  function(R, clusters) {

    nllh(R, clusters) + lambda * pen(R, clusters)

  }
}



##=========================== Positive condition on R ==========================

sub_theta <- function(R, clusters) {
  p <- sapply(clusters, length)           # Vector of cluster's size
  tilde_R <- - as.numeric(R %*% p)

  # Subtheta matrix :
  R %*% diag(p) + diag(tilde_R)

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
step_gradient <- function(gamma, weights, lambda, size_grid = 100) {
  # Initialization of functions
  dlog <- nloglike_grad_np(gamma)                     # Neg-lklh gradient part
  dpen <- penalty_grad(weights)                       # Penalty gradient part

  # Penalised negative log-likelihood
  nllh <- neg_likelihood_pen(gamma, weights, lambda)

  function(R, clusters) {
    # Initialization
    p <- sapply(clusters, length)           # Vector of cluster's size

    # Gradient matrix computation
    grad <- dlog(R, clusters) + lambda * dpen(R, clusters)

    # Grid line search for optimal gradient step
    # Grid line construction
    if (max(p) == 1) {
      s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters),
                     method = "Brent", lower = 0, upper = 1)$par
      while (!semi_def(sub_theta(R - s_opt * grad, clusters))) {
        s_opt <- 0.95 * s_opt
      }
      return(list(step = s_opt, gradient = grad))
    }
    s_max <- min(
      # Maximum step size to get positive matrix
      abs(((R %*% p) / (grad %*% p))[p > 1])
    )

    s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters),
                   method = "Brent", lower = 0, upper = min(s_max, 1))$par
    while (!semi_def(sub_theta(R - s_opt * grad, clusters))) {
      s_opt <- 0.95 * s_opt
    }

    # Returning results : size step and gradient matrix
    list(step = s_opt, gradient = grad)


    # Searching
    # s_opt <- 0
    # score <- nllh(R, clusters)
    #
    # for(i in 2:(length(s) - 1)){
    #   # Better step size
    #   if(score > nllh(R - s[i] * grad, clusters)){
    #     # Positive matrix checking
    #     if(semi_def(sub_theta(R - s[i] * grad, clusters))){
    #       s_opt <- s[i]
    #       score <- nllh(R - s_opt * grad, clusters)
    #     }
    #   }
    # }
  }
}

#' Function which merges clusters
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#' @param eps positive value : minimal tolerance for merging clusters
#' @param cost a function : Cost function of the optimisation
#'
#' @returns Returns, if merging, a list of the new clusters and the
#' corresponding R matrix, where the coefficient of the new clustered
#' is computing by averaging the coefficient of the two previous clusters.
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
merge_clusters <- function(R, clusters, eps = 1e-1, cost) {
  # Initialization
  D <- D_tilde2_r(R, clusters)               # Function of clusters distance
  K <- length(clusters)                      # Actual number of clusters

  # Computation of the distance matrix
  distance <- matrix(rep(Inf, K * K), nc = K)

  for (k in 1:(K - 1)) {
    for (l in (k + 1):K) {
      distance[k, l] <- D(k, l)
    }
  }

  # Search of the two potential clusters to merge
  index <- as.numeric(which(distance == min(distance), arr.ind = TRUE))
  k <- index[1]
  l <- index[2]

  # Checking uselessness of merging
  if (distance[k, l] > eps) {
    return(
      list(
        R = R,
        clusters = clusters
      )
    )
  }

  p <- sapply(clusters, length)           # Vector of cluster's size

  # Case when merging give only one cluster
  if (nrow(R) == 2) {
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

  return(
    list(
      R = R_new,
      clusters = new_clusters
    )
  )
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
get_cluster <- function(gamma, weights, lambda, ...) {
  L <- neg_likelihood_pen(gamma, weights, lambda)
  step <- step_gradient(gamma, weights, lambda, ...)
  function(R.init, it_max = 1000, eps_g = 1e-3) {
    # Initialization
    d <- nrow(gamma)
    R <- R.init
    clusters <- as.list(1:d)
    gradstep <- list(gradient = eps_g + 1)
    cpt <- 1
    while ((cpt < it_max) && (length(R) != 1) && (sum(gradstep$gradient**2) > eps_g)) {
      # Gradient step
      gradstep <- step(R, clusters)

      R <- R - gradstep$step * gradstep$gradient
      # Try for merging
      res.merge <- merge_clusters(R, clusters, cost = L, ...)

      if (length(res.merge$R) != length(R)) {
        R <- res.merge$R
        clusters <- res.merge$clusters
      }
      cpt <- cpt + 1
    }

    if (length(R) == 1) {
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
best_clusters <- function(data, chi, l_grid, it_max = 1000) {
  # Initialization
  Gamma_est <- graphicalExtremes::emp_vario(data)
  d <- ncol(data)
  R.init <- graphicalExtremes::Gamma2Theta(Gamma_est)

  # Exponential weights construction 
  D <- D_tilde2_r(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)){
    for (l in (k + 1):d){
      W[k, l] <- exp(-chi * D(k, l))
    }
  }
  W <- W + t(W)

  lambda <- l_grid[l_grid != 0]

  # for one grid lambda
  if (length(lambda) <= 1) {
    if (l_grid == 0) {
      L <- neg_likelihood(Gamma_est)
      return(
        list(
          R = R.init,
          clusters = as.list(1:d),
          nllh = L(R.init, as.list(1:d)),
          lambda = 0
        )
      )
    }
    Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = lambda)
    res_base <- Cluster_HR(R.init, it_max = it_max)
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
  res_base <- Cluster_HR(R.init, it_max = it_max)

  # Search of an optimal penalty
  for (i in 2:length(lambda)){
    Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = lambda[i])
    res <- Cluster_HR(R.init, it_max = it_max)
    if (res$nllh < res_base$nllh) {
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

##============================ Hierachy graphics ===============================
#' Compare two clusters.
#'
#' @param clusters1 The first cluster to compare.
#' @param clusters2 The second cluster to compare.
#'
#' @returns Return FALSE if the clusters have not the same number of cluster. Otherwise,
#'
#'
#' @examples
#'
compare_clusters <- function(clusters1, clusters2) {
  clust1 <- lapply(clusters1, sort)
  clust2 <- lapply(clusters2, sort)
  if (length(clust1) != length(clust2)) return(FALSE)
  all(
    sapply(
      clust1, function(x) {
        any(sapply(clust2, function(y) ifelse(length(x) == length(y), all(x == y), FALSE)))
      }
    )
  )
}


#' Detect the merge and give all information about cluster and lambda.
#'
#' @param solution_list Solution's list from best_clusters.
#'
#' @returns Return the list of the solution where a change occurs between the list of clusters.
#'
#' @examples
#'
detect_merge <- function(solution_list) {
  # Initialiser les états
  d <- sum(sapply(solution_list[[1]]$clusters, length))
  prev_clusters <- NULL
  event_list <- list(list(lambda = 0, clusters = 1:d))
  event_idx <- 2
  # Boucle pour détecter les événements
  options(warn = 1)
  for (i in 1:length(solution_list)) {
    sol <- solution_list[[i]]
    clusters <- sol$clusters

    if (!is.null(prev_clusters)) {
      if (!compare_clusters(prev_clusters, clusters)) {
        # On détecte une transition !
        event_list[[event_idx]] <- list(lambda = sol$lambda, clusters = clusters)
        event_idx <- event_idx + 1
      }
    } else {
      # Premier état
      event_list[[event_idx]] <- list(lambda = sol$lambda, clusters = clusters)
      event_idx <- event_idx + 1
    }

    prev_clusters <- clusters
  }
  options(warn = 0)
  event_list

}


#' Give the adjacency matrix coresponding to the right hierarchical according to the fusion
#' during optimization.
#'
#' @param event_list the list of solution from detect_fusion.
#' @param lambda_max The lambda which end the optimization.
#'
#' @returns A matrix which will be use for the construction of the dendrogram.
#'
get_adjacency_matrix <- function(event_list, lambda_max) {

  d <- sum(sapply(event_list[[1]]$clusters, length))

  M <- length(event_list)
  A <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(M - 1)) {
    prev <- event_list[[k]]$clusters
    nex <- event_list[[k + 1]]$clusters

    for (l in setdiff(nex, prev)) {
      for (i in prev) {
        if (length(intersect(l, i)) > 0) {
          for (j in setdiff(l, i)) {
            A[i, j] <- event_list[[k + 1]]$lambda
            A[j, i] <- event_list[[k + 1]]$lambda
          }
        }
      }
    }
  }

  A[A == 0] <- lambda_max

  diag(A) <- 0

  A

}

#' Plot the dendrogram for the optimization
#'
#' @param list_results A list of results optimization from best_clusters.
#'
#' @returns The dendrogram obtained from the simulation results for each lambda.
#'
gg_cluster <- function(list_results) {

  d <- sum(sapply(list_results[[1]]$clusters, length))
  lambda_max <- list_results[[length(list_results)]]$lambda
  options(warn = 1)
  event_list <- detect_merge(list_results)
  options(warn = 0)

  A <- get_adjacency_matrix(event_list, lambda_max)

  hclust_results <- hclust(as.dist(A), method = "average")

  hclust_results$label <- 1:d

  graph <- as_tbl_graph(hclust_results)

  # Dessiner l'arbre
  ggraph(graph, layout = "dendrogram", height = height) +
    geom_edge_elbow(linewidth = 1.5, alpha = 0.8, color = "darkorange2") +
    geom_node_text(aes(label = ifelse(leaf, label, "")), size = 4, vjust = 1.7, color = "grey40") +
    geom_node_point(color = "grey40", shape = 18, size = 3) +
    ylab(expression(lambda)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(), axis.line.y = element_line(color = "grey50"),
          plot.margin = margin(10, 10, 10, 10),
          axis.title.y = element_text(angle = 0, size = 15),
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.5),  # Couleur et taille des ticks
          axis.ticks.length = unit(0.1, "cm"))
}


#' Plot the average dendrogram for the replicate optimization
#'
#' @param replicates A list of results optimization replicates from best_clusters.
#'
#' @returns The dendrogram obtained from the simulation results for each lambda.
#'
average_hierarchy <- function(replicates) {

  d <- sum(sapply(replicates[[1]][[1]]$clusters, length))

  N <- length(replicates)

  lambda_max <- replicates[[1]][[length(replicates[[1]])]]$lambda

  list_detect <- lapply(replicates, detect_merge)

  Adj_list <- lapply(list_detect, \(.) get_adjacency_matrix(., lambda_max))

  A <- as.dist(Reduce("+", Adj_list) / N)

  hclust_results <- hclust(as.dist(A), method = "average")

  hclust_results$label <- 1:d

  graph <- as_tbl_graph(hclust_results)

  # Dessiner l'arbre
  ggraph(graph, layout = "dendrogram", height = height) +
    geom_edge_elbow(linewidth = 1.5, alpha = 0.8, color = "darkorange2") +
    geom_node_text(aes(label = ifelse(leaf, label, "")), size = 4, vjust = 1.7, color = "grey40") +
    geom_node_point(color = "grey40", shape = 18, size = 3) +
    ylab(expression(lambda)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(), axis.line.y = element_line(color = "grey50"),
          plot.margin = margin(10, 10, 10, 10),
          axis.title.y = element_text(angle = 0, size = 15),
          axis.ticks.y = element_line(color = "grey50", linewidth = 0.5),  # Couleur et taille des ticks
          axis.ticks.length = unit(0.1, "cm"))
}

##============================= Simulation study ===============================
#' Give the lambda for a list of optimization results from replication.
#'
#' @param list List of results from best_clusters with include_zero = FALSE.
#'
#' @returns A tibble giving the optimal lambda (non null) for each replication.
#'
#' @examples
extract_lambda <- function(list) {
  n <- length(list)
  lambda <- rep(NA, n)
  for (i in 1:n) {
    lambda[i] <- list[[i]]$lambda_optim
  }

  tibble(simulation = 1:n, lambda_opt = lambda)

}

#' Extract negative loglikelihood for best penalised results and with no penalty.
#'
#' @param list_pen List of results from best_clusters with include_zero = FALSE.
#' @param list_nopen List o results from best_clusters with one size grid equal
#' to 0.
#'
#' @returns A tibble of the negative loglikelihood in both situation, for each
#' replication.
#'
#' @examples
extract_nllh <- function(list_pen, list_nopen) {
  n <- length(list_pen)
  nllh_pen <- rep(NA, n)
  nllh_nopen <- rep(NA, n)
  for (i in 1:n) {
    nllh_pen[i] <- list_pen[[i]]$nllh
    nllh_nopen[i] <- list_nopen[[i]]$nllh
  }

  tibble(simulation = 1:n, nllh_pen = nllh_pen, nllh_nopen = nllh_nopen)

}

#' Computation of the Adjusted Rand Index for two clusters
#'
#' @param cluster1 a list of vector : the first cluster.
#' @param cluster2 a list of vector : the second cluster.
#'
#' @returns The ARI for the two clusters.
#'
#' @examples
ARI <- function(cluster1, cluster2) {
  U1 <- U_matrix(cluster1)
  U2 <- U_matrix(cluster2)
  # Contingency table
  N <- t(U1) %*% U2
  n <- sum(N)

  # total rows and columns
  a <- rowSums(N)
  b <- colSums(N)

  num <- sum(N * (N - 1) / 2) - (sum(a * (a - 1) / 2) * sum(b * (b - 1) / 2)) / (n * (n - 1) / 2)
  denom <- 0.5 * (sum(a * (a - 1) / 2) + sum(b * (b - 1) / 2)) - (sum(a * (a - 1) / 2) *
                                                                    sum(b * (b - 1) / 2)) / (n * (n - 1) / 2)

  # Adjusted Rand Index :
  num / denom
}

#' Computation of the Rand Index for two clusters
#'
#' @param cluster1 a list of vector : the first cluster.
#' @param cluster2 a list of vector : the second cluster.
#'
#' @returns The RI for the two clusters.
#'
#' @examples
RI <- function(cluster1, cluster2) {
  U1 <- U_matrix(cluster1)
  U2 <- U_matrix(cluster2)
  n <- sum(U1)
  a <- 0
  b <- 0

  for (i in 1:(n - 1)) {
    for (j in ((i + 1):n)) {
      if (sum(U1[i, ] * U1[j, ]) * sum(U2[i, ] * U2[j, ])) {
        a <- a + 1
      }else if ((1 - sum(U1[i, ] * U1[j, ])) * (1 - sum(U2[i, ] * U2[j, ]))) {
        b <- b + 1
      }
    }
  }
  # Rand index
  2 * (a + b) / n / (n - 1)
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
get_rand_index <- function(cluster_init, list) {
  n <- length(list)
  ari <- rep(NA, n)
  ri <- rep(NA, n)
  for (i in 1:n){
    ari[i] <- ARI(cluster_init, list[[i]]$clusters)
    ri[i] <- RI(cluster_init, list[[i]]$clusters)
  }

  tibble(simulation = 1:n, RI = ri, ARI = ari)
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
get_info_replicate <- function(list_pen, list_nopen, cluster_init) {
  d_lambda <- extract_lambda(list_pen)
  d_nllh <- extract_nllh(list_pen, list_nopen)
  d_RI <- get_rand_index(cluster_init, list_pen)

  # Final Data :
  d_lambda |>
    inner_join(d_nllh, by = join_by(simulation)) |>
    inner_join(d_RI, by = join_by(simulation)) |>
    mutate(nb_cluster = sapply(list_pen,  \(.) length(.$clusters)))
}




#' Plot some graph to summarise the results of replications
#'
#' @param data Tibble from get_info_replicate.
#' @param true_number_cluster The number of cluster from simulation.
#' @param see_title Boolean, if TRUE (default) display general title.
#'
#' @returns Plots several graphs to analyze simulation results :
#'        - a scatter plot of the Rand indexes, with their means.
#'        - the relation between lambda_opt and final number of clusters.
#'        - the repartition of the difference between penalised and non penalised
#'          negative log-likelihood.
#'        - the barplot of the optimal number of cluster.
#'
#' @examples
plot_results <- function(data, true_number_cluster, see_title = TRUE) {
  # Rand Index graph
  data |>
    pivot_longer(cols = ends_with("RI"),
                 values_to = "Values",
                 names_to = "Type") |>
    ggplot() +
    aes(y = Type, x = Values) +
    geom_vline(xintercept = 1, linetype = "dashed", col = "darkgreen") +
    geom_point(col ="grey60") +
    xlab("Index value") +
    ggtitle("Rand Index") +
    theme_ipsum(base_family = "serif") -> p

  data |>
    pivot_longer(cols = ends_with("RI"),
                 values_to = "Values",
                 names_to = "Type") |>
    group_by(Type) |> summarise(m = mean(Values)) -> data2

  p +
    geom_point(data = data2, aes(x = m, y = Type, col = "Mean Index"),
               size = 3, shape = 15, show.legend = T) +
    scale_color_manual(name = " ", values = c("Mean Index" = "darkorange")) -> p1

  # Relation between lambda_opt and number of clusters
  data |> ggplot() +
    aes(x = lambda_opt, y = nb_cluster, group = factor(nb_cluster)) +
    geom_boxplot(fill = "grey80", col = "darkorange", alpha = 0.9) +
    ggtitle(expression(paste("Relation between ", lambda["opt"], " and ", K["opt"]))) +
    labs(x = expression(lambda["opt"]), y = expression(K["opt"]))
    theme_ipsum(base_family = "serif") -> p2

  # Negative log-likelihood differences
  data |> ggplot() +
    aes(x = nllh_nopen - nllh_pen) +
    geom_density(aes(y = after_stat(density)), col = "grey40", fill = "darkorange",
                 alpha = 0.9) +
    coord_cartesian(x = c(-0.03, 0)) +
    labs(x = "Difference", y = "Density",
         title = "Negative log-likelihood difference") +
    theme_ipsum(base_family = "serif") -> p3

  # Barplot of the optimal number of cluster
  data |>
    ggplot(aes(x = nb_cluster,
               fill = ifelse(nb_cluster == true_number_cluster, "Right number", "Other"))) +
    geom_bar(alpha = 0.9) +
    scale_fill_manual(name = " ", values = c("Right number" = "darkorange")) +
    labs(x = expression(K["opt"]), y = "Number of simulations",
         title = "Final number of clusters") +
    theme_ipsum(base_family = "serif") -> p4

  merged_plot <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)

  if (see_title) {
    title <- ggdraw() +
      draw_label(paste("Optimization results for",
                       nrow(data), "replications"), fontface = 'bold', size = 16)

    return(
      plot_grid(title, merged_plot, ncol = 1, rel_heights = c(0.1, 1))
    )
  }
  return(
    merged_plot
  )
}


all_info <- function(cluster.init, list_res, lambda, one_sim = FALSE){

  if (one_sim) {
    data <- get_rand_index(cluster.init, list_res) |>
      mutate(l = lambda,
             nb_cluster = sapply(list_res, \(.) length(.$clusters)))

    return(data)
  }

  m <- length(list_res)
  data <- get_rand_index(cluster.init, list_res[[1]]) |>
    mutate(l = lambda, simulation = 1,
           nb_cluster = sapply(list_res[[1]], \(.) length(.$clusters)))
  for(i in 2:m){
    data_inter <- get_rand_index(cluster.init, list_res[[i]]) |>
      mutate(l = lambda, simulation = i,
             nb_cluster = sapply(list_res[[i]], \(.) length(.$clusters)))
    data |> add_row(data_inter) -> data
  }
  return(
    data
  )
}

plot_info <- function(results, one_sim = FALSE) {

  results |>
    group_by(l) |>
    summarise(ARI = mean(ARI)) |>
    ggplot() + aes(x = l, y = ARI) +
    geom_line(col = "grey50") +
    xlab(expression(lambda)) +
    geom_point(col = "darkorange2", shape = 20) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen") +
    ggtitle("ARI over lambda") + 
    theme_ipsum(base_family = "serif") -> p1

  results |>
    group_by(l) |>
    summarise(RI = mean(RI)) |>
    ggplot() +aes(x = l, y = RI) +
    geom_line(col = "grey50") +
    geom_point(col = "darkorange2", shape = 20) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "darkgreen") +
    xlab(expression(lambda)) +
    ggtitle("RI over lambda") +
    theme_ipsum(base_family = "serif") -> p2

  if (one_sim) {
    width <- (results$l[2] - results$l[1]) / 100
    results |>
      ggplot() + aes(x = l, y = nb_cluster) +
      geom_col(col = "darkorange1", width = width) +
      xlab(expression(lambda)) +
      ylab("K") +
      ggtitle("Number of cluster over lambda") +
      theme_ipsum(base_family = "serif") -> p3

  }else {
    results |>
      ggplot() + aes(x = l, y = factor(nb_cluster)) +
      geom_boxplot(aes(group = factor(nb_cluster)), col = "darkorange2",
                   fill = "grey80") +
      xlab(expression(lambda)) +
      ylab("K") +
      ggtitle("Number of cluster over lambda") +
      theme_ipsum(base_family = "serif") -> p3
  }


  # Graph structure
  top_row <- plot_grid(p1, p2, ncol = 2,rel_widths = c(0.5, 0.5))
  bottom_row <- plot_grid(NULL, p3, NULL, ncol = 3, rel_widths = c(0.2, 0.6, 0.2))
  final_plot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))

  return(
    final_plot
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
coord_function_step <- function(f, theta, h) {
  function(i, j) {
    d <- nrow(theta)
    if (i == j) {
      step <- diag(h * (1:d == i))
      return(
        f(theta + step)
      )
    }
    u_step <- matrix(ifelse((1:(d**2) %% d == i %% d) * (0:(d**2 - 1) %/% d + 1 == j), h, 0),
                     nc = d)

    step <- u_step + t(u_step) - diag(diag(u_step))

    # Value of the function on the step
    f(theta + step)
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
coord_diff_finite <- function(f, theta, h) {
  function(i, j) {
    f_h <- coord_function_step(f, theta, h)

    # Finite difference at coordinate (i,j)
    (f_h(i, j) - f(theta)) / h
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
gradient_diff_finite <- function(f, theta, h) {
  df_ij <- coord_diff_finite(f, theta, h)

  d <- nrow(theta)
  u_gradient <- matrix(rep(0, d * d), nc = d)
  for (i in 1:d) {
    for (j in i:d) {
      u_gradient[i, j] <- df_ij(i, j)
    }
  }
  gradient <- u_gradient + t(u_gradient) - diag(diag(u_gradient))

  # Gradient matrix
  gradient
}


#' Computation of the first element of the trace's gradient (obsolete)
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
dtrace_1 <- function(gamma) {
  function(clusters) {
    U <- U_matrix(clusters)

    # Derivative of the trace one :
    2 * t(U) %*% gamma %*% U
  }
}

#' Computation of the second element of the trace's gradient (obsolete)
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
dtrace_2 <-  function(gamma) {
  function(clusters) {
    p <- sapply(clusters, length)
    T.vector <- trace_vector(gamma, clusters)

    # Derivative of the trace two :
    p %*% t(T.vector) + T.vector %*% t(p) - diag(T.vector * p)

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
dtrace <- function(gamma) {
  dt1 <- dtrace_1(gamma)
  dt2 <- dtrace_2(gamma)
  function(clusters) {
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
row_f_gradient <- function(A, clusters) {
  U <- U_matrix(clusters)
  indic <- rep(1, length(clusters))
  p <- lapply(clusters, length)
  T.vector <- trace_vector(A, clusters)

  tUAU <- t(U) %*% A %*% U
  trAp <- tcrossprod(T.vector * p, indic)

  tUAU - (trAp + t(trAp)) + diag((1 + p) * T.vector - .5 * diag(tUAU))

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
