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
