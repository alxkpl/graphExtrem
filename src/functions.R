#' 3x3 variogram generator for a trivariate Husler-Reiss model
#'
#' @return A valid Husler-Reiss 3x3 variogram where the corresponding graphical model
#' doesn't have edge between nodes 1 and 2.
#'
#' @examples
#' random_Gamma12()
random_Gamma12 <- function(){
  b <- 5 * runif(1)
  c <- 5 * runif(1)
  a <- b + c 
  Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
  
  while(!(checkGamma(Gamma, returnBoolean = TRUE, alert = FALSE))){
    b <- 5 * runif(1)
    c <- 5 * runif(1)
    a <- b + c
    Gamma <- matrix(c(0, a, b, a, 0, c, b, c, 0), nr = 3)
  }
  return(Gamma)
}


#' 2-variate extremal coefficient of a bivariate Husler-Reiss distribution.
#'
#' @param gamma A real number: the value of the Husler-Reiss parameter.
#'
#' @return The value of the 2-variate extremal coefficient for a Husler-Reiss distribution
#' with parameter gamma. This value is \Lambda(1,1) = 2 phi(sqrt(gamma)/2), where 
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
#' @param Gamma A conditionnal negative definite matrix : Husler-Reiss parameter.
#'
#' @returnThe value of the 3-variate extremal coefficient for a Husler-Reiss distribution
#' with variogram Gamma.
#'
#' @examples
#' Gamma <- matrix(c(0, 5, 1, 5, 0, 4, 1, 4, 0), nr = 3)
#' lambda_2(Gamma)
lambda_2 <- function(Gamma){
  # Parameters of the variogram
  a <- Gamma[1, 2]
  b <- Gamma[1, 3]
  c <- Gamma[2, 3]
  
  # Computation of the correlation for the Gaussian distribution function
  rho_1 <- (a + b - c) / (2 * sqrt(a * b))
  rho_2 <- (a + c - b) / (2 * sqrt(a * c))
  rho_3 <- (b + c - a) / (2 * sqrt(b * c))
  
  return(
    pnorm2d(x = sqrt(a) / 2, y = sqrt(b) / 2, rho = rho_1)[1] +
      pnorm2d(x = sqrt(a) / 2, y = sqrt(c) / 2, rho = rho_2)[1] +
      pnorm2d(x = sqrt(b) / 2, y = sqrt(c) / 2, rho = rho_3)[1] 
    )
}


#' Trivariate extremal coefficient in a Husler-Reiss model.
#'
#' @param Gamma The HR variogram of a size-3 (sub-)vector . 
#'
#' @return The value of the trivariate coefficient in the corresponding HR 
#' random vector.
#'
#' @examples
#' chi_trivariate_HR(random_Gamma12())
chi_trivariate_HR <- function(Gamma){
  return(
    3 - theta(Gamma[1,2]) - theta(Gamma[1,3]) - theta(Gamma[2,3]) + lambda_2(Gamma)
    )
}



