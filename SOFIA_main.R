

#install.packages(c("nloptr", "grplasso","refund", "kernlab", "caret", "mgcv" , "MASS" ,"here" , "fda", "readxl"))

# 
# # Libraries
# 
# # Define the required libraries
# required_libraries <- c("nloptr", "grplasso", "refund", "kernlab", 
#                         "caret", "mgcv", "MASS", "here", "fda", "readxl", "splines")
# 
# # Install missing libraries
# for (lib in required_libraries) {
#   if (!requireNamespace(lib, quietly = TRUE)) {
#     install.packages(lib, dependencies = TRUE)
#   }
# }
# 
# # Load the libraries
# for (lib in required_libraries) {
#   library(lib, character.only = TRUE)
# }
# 
# # Confirm the loaded libraries
# cat("All required libraries are successfully installed and loaded.\n")


library("nloptr") 
library("grplasso")
library('refund')
library('kernlab')
library('caret')
library('mgcv')
library('refund')
library('MASS')
library('here')
library('fda')
library('readxl')

###############################################################


# Sobolev Kernel Generation Function

#' #@H: This part is taken from Parodi & Reimherr. 

#' @description This function generates a Sobolev kernel matrix and computes its eigenvalues and eigenvectors.
#' @param a Numeric. Lower bound of the interval.
#' @param b Numeric. Upper bound of the interval.
#' @param m Integer. Number of grid points for the kernel matrix.
#' @param sigma Numeric. Parameter controlling the shape of the kernel.
#' @param plot.eigen Logical. If TRUE, plots the cumulative sum of eigenvalues normalized by the total sum.
#' @return A list containing:
#'   \item{vectors}{A matrix whose columns are the eigenvectors of the kernel matrix.}
#'   \item{values}{A vector containing the eigenvalues of the kernel matrix.}
#' @examples
#' result <- sobolev_kernel_generation(0, 1, 100, 1.5, TRUE)
#' @export

#source("1.sobolev_kernel_generation.r")


sobolev_kernel_generation <- function( a, b, m, sigma, plot.eigen = FALSE )
{
  # Define the kernel function K1(t, s) based on Sobolev formulation
  K1<-function(t,s){
    if(t < s)
      return(sigma * cosh(sigma*(b-s))*cosh(sigma*(t-a))/sinh(sigma*(b-a)))
    else
      return(sigma * cosh(sigma*(b-t))*cosh(sigma*(s-a))/sinh(sigma*(b-a)))
  }
  
  # Generate m equally spaced points between 0 and 1
  pts<-seq(0,1,length=m)
  
  # Initialize an m x m matrix for the kernel values
  Sigma<-matrix(nrow=m,ncol=m)
  
  # Compute the symmetric kernel matrix
  for(i in 1:m){
    for(j in i:m){
      Sigma[i,j] = K1(pts[i],pts[j])
      Sigma[j,i] = Sigma[i,j]
    }
  }
  
  # Compute eigenvalues and eigenvectors of the kernel matrix
  E_Sig<-eigen(Sigma)
  VAL <- E_Sig$values # Extract eigenvalues
  V<-E_Sig$vectors  # Extract eigenvectors
  
  #plot the cumulative proportion of eigenvalues
  if (plot.eigen){ plot(cumsum(VAL)/sum(VAL))}
  
  # Return eigenvectors and eigenvalues as a list
  return(list(vectors = V, values = VAL )) 
  
}

######################################################################



# Matern Kernel Generation Function

#' @H: Modified from Mirshani, I have changed the parameter scale to be 
#' compatible with Gaussian kernel that is defined later.  

#' @description This function generates a Matern kernel matrix and computes its eigenvalues and eigenvectors.
#' @param m Integer. Number of grid points for the kernel matrix.
#' @param ord Integer. Order of the Matern kernel. Must be either 3 or 5.
#' @param domain Numeric vector. The set of points at which the kernel is evaluated.
#' @param sigma Numeric. Scale parameter controlling the shape of the kernel.
#' @param plot.eigen Logical. If TRUE, plots the cumulative sum of eigenvalues normalized by the total sum.
#' @return A list containing:
#'   \item{vectors}{A matrix whose columns are the eigenvectors of the kernel matrix.}
#'   \item{values}{A vector containing the eigenvalues of the kernel matrix.}
#' @examples
#' domain <- seq(0, 1, length = 100)
#' result <- Matern_kernel_generation(100, 3, domain, 1.5, TRUE)
#' @export



Matern_kernel_generation <- function(m, ord=3, domain, sigma, plot.eigen = FALSE )
{
  
  # # Transform sigma to match RBF behavior
  # if (ord == 3) {
  #   sigma_transformed <- (sqrt(3)/sigma )
  # } else if (ord == 5) {
  #   sigma_transformed <- (sqrt(5)/sigma )
  # } else {
  #   stop("Invalid order. Use ord = 3 or ord = 5.")
  # }
  
  
  # Transform sigma to match RBF behavior
  if (ord == 3) {
    sigma_transformed <- (sqrt(3/(2*sigma)) )
  } else if (ord == 5) {
    sigma_transformed <- (sqrt(5/(6*sigma) ))
  } else {
    stop("Invalid order. Use ord = 3 or ord = 5.")
  }
  
  
  # Define the Matern kernel function K(t, s)
  K <- function(t,s){
    if(ord == 3)
      return((1 + sqrt(3)*abs(t-s)/sigma_transformed) * exp(-sqrt(3)*abs(t-s)/sigma_transformed))
    if(ord == 5)
      return((1 + sqrt(5)*abs(t-s)/sigma_transformed + (5*(abs(t-s)^2)/(3*sigma_transformed^2))) * exp(-sqrt(5)*abs(t-s)/sigma_transformed))
  }
  
  # Compute the symmetric kernel matrix
  pts <- domain
  Sigma <- matrix(nrow=m, ncol=m)
  for(i in 1:m){
    for(j in i:m){
      Sigma[i,j] = K(pts[i], pts[j])
      Sigma[j,i] = Sigma[i,j]
    }
  }
  
  # Compute eigenvalues and eigenvectors of the kernel matrix
  E_Sig <- eigen(Sigma)
  VAL <- E_Sig$values   # Extract eigenvalues
  V <- E_Sig$vectors  # Extract eigenvectors
  
  if (plot.eigen) {plot(cumsum(VAL)/sum(VAL))}
  
  # Return eigenvectors and eigenvalues as a list
  return(list(vectors = V, values = VAL ))
}



#########################################################


# Kernel Generation Function

#' #@H: This part is taken from Parodi & Reimherr.

#' source("2.generation_kernel.r")

#' @description This function generates kernel matrices and computes their eigenvalues 
#' and eigenvectors for various kernel types.
#' @param type Character. The type of kernel to generate. Options are "sobolev", "exponential", 
#' "gaussian", "Matern5/2", and "Matern3/2".
#' @param parameter Numeric. The kernel parameter controlling the scale or smoothness.
#' @param domain Numeric vector. The set of points at which the kernel is evaluated.
#' @param thres Numeric. The threshold for selecting significant eigenvalues (default: 0.99).
#' @param return.derivatives Logical. If TRUE, returns the derivatives of eigenvectors.
#' @return A list containing:
#'   \item{eigenvect}{A matrix whose columns are the selected eigenvectors.}
#'   \item{eigenval}{A vector containing the selected eigenvalues.}
#'   \item{derivatives}{(Optional) A matrix of derivatives of eigenvectors if return.derivatives is TRUE.}
#' @examples
#' domain <- seq(0, 1, length = 100)
#' result <- generation_kernel("sobolev", 1.5, domain)
#' @export



generation_kernel <- function(type = 'sobolev', parameter = NULL, domain, thres = 0.99,
                              return.derivatives = FALSE)
{
  
  M_integ <- length(domain)/diff(range(domain)) # integration normalization factor
  
  # Validate kernel type
  if (!(type %in% c('sobolev', 'exponential', 'gaussian', 'Matern5/2', 'Matern3/2')))
  {
    stop ("error: not defined kernel, please define the set
              of eigenfunctions and eigenvectors manually")
  }
  
  # Validate parameter input
  if (length(parameter) != 1 )
  {
    stop ("please provide the parameters of the kernel. See the help page
              of the function for the parameter definition")
  }
  
  # Generate kernel eigenfunctions and eigenvalues for 5 types of kernels.
  if (type =='sobolev')
  {
    kernel_def <- sobolev_kernel_generation(a = domain[1], b = domain[length(domain)],
                                            m = length(domain), sigma = parameter,
                                            plot.eigen = FALSE)
    kernel_def$values <- kernel_def$values/M_integ
    kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)
  }
  if (type =="Matern5/2")
  {
    kernel_def <- Matern_kernel_generation(m = length(domain), ord = 5, domain = domain,
                                           sigma = parameter, plot.eigen = FALSE)
    kernel_def$values <- kernel_def$values/M_integ
    kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)
  }
  if (type =="Matern3/2")
  {
    kernel_def <- Matern_kernel_generation(m = length(domain), ord = 3, domain = domain,
                                           sigma = parameter, plot.eigen = FALSE)
    kernel_def$values <- kernel_def$values/M_integ
    kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)
  }
  if (type == 'exponential')
  {
    rbfkernel <- laplacedot(sigma = sqrt(parameter/2))
    mat <- kernelMatrix(rbfkernel, domain)
    kernel_def <- list(vectors = eigen(mat)$vectors*sqrt(M_integ),
                       values = eigen(mat)$values/M_integ)
  }
  if (type == 'gaussian')
  {
    rbfkernel <- rbfdot(sigma = parameter)
    mat <- kernelMatrix(rbfkernel, domain)
    kernel_def <- list(vectors = eigen(mat)$vectors*sqrt(M_integ),
                       values = eigen(mat)$values/M_integ)
  }
  
  # Select significant eigenvalues based on threshold
  num_eigen <- which((cumsum(kernel_def$values)/sum(kernel_def$values))>thres)[1]
  eigen_chosen <- 1:num_eigen
  
  # Extract selected eigenvalues and eigenvectors
  autoval <- kernel_def$values[eigen_chosen]
  autovett <- kernel_def$vectors[,eigen_chosen]
  
  # Compute derivatives if requested
  if (return.derivatives == TRUE)
  {
    diff_autovett <- apply(autovett, 2, diff)/(domain[2]-domain[1])
    return(list(eigenvect = autovett, eigenval = autoval,
                derivatives = diff_autovett))
  }else{
    return(list(eigenvect = autovett, eigenval = autoval))
  }
  
}

###########################################################



# Periodic Kernel Generation Function

#' #@H: This part is taken from Parodi & Reimherr.
#  source("3.generation_kernel_periodic.r")

#' @description This function generates a periodic kernel matrix and computes its eigenvalues and eigenvectors.
#' @param period Numeric. The period of the kernel, given as a proportion of the domain.
#' @param parameter Numeric. The kernel parameter controlling the smoothness.
#' @param domain Numeric vector. The set of grid points at which the kernel is evaluated.
#' @param thres Numeric. The threshold for selecting significant eigenvalues (default: 0.99).
#' @param return.derivatives Logical. If TRUE, returns the derivatives of eigenvectors.
#' @return A list containing:
#'   \item{eigenvect}{A matrix whose columns are the selected eigenvectors.}
#'   \item{eigenval}{A vector containing the selected eigenvalues.}
#'   \item{derivatives}{(Optional) A matrix of derivatives of eigenvectors if return.derivatives is TRUE.}
#' @examples
#' domain <- seq(0, 1, length = 100)
#' result <- generation_kernel_periodic(0.5, 1.5, domain)
#' @export


generation_kernel_periodic <- function(period = NULL, parameter = NULL, domain, thres = 0.99,
                                       return.derivatives = FALSE)
{
  M_integ <- length(domain)/diff(range(domain)) #integration normalization factor
  
  # Validate period input
  if (length(period) != 1 )
  {
    stop ("please provide the period of the kernel.") # period is as proportion of the domain
  }
  # Validate parameter input
  if (length(parameter) != 1 )
  {
    stop ("please provide the parameters of the kernel. See the help page
              of the function for the parameter definition")
  }
  
  # Define periodic distance-based kernel function
  dist_periodic <- function(x, y, sigma, p)
  {
    sigma^2* exp(- (2*sin(pi*abs(x-y)/(p))^2)/(sigma))
  }
  
  # Generate kernel matrix
  generate_matrix <- function(domain, sigma, p){
    sapply(domain, function(x){
      sapply(domain, function(y) {dist_periodic(x,y, sigma,p)})
    })
  }
  
  # Compute kernel matrix
  kernel_matrix <- generate_matrix(domain, sigma = parameter, p = period)
  
  kernel_def <- eigen(kernel_matrix)
  kernel_def$values <- kernel_def$values/M_integ
  kernel_def$vectors <- kernel_def$vectors*sqrt(M_integ)
  
  # Select significant eigenvalues based on threshold
  num_eigen <- which((cumsum(kernel_def$values)/sum(kernel_def$values))>thres)[1]
  eigen_chosen <- 1:num_eigen
  
  # Extract selected eigenvalues and eigenvectors
  autoval <- kernel_def$values[eigen_chosen]
  autovett <- kernel_def$vectors[,eigen_chosen]
  
  # Compute derivatives if requested
  if (return.derivatives == TRUE)
  {
    diff_autovett <- apply(autovett, 2, diff)/(domain[2]-domain[1])
    return(list(eigenvect = autovett, eigenval = autoval,
                derivatives = diff_autovett))
  }else{
    return(list(eigenvect = autovett, eigenval = autoval))
  }
  
}


#########################################################

# Projection onto Basis of H = L^2

#' #@H: This part is modification of Parodi & Reimherr. 

#' @description This function projects a data matrix onto an eigenvector basis in L^2 space (It you need other spaces modify that).
#' @param X Numeric matrix of dimensions N x T. The data to be projected, 
#' where N is the number of observations, and T is the number of time points.
#' @param eigenvect Numeric matrix of dimensions T x J. The eigenvector basis, where J is the number of basis functions.
#' @param M_integ Numeric. Integration normalization factor.
#' @return A matrix of dimensions N x J, representing the projected data in the basis of H.
#' @examples
#' X <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' eigenvect <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
#' projection <- projection_basis(X, eigenvect, M_integ = 1)
#' @export



projection_basis <- function(X, eigenvect, M_integ)
{
  N <- dim(X)[1] # Number of observations
  
  # Initialize matrix to store inner products
  X_mat_inprods <- matrix(NA, N, dim(eigenvect)[2]) #attention: ours is the transpose of Parodi
  
  # Compute inner products (projection of X onto eigenvectors)
  for (i in 1:N)
  {
    for(j in 1:dim(eigenvect)[2])
    {
      X_mat_inprods[i,j] <- sum(eigenvect[,j] * X[i,]) / M_integ
    }
  }
  return (X_mat_inprods)
}

#############################################################


# Projection onto Space K


#' #@H: This function is by me. 

#' @description This function projects a matrix of functions (rows) onto a kernel-based eigenvector space.
#' @param X Numeric matrix of dimensions N x T, where N is the number of observations, and T is the number
#'  of time points.
#' @param eigenvect Numeric matrix of dimensions T x J. The eigenvector basis of the kernel space, 
#' where J is the number of basis functions.
#' @param eigenvalues Numeric vector of length J. The corresponding eigenvalues for the eigenvectors.
#' @param M_integ Numeric. Integration normalization factor.
#' @return A matrix of dimensions N x J, representing the projected data in the kernel space.
#' @examples
#' X <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' eigenvect <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
#' eigenvalues <- runif(10, 0.1, 1)
#' projection <- projection_K(X, eigenvect, eigenvalues, M_integ = 1)
#' @export


projection_K <- function(X, eigenvect, eigenvalues, M_integ) {
  
  N <- dim(X)[1] # Number of observations
  J <- dim(eigenvect)[2] # Number of basis functions
  
  # Initialize matrix to store projections
  X_mat_inprods <- matrix(NA, N, J)
  
  # Compute projection using the inner product defined by K
  for (i in 1:N) {
    for (j in 1:J) {
      # Projecting using the inner product defined by K
      X_mat_inprods[i, j] <- sum(eigenvect[, j] * X[i, ]) / (M_integ * sqrt(eigenvalues[j]))
    }
  }
  return (X_mat_inprods)
}

###########################################################


# Projection onto the Domain

#' #@H: This part is taken Parodi & Reimherr.

#' @description This function projects a matrix of coefficients onto the original time grid using eigenfunctions.
#' @param B Numeric matrix of dimensions J x N. The coefficients in the eigenfunction basis, 
#' where J is the number of basis functions and N is the number of observations.
#' @param eigenfun Numeric matrix of dimensions T x J. 
#' The eigenfunctions evaluated at the time grid, where T is the number of time points.
#' @return A matrix of dimensions N x T, representing the projected data over the time grid.
#' @examples
#' B <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' eigenfun <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
#' projection <- projection_domain(B, eigenfun)
#' @export

# Project coefficients onto the time grid
projection_domain <- function(B, eigenfun)
{
  B_mat <- eigenfun %*% B
  return(t(B_mat)) # Return the transposed result to match expected dimensions (N x T)
  
}

####################################################



# Projection from K to the Domain

#' #@H: This function is by me. 

#' @description This function projects basis coefficients from the space K back to the original time domain.
#' @param B Numeric matrix of dimensions J x N. The coefficients in the eigenfunction basis from K,
#'  where J is the number of basis functions and N is the number of observations.
#' @param eigenfun Numeric matrix of dimensions T x J. The eigenfunctions evaluated at the time grid, 
#' where T is the number of time points.
#' @param eigenvalues Numeric vector of length J. The corresponding eigenvalues for the eigenvectors.
#' @return A matrix of dimensions N x T, representing the reconstructed functions over the time grid.
#' @examples
#' B <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' eigenfun <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
#' eigenvalues <- runif(10, 0.1, 1)
#' projection <- projection_domain_K(B, eigenfun, eigenvalues)
#' @export


projection_domain_K <- function(B, eigenfun, eigenvalues) {
  
  B <- as.matrix(B)
  eigenfun <- as.matrix(eigenfun)
  eigenvalues <- as.numeric(eigenvalues)
  
  # Scale coefficients using square root of eigenvalues
  B_scaled <- B * sqrt(eigenvalues)
  
  # Project scaled coefficients onto the time grid
  B_mat <- eigenfun %*% B_scaled
  
  # Return the transposed result to match expected dimensions (N x T)
  return(t(B_mat))
}

################################################




# Norm of Beta Matrix in K

#' #@H: Taken from Parodi& Reimherr.

#' @description This function computes the K-norm of a matrix, where each column represents
#'  the projection of a function onto a basis. This function is only used while updating weights omega.
#'  For each column, it does not give the norm K of the column! 
#'  It actually compute the norm K of the function whose projection on basis is the column.
#' @param M Numeric matrix of dimensions J x N, where J is the number of basis functions 
#' and N is the number of observations.
#' @param eigenval Numeric vector of length J. The eigenvalues associated with the eigenvectors.
#' @return A numeric vector of length N, where each element is the K-norm of the corresponding column in M.
#' @examples
#' M <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' eigenval <- runif(10, 0.1, 1)
#' norms <- norm_matrix_K(M, eigenval)
#' @export



norm_matrix_K <- function(M, eigenval)
{
  
  # Ensure M is a matrix
  if (is.null(dim(M)))
  {
    M <- matrix(M, length(M), 1)
  }
  
  # Compute the K-norm for each column
  norm <- apply(M, 2, function(x){sqrt(sum(x^2/eigenval))} )
  return(norm)
}

###########################################################



# Norm of a matrix in H, 


#source("5.norm_matrix_H.r")

#' #@:H The same as "norm_matrix_H" in Parodi& Reimherr, but I think made a mistake and considered 
#' rows instead of columns. I fixed that.

#' @description This function computes the H-norm of each column of a given matrix.
#' @param M Numeric matrix of dimensions J x N, where J is the number of basis functions 
#' and N is the number of observations.
#' @return A numeric vector of length N, where each element is the H-norm of the corresponding column in M.
#' @examples
#' M <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' norms <- norm_matrix_H_vec(M)
#' @export


norm_matrix_H_vec <- function (M)
{
  # Ensure M is a matrix
  if (is.null(dim(M)))
  {
    M <- matrix(M, length(M),1)
  }
  
  # Compute the H-norm for each column
  norm <- apply(M,2, function(x){sqrt(sum(x^2))} )
  return(norm)
}



############################################################


# H^p-Norm of a Matrix

#' #@H: This function is defined by me to find the Frobenius norm of a matrix.
#'  I want to use this function to compute the norm of (B_old - B) in the main function. 
#'  In Parodi they have used the norm of K, I used norm of H. 
#'  #@H: Compare with the function "square_norm_H_matrix" in C++ file. 
#'  
#' @description This function computes the Frobenius norm (H^p-norm) of a matrix.
#' @param M Numeric matrix. The input matrix for which the norm is computed.
#' @return A single numeric value representing the H^p-norm of the matrix.
#' @examples
#' M <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' norm_value <- norm_matrix_H(M)
#' @export



norm_matrix_H <- function (M){
  
  # Initialize norm value
  norm_matrix <- 0
  
  # Compute the squared sum of all elements
  for (i in 1:ncol(M)) {
    norm_matrix <- norm_matrix + sum((M[, i] * M[, i]))
  }
  
  # Return the square root of the sum (Frobenius norm)
  return(sqrt(norm_matrix))
}


############################################################


# K^p-Norm of a Matrix

#' #@H: It is used in determining the norm of (B_old - B).
#' 
#' @description This function computes the K^p-norm of a matrix, summing the squared K-norms of its columns.
#' @param B_here Numeric matrix. The input matrix for which the K^p-norm is computed.
#' @param tau Numeric vector. The eigenvalues of the kernel associated with the matrix.
#' @return A single numeric value representing the K^p-norm of the matrix.
#' @examples
#' B_here <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' tau <- runif(10, 0.1, 1)
#' norm_value <- norm_K_matrix(B_here, tau)
#' @export


norm_K_matrix <- function(B_here, tau) {
  
  norm_matrix <- 0
  
  # Compute the squared sum of all elements, weighted by tau
  for (i in 1:ncol(B_here)) {
    norm_matrix <- norm_matrix + sum((B_here[, i] * B_here[, i]) / tau)
  }
  
  return(sqrt(norm_matrix))
}



########################################################


#norm in K of K(x)



#' @description Computes the norm in K of K(x), given the components in the kernel basis 
#' and the eigenvalues of the kernel.
#' @param B Numeric matrix or vector. The components of K(x) in the kernel basis.
#' @param tau Numeric vector. The eigenvalues of the kernel.
#' @return A single numeric value representing the norm in K of K(x).
#' @examples
#' B <- rnorm(10)
#' tau <- runif(10, 0.1, 1)
#' norm_value <- norm_K_Kx(B, tau)
#' @export


norm_K_Kx <- function(B, tau) {
  sqrt(sum(B * B * tau))
}

#######################################################


# Nonlinear Optimization for Beta Norm Estimation (NLOPT_LN_COBYLA)



#' #@H: This part is taken from Parodi & Reimherr. have changed a part of that in "tot" part. 
#' @description This function estimates the norm of beta coefficients using nonlinear optimization with the COBYLA algorithm.
#' @param lambda Numeric. The regularization parameter.
#' @param omega Numeric. A weight parameter in the optimization.
#' @param B Numeric matrix. The matrix of beta coefficients.
#' @param tau Numeric vector. The eigenvalues associated with the kernel.
#' @return A numeric value representing the estimated norm.
#' @examples
#' lambda <- 0.1
#' omega <- 1
#' B <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' tau <- runif(10, 0.1, 1)
#' norm_estimate <- estimation_norm_COBYLA(lambda, omega, B, tau)
#' @export


estimation_norm_COBYLA <- function(lambda, omega, B, tau)
{
  
  # Define the optimization function  
  optimization_function <- function(x, lambda, omega, B, tau)
  {
    numerat <- (B^2) * (tau)
    denom <- (tau * x + lambda * omega)^2
    
    # This one seems to be more reasonable 
    tot <- 1 - (sum(numerat/denom))
    
    # # The following is the one in FLAME
    #tot <- 1 - 1 / (sum(numerat/denom))
    
    return(abs(tot))
  }
  
  # Set up optimization options for COBYLA
  opts= list("algorithm"= "NLOPT_LN_COBYLA", "xtol_rel"=1.0e-16)
  
  # Perform optimization using the COBYLA algorithm
  ott_model <- nloptr(0, optimization_function, opts=opts, lb=0, omega=omega,lambda=lambda, B=B, tau=tau);
  
  return(ott_model$solution)
  
}

##########################################################




# #Definition of the vector of all the possible values of lambda



#' #@H: Translated from C++.
#' @description This function generates a sequence of lambda values using a logarithmic scale, 
#' based on the maximum norm of the projected predictors.
#' @param ratio_lambda Numeric. The ratio between the smallest and largest lambda values.
#' @param num_lambda Integer. The number of lambda values to generate.
#' @param omega Numeric vector. The weight coefficients associated with each predictor.
#' @param tau Numeric vector. The eigenvalues of the kernel.
#' @param B_ls Numeric matrix. The matrix of beta functions in the kernel space, 
#' where columns represent different betas.
#' @return A numeric vector of length `num_lambda`, representing the sequence of lambda values.
#' @examples
#' ratio_lambda <- 0.01
#' num_lambda <- 100
#' omega <- runif(10, 0.1, 1)
#' tau <- runif(10, 0.1, 1)
#' B_ls <- matrix(rnorm(10 * 10), nrow = 10, ncol = 10)
#' lambda_values <- definition_lambda(ratio_lambda, num_lambda, omega, tau, B_ls)
#' @export



definition_lambda <- function(ratio_lambda, num_lambda, omega, tau, B_ls) {
  
  num_pred <- ncol(B_ls) # Number of predictors (or betas)
  lambda_max_vec <- numeric(num_pred) # Initialize vector for storing max lambda values per predictor
  
  # Compute max lambda value for each predictor
  for (j in 1:num_pred) {
    B_temp <- B_ls[, j]
    lambda_max_vec[j] <- norm_K_Kx(B_temp, tau) / omega[j]
  }
  
  # Define lambda max and min values. 
  lambda_max <- 2 * max(lambda_max_vec)
  lambda_min <- lambda_max * ratio_lambda
  
  # Compute logarithmic range for lambda values
  tau_max <- log10(lambda_max)
  tau_min <- log10(lambda_min)
  
  # Generate equidistant lambda values in log scale
  tau_vect <- numeric(num_lambda)
  for (i in 1:num_lambda) {
    tau_vect[i] <- tau_max - (i-1) * (tau_max - tau_min) / (num_lambda-1)
  }
  
  lambda <- 10^tau_vect
  
  
  # #uncomment these lines to define an equispaced vector (non logarithmic)
  # lambda <- numeric()
  # for (i in 1:num_lambda)
  # {
  #   lambda[i] = lambda_max - (i-1)*(lambda_max-lambda_min)/(num_lambda-1);
  # }
  
  
  # #@H: another exponential way:
  # exponent <- 2 # Example value, can be adjusted as needed
  # tau_vect <-numeric(num_lambda)
  # for (i in 1:num_lambda) {
  #   ratio <- (num_lambda - i) / (num_lambda - 1)
  #   tau_vect[i] <- tau_min + (tau_max - tau_min) * ratio^exponent
  #   
  # }
  # lambda <- 10^tau_vect
  
  
  return(lambda)
}

######################################################




# Main Function for Variable Selection and Coefficient Estimation


#' @description This function performs variable selection and coefficient estimation for functional predictors.
#'  predictors are projected onto the space H ( They are not on the time grid).
#' @param proj_X List of length p = number of predcitors. Each element is a numeric matrix of dimensions N x J, 
#' representing the projection of each functional predictor into the space H.
#' @param Y Numeric matrix of dimensions N x 1. The response variable.
#' @param Beta Numeric matrix of dimensions J x p (optional). Initial values for the coefficient matrix. 
#' Defaults to NULL. If you want to start with specific beta, start with that one.
#' @param domain Numeric vector of length T. The domain of functional predictors.
#' @param tau Numeric vector of length J. The eigenvalues of the kernel matrix.
#' @param eigenvectors Numeric matrix of dimensions T x J. The eigenvectors of the kernel.
#' @param omega_vect Numeric vector of length p. The penalty weights for each predictor.
#' @param N_iter Integer. The maximum number of iterations.
#' @param computation_norm Function. The function used to compute the norm of estimated betas.
#' @param thres Numeric. The convergence threshold for stopping iterations.
#' @param non_zeros_pred Integer. The maximum number of non-zero predictors used for kill switch.
#' @param lambda_start Numeric vector (optional). Predefined sequence of lambda values (it you want to give the lamdas manually).
#'  Defaults to NULL. If it remains null the code will select lambdas automatically. 
#' @param ratio_lambda Numeric (optional). The ratio for computing lambda values. Defaults to 0.
#' @param num_lambda Integer (optional). The number of lambda values to compute. Defaults to 0.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @return A list containing:
#'   \item{Beta}{A numeric matrix of dimensions J x p, representing the estimated coefficient matrix.}
#'   \item{Pred}{A numeric vector of selected predictor indices.}
#'   \item{Number_non_zeros}{Integer. The number of non-zero predictors.}
#'   \item{estimated_lambda}{Numeric. The selected lambda value.}
#'   \item{Lambda_vect}{Numeric vector of computed lambda values.}
#' @examples
#' proj_X <- list(matrix(rnorm(100 * 10), nrow = 100, ncol = 10))
#' Y <- matrix(rnorm(100), nrow = 100, ncol = 1)
#' domain <- seq(0, 1, length = 10)
#' tau <- runif(10, 0.1, 1)
#' eigenvectors <- matrix(rnorm(10 * 10), nrow = 10, ncol = 10)
#' omega_vect <- runif(10, 0.1, 1)
#' result <- definition_beta(proj_X, Y, domain = domain, tau = tau, eigenvectors = eigenvectors,
#'                           omega_vect = omega_vect, N_iter = 100, computation_norm = function(a, b, c, d) 1,
#'                           thres = 1e-6, non_zeros_pred = 1, verbose = TRUE)
#' @export




definition_beta <- function(proj_X, Y, Beta = NULL, domain, tau, 
                            eigenvectors, omega_vect, N_iter, computation_norm,
                            thres, non_zeros_pred, lambda_start = NULL,
                            ratio_lambda = 0, num_lambda = 0, verbose) {
  
  X <- proj_X # List of projected predictors, each a matrix of N x J
  num_pred = length(X) # Number of predictors
  num_data = dim(Y)[1] # Number of observations
  
  M_integ <-  length(domain)/diff(range(domain)) # Normalization for integration
  
  # Initialize coefficient matrix B (J x p)
  if (is.null(Beta)) {
    B <- matrix(0, nrow= length(tau), ncol=length(X))
  } else {
    B <- Beta
  }
  
  # Initialize least-squares estimation matrix B_ls (J x p)
  B_ls <- matrix(0, nrow= length(tau), ncol=length(X))
  
  # Initialize residual matrix E (N x 1)
  E <- Y 
  
  # Compute initial estimates of coefficients using least squares
  for (j in 1:num_pred) {
    B_ls[,j] <- (t(X[[j]]) %*%Y ) / num_data
  }
  
  # Generate lambda values if not provided
  if (is.null(lambda_start)) {
    lambda <- definition_lambda(ratio_lambda, num_lambda, omega_vect, tau, B_ls)
  } else {
    lambda <- lambda_start
  }
  
  #num_zeros_pred = 0
  E_glob <- matrix(0, nrow = nrow(E), ncol = ncol(E)) # global residuals (N x 1)
  B_old <- matrix(nrow = nrow(B), ncol = ncol(B)) # Beta of the previous iteration
  lambda_iter <- NULL
  vector_non_zeros_pred <- numeric(num_pred) #Indeces of selected predcitors
  number_lambda_computed <- length(lambda) 
  
  
  # Loop over lambda values
  for (l in 1:length(lambda)) {
    lambda_iter <- lambda[l]
    if (verbose) { message("Lambda = ", lambda_iter)}
    
    for (k in 1:N_iter){
      if (verbose) { message("*", appendLF = FALSE) }
      
      
      #@Parodi: E_glob <- Y - B %*% t(X)
      
      
      for (n in 1:num_data) {
        inner_product_sum <- 0
        for (i in 1:num_pred) {
          inner_product_sum <- inner_product_sum + X[[i]][n, ] %*% B[,i]
        }
        
        E_glob[n,] <- Y[n,] - inner_product_sum
      }
      
      B_old <- B
      num_zeros_pred <- 0 # number of NOT selected predcitors
      
      for (j in 1:num_pred) {
        omega <- omega_vect[j]
        
        X_j <- X[[j]]
        B_j <- B_old[, j]
        
        # proj_X_j <- projection_basis(X_j, eigenvectors , M_integ)
        
        for (i in 1:num_data) {
          #E [i,1] <- E_glob[i,1] + integral_over_domain(projection_domain(B_j, eigenvectors) * X_j[i,], domain)
          E [i,1] <- E_glob[i,1] + X_j[i,]%*% B_j
        }
        
        # #@p: E <- E_glob + B_j %*% t(X_j)
        
        
        # B_tilde <- proj_X_j  %*% E  / num_data
        B_tilde <- t(X_j)  %*% E  / num_data
        
        if (norm_K_Kx(B_tilde, tau) <= lambda_iter * omega) {
          
          B_j <- rep(0, length(B_j))
          
          num_zeros_pred <- num_zeros_pred + 1
        } else {
          vector_non_zeros_pred[j - num_zeros_pred] <- j
          norm_B_j <- 0
          norm_B_j <- computation_norm(lambda_iter, omega, B_tilde, tau)
          B_j <- (tau * B_tilde * norm_B_j) / (tau * norm_B_j + lambda_iter * omega)
        }
        B[, j] <- B_j
        
      } # End of p-loop
      
      error_beta_matrix <- B_old - B
      
      # # @H : We can use norm of H istead of K.
      
      # if (norm_matrix_H(error_beta_matrix) < thres) {
      if (norm_K_matrix(error_beta_matrix, tau) < thres) {
        
        # cat("\nNot significant change of Beta => Stop iteration at iter =", k + 1, "\n")
        # cat("norm_K_Kx(errr, tau) =", "\n")
        # print(norm_K_Kx(norm_K_matrix(error_beta_matrix, tau), tau))
        #k <- N_iter
        break
      }
    } # end of the lambda loop (N)
    
    if (num_zeros_pred != num_pred){
      vector_non_zeros_pred <- vector_non_zeros_pred [1:(num_pred - num_zeros_pred)]
    }  else{
      vector_non_zeros_pred <- NULL}
    
    if (verbose) {
      cat("\nnumber of non zero predictors fitted with lambda", lambda_iter,
          "is", num_pred - num_zeros_pred, "\n" )  }
    
    if (num_pred - num_zeros_pred >= non_zeros_pred) {
      # Number of non-zero predictors reached
      number_lambda_computed <- l 
      cat("Number of non-zero predictors reached!", num_pred - num_zeros_pred , "\n")
      
      break  
    } 
  }#end of the lambda loop (l)
  
  
  lambda_subset <- lambda[1 : number_lambda_computed] 
  
  
  result <- list(
    Beta = B, 
    Pred = vector_non_zeros_pred,
    Number_non_zeros = num_pred - num_zeros_pred,
    estimated_lambda = lambda_iter,
    Lambda_vect = lambda_subset
  )
  
  return(result)
  
}


###################################################################


# Cross validation


definition_beta_CV <- function(proj_X_train, Y_train, proj_X_test, Y_test, 
                               domain, tau, eigenvectors, omega_vect, N_iter,
                               computation_norm, thres, non_zeros_pred, 
                               lambda, n_fold, verbose) {
  
  
  #M_integ <-  length(domain)/diff(range(domain))
  
  Y <- Y_train
  X <- proj_X_train
  
  num_data <- dim(Y)[1]
  num_pred <-length(X)
  
  B <- matrix(0, nrow= length(tau), ncol=num_pred)
  
  Beta_best <- NULL
  
  E <- Y  # initialization of E
  E_glob <- matrix(0, nrow = nrow(E), ncol = ncol(E))
  B_old <- B
  
  lambda_iter <- NULL
  vector_non_zeros_pred <- numeric(num_pred) 
  num_zeros_pred <- 0
  
  number_lambda_computed <- length(lambda) 
  error_lambda <-numeric(number_lambda_computed)
  error_def <- 100000000
  
  
  for (l in 1:length(lambda)) {
    
    lambda_iter <- lambda[l]
    
    if (verbose) { cat("Lambda =", lambda_iter, "\n")   }
    
    for (k in 1:N_iter) {
      #if (verbose) {  cat("*", appendLF = FALSE)  }
      
      # #This part is by Me
      # 
      # S <- matrix(0, nrow = num_data, ncol = 1)
      # for (i in 1:num_data) {
      #   for (j in 1:num_pred) {
      #     S[i,1] <- S[i,1] + integral_over_domain (projection_domain(B[,j], eigenvectors) * X[[j]][i,], domain)   }
      #   E_glob[i,1] <- Y[i,1] - S[i,1]
      # }
      
      
      # proj_X <- list() 
      # 
      # for (i in 1:length(X)) {
      #   proj_X[[i]] <- projection_basis(X[[i]], eigenvectors, M_integ) # takes a list and gives a list
      #    
      #     }
      
      
      
      for (n in 1:num_data) {
        inner_product_sum <- 0
        for (i in 1:num_pred) {
          inner_product_sum <- inner_product_sum + X[[i]][n, ] %*% B[,i]
        }
        
        E_glob[n,] <- Y[n,] - inner_product_sum
      } 
      
      
      B_old <- B
      
      # Initialize variables for tracking zero predictors
      num_zeros_pred <- 0
      vector_non_zeros_pred <- matrix(0, nrow =length(num_pred), ncol = 1 ) 
      
      for (j in 1:length(X)) {
        omega <- omega_vect[j]
        
        # X_j <- X[[j]]
        
        X_j <- X[[j]]
        B_j <- B[, j]
        
        for (i in 1:num_data) {
          #E [i,1] <- E_glob[i,1] + integral_over_domain(projection_domain(B_j, eigenvectors) * X_j[i,], domain)
          # proj_X_j <- projection_basis(X_j, eigenvectors , M_integ)
          # E [i,1] <- E_glob[i,1] + proj_X_j[,i]%*% B_j
          E [i,1] <- E_glob[i,1] + X_j[i,]%*% B_j
        }
        
        #E <- E_glob + B_j %*% t(X_j)
        
        # B_tilde <- proj_X_j  %*% E / num_data
        B_tilde <- t(X_j)  %*% E  / num_data
        # B_tilde <- projection_basis((t(E) %*% X_j)/ num_data, eigenvectors, M_integ)
        
        # B_tilde <- E %*% X_j / num_data
        
        
        # Check if the norm condition is met
        if (norm_K_Kx(B_tilde, tau) <= lambda_iter * omega) {
          B_j <- numeric(length(B_j))
          num_zeros_pred <- num_zeros_pred + 1
        } else {
          vector_non_zeros_pred[j - num_zeros_pred] <- j
          
          norm_B_j <- computation_norm(lambda_iter, omega, B_tilde, tau)
          
          B_j <- tau * B_tilde * norm_B_j / (tau * norm_B_j + lambda_iter * omega)
        }
        
        B[, j] <- B_j
      } # End of the predictors loop.
      
      error_beta_matrix <- B_old - B
      #if (verbose) {cat( "\n" ,"error_beta_matrix =", error_beta_matrix, "\n")   }
      #vector_non_zeros_pred <- matrix(0, nrow =non_zeros_pred, ncol = 1 )
      
      # @H : We can use the norm of H instead of norm of K.
      if (norm_K_matrix(error_beta_matrix, tau) < thres) {
        #if (norm_matrix_H(error_beta_matrix) < thres) {  
        
        #k <- N_iter 
        break
      }
    } #End of N-loop
    
    if (verbose) {
      cat("\nnumber of non zero predictors fitted with lambda ", lambda_iter,
          " is ", num_pred - num_zeros_pred, "\n")
    }
    
    # Calculate and store the error
    
    
    Y_pred_subset <- matrix(0, nrow = dim(Y_test)[1], 1)
    
    # #@H: one way is to prject B on the domain and compute the integral. Uncomment if you want this approach: 
    
    #proj_B <- projection_domain(B , eigenvectors)
    # S <- matrix(NA, nrow = num_pred, length(T_domain))
    #   for (i in 1:dim(Y_test)[1]) {
    #     for (j in 1:num_pred) {
    #       S[j,] <- proj_B[j,] * X_test[[j]][i,]
    #       Y_pred_subset[i,] <-  integral_over_domain (S, T_domain)
    #     }
    #   }
    
    for (n in 1:dim(Y_test)[1]) {
      Y_pred_subset[n,] <- 0
      for (j in 1:num_pred){
        Y_pred_subset[n,] <- Y_pred_subset[n,] + B [,j] %*% proj_X_test[[j]][n,]
      }
    }
    #Y_pred_subset <- B %*% t(X_test)
    
    
    difference <- Y_pred_subset - Y_test
    
    # #@H: Similar to Gertheiss
    error_lambda[l] = (sum((difference) ^ 2))/dim(Y_test)[1]
    
    # #  RMSE (Root Mean Square Error)
    # square_error <- (difference)^2
    # mean_square_error <- mean(square_error)
    # error_lambda[l]<- sqrt(mean_square_error)
    
    
    # # adjusted RSS
    #error_lambda[l] = sqrt(squared_euclidean_norm(difference))/(num_pred - num_zeros_pred + 1)
    
    # # AIC
    # error_lambda[l] = 2 * (num_pred - num_zeros_pred) + num_data * log(squared_euclidean_norm(difference)/num_data)
    
    #error_lambda[l] <- norm_K_matrix(difference, tau)
    
    if (verbose) { cat("error_Lambda =", error_lambda[l], "\n") }
    
    
    # Update the best result if necessary
    if (l!=1){
      
      if (error_lambda[l] < error_def ) {
        Beta_best <- B
        error_def = error_lambda[l]
      }else
      {
        #std::cout<<"not updated beta"<<std::endl;
        
      }
    }
    
  } # End of lambda loop
  
  
  
  # Return the results in a list
  result <- list(
    Beta = Beta_best,
    error = error_lambda
  )
  
  return(result)
}


############################################################



# estimation beta

estimation_beta <- function (X, Y, domain, eigenval, eigenvectors, NoI, thres,
                             number_non_zeros, ratio_lambda, lambda_start, 
                             number_lambda, proportion_training_set, n_fold, ordered,
                             verbose = verbose) { 
  
  function_computation <- estimation_norm_COBYLA
  
  weights = rep(1, length(X))
  
  # First step: all the weights are set to 1
  
  num_lambda_NONad <- number_lambda
  
  #num_lambda_NONad <- round(number_lambda/5) # for the first step number of lambdas are reduced.
  
  
  if (verbose) print('Non-adaptive step')
  
  M_integ <- length(domain)/diff(range(domain))
  
  lambda <- numeric()
  
  
  proj_X <- list() 
  for (i in 1:length(X)) {
    proj_X[[i]] <- projection_basis(X[[i]], eigenvectors, M_integ) # takes a list and gives a list
    #proj_X[[i]] <- scale(proj_X[[i]], center = TRUE, scale = TRUE)
  }
  
  estimation_first <- definition_beta(proj_X, Y, Beta = NULL, domain, 
                                      tau = eigenval, eigenvectors = eigenvectors,
                                      omega_vect = weights,
                                      N_iter =NoI,
                                      computation_norm = function_computation,
                                      thres = thres,
                                      non_zeros_pred = number_non_zeros,
                                      lambda_start = lambda_start,
                                      ratio_lambda = ratio_lambda,
                                      num_lambda = num_lambda_NONad,
                                      verbose = verbose)
  
  lambda_first <- estimation_first$Lambda_vect
  #cat("estimation_first prdictors =", estimation_first$pred, "\n" )
  
  cat("lambda_ first = ", lambda_first, "\n")
  
  if (verbose) print('Validation on the test set: identification of lambda')
  
  
  # # definition of X_train, Y_train, X_test, Y_test
  subset <- c(rep(2, floor(dim(Y)[1]*proportion_training_set)),
              rep(1, floor(dim(Y)[1]*(1-proportion_training_set))))
  
  if (ordered == FALSE)
  {
    set.seed(123)
    random_groups <- subset[sample(1:dim(Y)[1])] # definition of the
    # training and test set 
  }
  else {random_groups <- subset}
  
  i <- 1
  
  left_out <- which(random_groups==i)
  
  
  proj_X_train <- lapply(proj_X, function(x) x[-left_out,  , drop = FALSE])
  Y_train <- Y[-left_out, , drop = FALSE]
  proj_X_test <- lapply(proj_X, function(x) x[ left_out, , drop = FALSE])
  Y_test <- Y[left_out, , drop = FALSE]
  
  
  # fitting of the model with the proportion_training_set% of the data and
  # computation of the CV error
  
  estimation_first_CV <- definition_beta_CV(proj_X_train = proj_X_train,
                                            Y_train = Y_train,
                                            proj_X_test = proj_X_test,
                                            Y_test = Y_test,
                                            domain = domain, 
                                            tau =eigenval, 
                                            eigenvectors =eigenvectors, 
                                            omega_vect = weights, 
                                            N_iter = NoI,
                                            computation_norm = function_computation,
                                            thres = thres, 
                                            non_zeros_pred = number_non_zeros,
                                            lambda = lambda_first, 
                                            n_fold = n_fold,
                                            verbose = verbose)
  
  lambda_first_selected <- which.min(estimation_first_CV$error)
  # optimum lambda. It minimizes the CV error
  
  if (verbose) {print(paste("lambda_first_selected = ", lambda_first_selected) )}
  
  if (verbose) { print(paste("Non adaptive step: final estimation with the optimum lambda")) }
  
  
  estimation_first_definite <- definition_beta(proj_X, Y, 
                                               Beta = NULL, #estimation_first_CV$Beta,
                                               domain,
                                               eigenval, eigenvectors, weights, NoI,
                                               function_computation, thres,
                                               number_non_zeros, 
                                               lambda_first[lambda_first_selected],
                                               verbose = verbose)
  
  
  
  beta_no_adaptive <- estimation_first_definite$Beta
  predictors_no_adaptive <- estimation_first_definite$Pred
  errors_non_adaptive <- estimation_first_CV$error
  
  if (is.null(dim(beta_no_adaptive))) {
    beta_2 <- matrix(beta_no_adaptive, nrow = length(beta_no_adaptive), ncol = 1)
  }
  # Compute predictions and MSE for the non-adaptive step
  if (length(predictors_no_adaptive) > 0) {
    pred_Y_non_adaptive <- sapply(1:dim(Y)[1], function(i) {
      sum(sapply(1:length(predictors_no_adaptive), function(j) {
        as.numeric(proj_X[[predictors_no_adaptive[j]]][i, , drop = FALSE] %*% 
                     beta_no_adaptive[, predictors_no_adaptive[j], drop = FALSE])
      }))
    })
    mse_non_adaptive <- mean((Y - pred_Y_non_adaptive)^2)
  } else {
    mse_non_adaptive <- NA
  }
  
  #if(verbose) {print (paste("beta_selected", matrix(estimation_first_definite$Beta[,estimation_first_definite$Pred],
  #                        length(estimation_first_definite$Beta[,estimation_first_definite$Pred]), 1)))}
  
  if(verbose) {print (paste("length(estimation_first_definite$Pred) = ", 
                            length(estimation_first_definite$Pred)))} 
  
  if(verbose) {print (paste("Adaptive Step: update of the estimation of the",
                            length(estimation_first_definite$Pred),
                            "non zeros predictors identified."))}
  
  # isolation of the significant predictors fitted by the Non-Adaptive step.
  
  if (length(estimation_first_definite$Pred)==0)
  {
    print('No significant predictors indentified.')
    result <- list(beta=NULL,
                   beta_no_adaptive=NULL,
                   predictors=NULL,
                   predictors_no_adaptive=NULL,
                   mse_non_adaptive = NULL)
    return(result)
  }else{
    X2_data <- list()
    if(length(estimation_first_definite$Pred)==1)
    {
      beta_selected <- matrix(estimation_first_definite$Beta[,estimation_first_definite$Pred],
                              length(estimation_first_definite$Beta[ ,estimation_first_definite$Pred]), 1)
      X2_data[[1]] <- X[[estimation_first_definite$Pred]]
      
    }else
    {
      beta_selected <- estimation_first_definite$Beta[,estimation_first_definite$Pred]
      X2_data <- X[estimation_first_definite$Pred]
    }
    
    
    
    ######### Adaptive step ######## 
    
    
    ## defintion of the weights
    
    weights_new <- 1/norm_matrix_K(beta_selected, eigenval)
    
    #weights_new <- 1 / norm_matrix_H_vec(beta_selected)
    
    scales_adaptive <- list()
    
    proj_X2 <- list() 
    for (i in 1:length(X2_data)) {
      proj_X2[[i]] <- projection_basis(X2_data[[i]], eigenvectors, M_integ) # takes a list and gives a list
      proj_X2[[i]] <- scale(proj_X2[[i]], center = TRUE, scale = TRUE)
      scales_adaptive[[i]] <- attr(proj_X2[[i]], "scaled:scale")
      #scales_adaptive[[i]] <- apply(proj_X2[[i]], 2, sd)
    }
    
    estimation_second <- definition_beta (proj_X2, Y, Beta = NULL, domain, 
                                          tau = eigenval, eigenvectors, weights_new, NoI,
                                          function_computation, thres, 
                                          number_non_zeros, lambda_start = lambda_start, 
                                          ratio_lambda, number_lambda,  verbose = verbose)
    
    lambda_second <- estimation_second$Lambda_vect
    
    # Cross Validation procedure
    if (verbose) print('Validation on the test set: identification of lambda')
    
    
    
    flds <- createFolds(1:dim(Y)[1], k = n_fold, list = TRUE, returnTrain = FALSE)
    
    error_f <- rep(0,length(lambda_second))
    
    
    for( f in 1:n_fold){
      
      
      left_out <- flds[[f]]
      
      proj_X2_train <- lapply(proj_X2, function(x) x[-left_out,  , drop = FALSE])
      Y_train <- Y[-left_out, , drop = FALSE]
      proj_X2_test <- lapply(proj_X2, function(x) x[ left_out, , drop = FALSE])
      Y_test <- Y[left_out, , drop = FALSE]
      
      estimation_second_CV <- definition_beta_CV (proj_X_train = proj_X2_train, 
                                                  Y_train = Y_train,
                                                  proj_X_test = proj_X2_test,
                                                  Y_test = Y_test,
                                                  domain = domain,
                                                  tau= eigenval, 
                                                  eigenvectors= eigenvectors, 
                                                  omega_vect = weights_new, 
                                                  N_iter =NoI,
                                                  computation_norm = function_computation,
                                                  thres = thres, 
                                                  non_zeros_pred = number_non_zeros,
                                                  lambda = lambda_second, 
                                                  n_fold = n_fold,
                                                  verbose = verbose)
      
      E_se <- estimation_second_CV$error
      
      #print(estimation_second_CV$error)
      
      E_se [E_se ==0] <- Inf
      error_f <- error_f + E_se
      
      print(paste("fold",f,"   ,"))
      
    }
    
    Mean_Error_CV <- error_f / n_fold
    ind_min_cv <- which(Mean_Error_CV == min(Mean_Error_CV[Mean_Error_CV > 0]))[1]
    #ind_min_cv <- which(Mean_Error_CV == min(Mean_Error_CV))[1]
    #min_pred_CV <- min(Mean_Error_CV)
    lambda_second_selected <- lambda_second[ind_min_cv]
    
    
    #lambda_second_selected <- which.min(estimation_second_CV$error)
    
    if (lambda_second_selected == lambda_second[length(lambda_second)])
    {
      warning('Last value of the gird of lambda selected by Cross Validation. Change the grid for lambda and run FLAME again!')
    }
    
    # print(lambda_second_selected)
    # plot(lambda_first,estimation_first_CV$error, pch=19, main='First step: selection of the value of lambda')
    # points(lambda_first[lambda_first_selected], estimation_first_CV$error[lambda_first_selected], pch=19, col=2)
    
    if(verbose)
    {
      print(paste("Non adaptive step: final estimation with the optimum lambda")) }
    
    estimation_second_definite <- definition_beta(proj_X2, Y,
                                                  Beta = NULL, #estimation_second_CV$Beta,  
                                                  domain,  
                                                  eigenval,eigenvectors, weights_new,
                                                  NoI, function_computation,
                                                  thres, number_non_zeros,
                                                  lambda_second[ind_min_cv],
                                                  ratio_lambda ,num_lambda ,  
                                                  verbose = TRUE)
    # }
    
    estimation_second <- estimation_second_definite
    
    predictors_2=estimation_second$Pred # final set of
    # predictors different from 0 (among the ones isolated in the
    # First: Non-Adaptive step)
    beta_2=estimation_second$Beta[,predictors_2] # estimated betas,
    # it is the matrix of the coefficients of the basis expansion
    # of the betas (with respect to the basis defined by the eigenvectors
    # of the kernel)
    
    predictor_def <- estimation_first_definite$Pred[predictors_2] # index of
    # the non-zero predictors computed in the estimation
    
    
    if (is.null(dim(beta_2))) {
      beta_2 <- matrix(beta_2, nrow = length(beta_2), ncol = 1)
    }
    
    # Compute predictions and MSE for the adaptive step
    if (length(predictors_2) > 0) {
      pred_Y_adaptive <- sapply(1:dim(Y)[1], function(i) {
        sum(sapply(1:length(predictors_2), function(j) {
          # Ensure dimensions are properly matched
          as.numeric(proj_X2[[j]][i, , drop = FALSE] %*% 
                       beta_2[, j,drop = FALSE])
        }))
      })
      mse_adaptive <- mean((Y - pred_Y_adaptive)^2)
    } else {
      mse_adaptive <- NA
    }
    
    beta_def <- matrix(0, dim(estimation_first_definite$Beta)[1],
                       dim(estimation_first_definite$Beta)[2])
    beta_def_std <- matrix(0, dim(estimation_first_definite$Beta)[1],
                           dim(estimation_first_definite$Beta)[2])
    beta_def[, predictor_def] <- beta_2 # final matrix of the estimated betas
    # (still as coefficients of the kernel basis)
    
    
    
    for (j in predictor_def) {
      beta_def_std[, j] <-beta_def[, j]
      beta_def[, j] <- beta_def[, j] / scales_adaptive[[which(predictor_def == j)]]
    }
    
    
    if (is.null(beta_def)) {
      beta_on_time_grid <- matrix(0, nrow= length(domain), ncol= length(predictors_2))
      beta_std_on_time_grid <- matrix(0, nrow= length(domain), ncol= length(predictors_2))
    }else {
      beta_on_time_grid <- t(projection_domain(beta_def, eigenvectors))
      #beta_on_time_grid <- t(projection_domain_K(beta_def, eigenvectors, eigenval)) # this one consideres @B in K
      beta_std_on_time_grid <- t(projection_domain(beta_def_std, eigenvectors))
    }
    
    
    if (verbose) {print(paste("Total number of non zeros predictor estimated is ", length(predictors_2)))}
    
    result <- list(beta_no_adaptive=estimation_first_definite$Beta,
                   predictors_no_adaptive=estimation_first_definite$Pred,
                   errors_non_adaptive = data.frame(lambda = lambda_first, error = errors_non_adaptive),
                   beta_adaptive=beta_def,
                   predictors_adaptive =predictor_def,
                   errors_adaptive = data.frame(lambda = lambda_second, error = Mean_Error_CV),
                   MSEY_nonadaptive = mse_non_adaptive, 
                   MSEY_adaptive = mse_adaptive,
                   beta_on_time_grid= beta_on_time_grid,
                   beta_std_on_time_grid= beta_std_on_time_grid
    )
    return(result)
    
  }
}


#########################################################





#' @title SOFIA: Scalar-On-Functional regression via Integrated Adaptive group penalty
#'
#' @description #' This is the main function that runs the SOFIA procedure 
#' 
#' @param X is a list of lists, where each element has:
#'   \itemize{
#'     \item \code{$time_domain}: numeric vector of length \eqn{T}, the common grid.
#'     \item \code{$data}: numeric matrix of dimension \eqn{N \times T}, with rows
#'     equal to subjects and columns to time points. If the dimensions are not
#'     \eqn{N \times T}, the function attempts to transpose the matrix.
#'   }
#' @param Y Numeric vector or matrix of dimension \eqn{N \times 1}.  The scalar response.
#' @param type_kernel Character. Type of kernel used to define the basis.
#'   One of \code{"sobolev"}, \code{"exponential"}, \code{"gaussian"},
#'   \code{"Matern5/2"}, \code{"Matern3/2"}, or \code{"periodic"}.
#' @param param_kernel Numeric. Kernel parameter controlling smoothness / scale.
#' @param thres_eigen Numeric. Cumulative proportion of explained variance used to truncate the kernel eigenexpansion
#'   (default \code{0.99}).
#' @param period_kernel Numeric (optional). Period for the periodic kernel.
#' @param NoI Integer. Maximum number of coordinate-descent iterations for the
#'   inner optimization.
#' @param thres_CD Numeric. Convergence threshold for the change in the
#'   coefficient matrix (in K-norm) in the coordinate descent.
#' @param lambda_start Numeric vector (optional). A pre-specified grid of
#'   \eqn{\lambda} values. If \code{NULL}, the function builds its own grid
#'   using \code{definition_lambda}.
#' @param number_non_zeros Integer (optional). Kill-switch parameter: maximum
#'   number of non-zero predictors allowed in the model. If \code{NULL}, it is
#'   set equal to the total number of predictors.
#' @param ratio_lambda Numeric. Ratio between smallest and largest \eqn{\lambda}
#'   in the automatically generated grid (default \code{0.01}).
#' @param number_lambda Integer. Number of \eqn{\lambda} values on the grid
#'   (default \code{100}).
#' @param proportion_training_set Numeric in \eqn{(0,1)}. Proportion of the data
#'   used for the training set in the first (non-adaptive) CV step
#'   (default \code{0.75}).
#' @param n_fold Integer. Number of folds used in the adaptive CV step
#'   (default \code{5}).
#' @param ordered Logical. If \code{TRUE}, the training / test split for the
#'   non-adaptive CV uses the original ordering of subjects. If \code{FALSE},
#'   the split is randomized (default \code{FALSE}).
#' @param verbose Logical. If \code{TRUE}, prints progress information.
#'
#' @return A list with the following elements:
#' \item{beta} a list with:
#'   \itemize{
#'     \item \code{time_domain}: grid over which the coefficient functions are evaluated.
#'     \item \code{data}: matrix of dimension \eqn{T \times p}, where each column
#'     is an estimated coefficient function on the grid.
#'   }
#' \item{beta_no_adaptive}{Matrix of dimension \eqn{J \times p} containing the
#'   estimated coefficients in the kernel basis from the non-adaptive step.}
#' \item{beta_std_on_time_grid}{Matrix of dimension \eqn{T \times p} with the
#'   standardized coefficient functions on the time grid.}
#' \item{predictors_adaptive}{Indices of predictors selected in the adaptive step.}
#' \item{predictors_no_adaptive}{Indices of predictors selected in the non-adaptive step.}
#' \item{errors_non_adaptive}{Data frame with columns \code{lambda} and
#'   \code{error} (CV errors for the non-adaptive step).}
#' \item{errors_adaptive}{Data frame with columns \code{lambda} and
#'   \code{error} (CV errors for the adaptive step).}
#' \item{selected_lambda_non_adaptive}{Value of \eqn{\lambda} minimizing the
#'   non-adaptive CV error.}
#' \item{selected_error_non_adaptive}{Minimum non-adaptive CV error.}
#' \item{selected_lambda_adaptive}{Value of \eqn{\lambda} minimizing the
#'   adaptive CV error.}
#' \item{selected_error_adaptive}{Minimum adaptive CV error.}
#' \item{MSEY_nonadaptive}{Mean squared prediction error on \eqn{Y} for the
#'   non-adaptive model.}
#' \item{MSEY_adaptive}{Mean squared prediction error on \eqn{Y} for the
#'   adaptive model.}
#'



SOFIA <- function(X, Y, type_kernel = 'sobolev', param_kernel = 8,
                  thres_eigen = 0.99, period_kernel = NULL,
                  NoI = 10, thres_CD = 0.01, lambda_start = NULL, 
                  number_non_zeros = NULL, ratio_lambda = 0.01,
                  number_lambda = 100, proportion_training_set = 0.75, n_fold = 5,
                  ordered = FALSE,
                  verbose = FALSE){
  
  X_full <-list()
  #X_full <- replicate(0, length(X), simplify = FALSE)
  T_domain <- list()
  M_integ <- NULL
  
  for (i in 1:length(X)){
    if (class(X[[i]])[1] =='fd')    {
      T_domain[[i]] = sort(c(X[[i]]$basis$rangeval, X[[i]]$basis$params))
      X_full[[i]] <- t(eval.fd(T_domain, X[[i]], Lfdobj=0, returnMatrix=TRUE))
    } else
    {
      if (sum(names(X[[i]]) == c('time_domain')) + sum(names(X[[i]]) == c('data')) == 2)
      {
        
        T_domain[[i]] <- X[[i]]$time_domain
        X_full[[i]] <- X[[i]]$data
        if (dim(X_full[[i]])[2] != length(T_domain[[i]]))
        {
          if (dim(X_full[[i]])[1] != length(T_domain[[i]]))
          {
            warning('The data matrix should be Nobs x Tdomain, so it is transposed.')
            X_full[[i]] <- t(X_full[[i]])
          }else
          {
            stop('Number of columns different from the length of the domain.' )
          }
        }
      }
    }
    
  }
  
  M_integ  <- length(T_domain[[1]])/diff(range(T_domain[[1]]))
  
  if (type_kernel == 'periodic')
  {
    if (is.null(period_kernel))
    {
      stop('periodic kernel chosen, but no period provided')
    }
    
    kernel_here <- generation_kernel_periodic(period = period_kernel,
                                              parameter = param_kernel,
                                              domain = T_domain[[1]],
                                              thres = thres_eigen,
                                              return.derivatives = FALSE)
    eigenval <- kernel_here$eigenval
    eigenvect <- -kernel_here$eigenvect
    derivatives <- kernel_here$derivatives
  }else
  {
    if ((type_kernel != "exponential") & (type_kernel != "sobolev") & (type_kernel != "gaussian") &
        (type_kernel != "Matern5/2")& (type_kernel != "Matern3/2"))
    {
      stop('provide a valid value for the type of kernel')
    }
    
    kernel_here <- generation_kernel(type = type_kernel,
                                     parameter = param_kernel,
                                     domain = T_domain[[1]],
                                     thres = thres_eigen,
                                     return.derivatives = FALSE)
    eigenval <- kernel_here$eigenval
    eigenvect <- kernel_here$eigenvect
    derivatives <- kernel_here$derivatives
  }
  
  if (verbose)
  {
    cat(paste('Estimation: identification of the contribution of the ',
              length(X_full), ' predictors on the ', dim(Y)[1], ' observations', sep =''))} 
  
  for (i in 1:length(X_full)) {
    # check dimentions
    if (dim(X_full[[i]])[1] != dim(Y)[1])
    {
      stop('Number of rows of X doesn\'t coincide with the number of individuals')
    }
  }
  if (is.null(number_non_zeros))
  {
    number_non_zeros = length(X)
  }
  
  
  Our_est <- estimation_beta (X = X_full, 
                              Y = Y,
                              domain = T_domain[[1]],
                              eigenval = eigenval, # basis
                              eigenvectors = eigenvect,
                              NoI = NoI, # max. num. iterations coordinate descent
                              thres = thres_CD, # thres for the CV
                              number_non_zeros = number_non_zeros, # kill switch parameter
                              lambda_start = lambda_start,
                              ratio_lambda = ratio_lambda, # ratio for the min. lambda
                              number_lambda = number_lambda, # num. elements of the grid for lambda
                              proportion_training_set = proportion_training_set, # training set
                              n_fold = n_fold,
                              ordered = ordered,
                              verbose = FALSE ) # no show of all the iterations
  
  # if (is.null(Our_est$beta_adaptive)) {
  #   beta_on_time_grid <- matrix(0, nrow= length(T_domain[[1]]), ncol= length(X_full))
  # }else {
  #   #beta_on_time_grid <- t(projection_domain(Our_est$beta_adaptive, eigenvect))
  #   beta_on_time_grid <- t(projection_domain_K(Our_est$beta_adaptive, eigenvect, eigenval)) # this one consideres @B in K
  # }
  
  
  beta_on_time_grid <- Our_est$beta_on_time_grid
  beta_std_on_time_grid <- Our_est$beta_std_on_time_grid
  
  if (class(X[[1]]) == 'fd')
  {
    
    basis_beta <- create.bspline.basis(norder=4,
                                       rangeval = range(T_domain),
                                       nbasis = 10)
    
    beta <- Data2fd(T_domain, beta_on_time_grid,
                    basisobj = basis_beta, nderiv = 0)
  }else
  {
    beta <- vector('list', 2)
    beta[[1]] <- T_domain[[1]]
    beta[[2]] <- beta_on_time_grid
    names(beta) <- c('time_domain', 'data')
    
    min_error_row_non_adaptive <- Our_est$errors_non_adaptive[which.min(Our_est$errors_non_adaptive$error), ]
    min_error_lambda_non_adaptive <- min_error_row_non_adaptive$lambda
    min_error_non_adaptive <- min_error_row_non_adaptive$error
    
    min_error_row_adaptive <- Our_est$errors_adaptive[which.min(Our_est$errors_adaptive$error), ]
    min_error_lambda_adaptive <- min_error_row_adaptive$lambda
    min_error_adaptive <- min_error_row_adaptive$error
    MSEY_nonadaptive <- Our_est$MSEY_nonadaptive
    MSEY_adaptive <- Our_est$MSEY_adaptive
  }
  
  return(list(
    beta = beta, 
    beta_no_adaptive = Our_est$beta_no_adaptive,
    beta_std_on_time_grid= Our_est$beta_std_on_time_grid,
    predictors_adaptive = Our_est$predictors_adaptive,
    predictors_no_adaptive = Our_est$predictors_no_adaptive,
    errors_non_adaptive = Our_est$errors_non_adaptive,
    errors_adaptive = Our_est$errors_adaptive,
    selected_lambda_non_adaptive=min_error_lambda_non_adaptive,
    selected_error_non_adaptive=min_error_non_adaptive,
    selected_lambda_adaptive=min_error_lambda_adaptive,
    selected_error_adaptive=min_error_adaptive,
    MSEY_nonadaptive = Our_est$MSEY_nonadaptive,
    MSEY_adaptive = Our_est$MSEY_adaptive
    
  ))
}

##############################################################



