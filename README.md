



# SOFIA: Scalar-On-Functional Regression via Integrated Adaptive Group Penalty

This repository contains R code implementing  SOFIA: Scalar-On-Functional Regression via Integrated Adaptive Group Penalty, 
proposed in the paper "Selection of functional predictors and smooth coefficient estimation for scalar-on-function regression models", by Hedayat fathi, Marzia A. Cremona, and Federico Severino 
(https://arxiv.org/abs/2506.17773).


The method is based on projecting functional predictors into a RKHS and performing adaptive group-penalized coordinate descent for sparse estimation.


# 1. Introduction

SOFIA (Scalar-On-Functional regression via Integrated Adaptive penalty) provides a framework for analyzing models of the form:

$$
Y_i = \sum_{j=1}^p \int X_{ij}(t) \beta_j(t) dt + \varepsilon_i.
$$

where

- $Y_i$  is a scalar response,
- $X_{ij}(t)$ are functional predictors observed on a common grid,
- $\beta_j(t)$ are unknown coefficient functions,
- $\varepsilon_i$ is an error term.


Each coefficient function $\beta_j$ is represented in the eigenbasis of a kernel. 

The method is adaptive. Cross-validation is used to choose $\lambda$.

SOFIA is designed for researchers working in functional data analysis, statistics, machine learning, econometrics, and biostatistics.


---

## 2. Functions

### Kernel and Eigenbasis Construction

The code provides several kernel options and their corresponding eigen-expansions:

- Sobolev kernel  
- Exponential kernel  
- Gaussian kernel  
- Matern $3/2$ and $5/2$ kernels  
- Periodic kernel  

For each kernel, the code:

- constructs the kernel matrix on a given grid,
- computes its eigenvalues and eigenvectors,
- rescales eigenpairs to be compatible with numerical integration.

### Projection Operators

The code includes functions to:

- project functional data `X` onto an eigenbasis of a kernel (e.g. via `projection_basis()`),
- project into the kernel space `K` using eigenvalues (via `projection_K()`),
- reconstruct functions on the original domain from basis coefficients (via `projection_domain()` and `projection_domain_K()`).

### Optimization

Several norms are implemented:

The function `estimation_norm_COBYLA()` uses the `nloptr` package with the algorithm `"NLOPT_LN_COBYLA"` to solve a nonlinear optimization problem used to estimate norms appearing in the penalty.

### Penalized Regression and Variable Selection

The core estimation is based on penalized regression:

$$ \min_{\beta \in K} \frac{1}{2N} \lVert Y - X\beta \rVert^2 + \lambda \sum_{j=1}^p \omega_j \lVert \beta_j \rVert_K $$

where:

- $\lVert\beta_j \rVert_K$ is the norm of the coefficient function in the space induced by $K$,
- $\omega_j$ are weights (equal to 1 in the non-adaptive step, updated in the adaptive step),
- $\lambda$ is the tuning parameter.

The function

- `definition_beta()`

implements a coordinate-descent scheme over predictors for a given grid of `lambda` values and returns:

- estimated coefficients in the kernel basis,
- selected predictors,
- the chosen `lambda`.

An additional function

- `definition_beta_CV()`

performs cross-validation for a given lambda grid.

### Cross-Validation and Adaptive Step

The function

- `estimation_beta()`

implements the full two-step procedure:

1. **Non-adaptive step**  
   - All weights  $\omega_j = 1$.  
   - A grid of `lambda` values is defined (using `definition_lambda()` if not supplied).  
   - A training/test split is used to select the best `lambda`.  

2. **Adaptive step**  

   - Weights are updated as
     $$\omega_j = 1/ \lVert \hat{\beta}_j\lVert_K,$$
     where $\hat{\beta}_j$ is the estimate from the non-adaptive step.  
   - Only predictors selected in the first step are kept.  
   - Cross-validation is used to choose `lambda` in the adaptive step.

The output includes:

- selected predictors (non-adaptive and adaptive),
- coefficients in the kernel basis,
- coefficients reconstructed on the time grid,
- cross-validation error curves,
- prediction mean squared errors (MSE) for both stages.

---

## Examples

A simulation example `Simulation_SOFIA_n500_p10_active5_Gauss_8.r` is provided. This script reproduces one of the main simulations in the paper (n=500, p=10, 5 active predictors)

We use Gaussian kernel with parameter 8.










