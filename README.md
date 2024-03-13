# normalCowboy
Guassian LASSO with Julia and Python API

The "chunky" version divides the covariance matrix to be fit into sections and fits seperately (fitting to residuals). This means that many interactions are fit many times, but with smaller semi-definite programs and so this is often faster. However, this will not give an equivalent result because it enforces the LASSO penalty uniformly on each chunk, rather than on the matrix as a whole.
