Algorithm
-------
sglasso: Smoothed group Lasso.


Maintainer
-------
Jin Liu   <jin.liu@duke-nus.edu.sg>


Publication
-------
Liu, J., Huang, J., Ma, S., & Wang, K. (2013). Incorporating group correlations in genome-wide association studies using smoothed group Lasso. Biostatistics, 14(2), 205-219.


Description
-------
SGL implements penalization method for group variable selection which can properly accommodate the correlation between adjacent groups. This method is based on a combination of the group Lasso penalty and a quadratic penalty on the difference of regression coefficients of adjacent groups. It encourages group sparsity and smoothes regression coefficients for adjacent groups. Canonical correlations are applied to the weights between groups in the quadratic difference penalty. 


Usage
-------
source("CallC_m_rescaled.r")
  - sglasso_sp() produces the solution path.

srouce("sp_scaled.r")
  - sglasso() solves for a fixed tuning.
