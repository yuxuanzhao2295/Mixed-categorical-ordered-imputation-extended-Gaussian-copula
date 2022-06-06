# Preparation

## Working path

Experimental codes are implemented in R. For codes to work properly, set the working directory as below:

```setwd("~/Mixed-categorical-ordered-imputation-extended-Gaussian-copula/Categorical_EGC")```

 with `~` replaced by the relative path in your place.

 ## Software

 The imputation algorithm of the extended Gaussian  copula is implemented in  `func_EGC.R`. Its implementation depends on the R package `gcimputeR` and `rootSolve`. `rootSolve` can be installed from CRAN. `gcimputeR` can be installed from Github:
  ```
library(devtools)
install_github("udellgroup/gcimputeR")
 ```
 
 The following R packages are also required for experiments: `missForest`, `missMDA`, `softImpute`, `mice`, `purrr`. All of them can be installed directly from CRAN.
 
 As discussed in Sec 2.1 of the supplement, EGC uses a Python implementation to estimate the marginal and a R implementation to estimate the copula correlation. In this repo, we use a complete R implementation of EGC for simplicity, which is slower. A complete Python implementation of EGC will be available soon.
 
# Replication Instruction

## Experiments and Results
Main results can be replicated using file `main_sim.R` and `main_realdata.R`. 

 
 