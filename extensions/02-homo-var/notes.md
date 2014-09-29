## Working notes

R code:

* `mkMerMod` in `utilities.R` just copies the opt. output to the class `*merMod` object.
* Searching for code that extracts the relative covariance matrix:
  * `print.merMod` function in `lmer.R` prints the model object
    * `VarCorr.merMod` function in `lmer.R` extracts the rel. cov. matrix
    * `sigma.merMod` function in `lmer.R` extracts the residual sigma, which is further passed to that `VarCorr.merMod` function.
    
  
Cpp code:

* `merPredD::setTheta` function in `predModule.cpp` does the `mapping` job:
  * update theta
  * update Lambdat
  
