% Homogeneous covariance among random effects terms
% Andrey Ziyatdinov
% 14/09/2014







## About

This report aims to reproduce an example in the article 
[Fitting Linear Mixed-Effects Models using lme4](http://arxiv.org/abs/1406.5823), Appendix A.1.

The model under study has a form:

```
respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2)
```

The model has two random terms, both including a random intercept and slopes.
Each random term contains a 2x2 covariance matrix.
The current version of `lme4` does not allows a model, where random terms share the same covariance matrices.
This example shows such a model thanks to the modular imlementation of the `lme4` package.

The `lme4` package version:


```r
library(lme4)
packageVersion("lme4")
```

```
## [1] '1.1.8'
```

### Include


```r
library(plyr)
```


## Summary

Steps to implement an homogeneous extension of `lmer`:

* Bear in mind two models, a seed model `Hetero` and a model `Homo` to be extended.
  * On the first phase, you will modify the seed object.
    The final results comes from the optimization Step 3,
    which holds a solution, but it is not intended for exploring the model output (LRT tests, df, etc).
  * On the second phase, you will need to modify the optimized model somehow back, 
    in order to emulate the `lmer` structure, as given in the `Hetero` object.
* For a better understanding of manipulation with the model parameters,
  write down the model in a matrix form and get know all the model paraemeters.
  * For the `Hetero` model: two random effect terms (`n_trms` variable),
    each term has two variables `pRE`, each term has three parameters types are `sigma_1_i`, `sigma_2_i`, `rho_i`,
    where `i = 1, 2`. 
  * For the `Homo` model: two random effect terms (`n_trms` variable),
    each term has two variables `pRE`, two term share three model parameters types are `sigma_1`, `sigma_2`, `rho`.
* Fit the `Hetero` model and pass through the relevant matrices, looking at their structure and their dimensions.
* The **key** trick goes on Step 1 of parsing the formula.
  * See section 2.3 in the `lme4` article.
  * The changes affect three variables stored in `lfHomo$reTrms`:
      * `theta`: vector now has 3 values (not 6)
      * `lower`: vector now has 3 values (not 6)
      * `Lind`: vector now has 3 unique values (not 6)
  * The last modification `Lind` defines a new mapping operation,
    when updating the covariance factor.
* The Steps 2 and 3 (deviance function and optimization) 
  are performed as for the `Hetero` model.
* The Step 4 (output) includes two parts.
  * A new deviance function `devf0` is created for the `Hetero` model.
    The solution (`theta` values saved the `opt` object) 
    is passted to the environment of the `devf0` funcion.
  * A special S3 class `homoLmerMod` is created
    with a corrected computation `logLik`.
  

## Data simulation


```r
set.seed(1)

nGrps <- 50

explVar <- data.frame(explVar1 = rnorm(nGrps^2), 
  explVar2 = rnorm(nGrps^2))

groupFac <- expand.grid(groupFac1 = as.factor(1:nGrps),
  groupFac2 = as.factor(1:nGrps))

randomIntercept <- expand.grid(randomIntercept1 = rnorm(nGrps),
  randomIntercept2 = rnorm(nGrps))
  
rnmdSlope <- expand.grid(randomSlope1 = rnorm(nGrps),
  randomSlope2 = rnorm(nGrps)) - randomIntercept
  
linearPredictor <- apply(randomIntercept + rnmdSlope*explVar, 1, sum)

residError <- rnorm(nGrps^2)

respVar <- linearPredictor + residError

dat <- data.frame(respVar, explVar, groupFac)

dim(dat)
```

```
## [1] 2500    5
```

```r
head(dat)
```

```
##   respVar explVar1 explVar2 groupFac1 groupFac2
## 1 -1.9019  -0.6265  -1.8055         1         1
## 2  2.2711   0.1836  -0.6780         2         1
## 3 -3.1781  -0.8356  -0.4734         3         1
## 4 -3.0882   1.5953   1.0274         4         1
## 5 -0.2488   0.3295  -0.5974         5         1
## 6 -3.6511  -0.8205   1.1598         6         1
```

## Default model


```r
mod1 <- lmer(respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2), dat)
mod1
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: respVar ~ 1 + (explVar1 | groupFac1) + (explVar2 | groupFac2)
##    Data: dat
## REML criterion at convergence: 7908
## Random effects:
##  Groups    Name        Std.Dev. Corr 
##  groupFac1 (Intercept) 0.92          
##            explVar1    1.47     -0.80
##  groupFac2 (Intercept) 1.03          
##            explVar2    1.26     -0.74
##  Residual              1.01          
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##     -0.0487
```

The `theta` model parameters are (expected to be of length `6`):

```
round(mod1@theta, 2)
```

The length of the vector of the random effects is (expected to be: 2 terms x 50 group levels x (slope + intercept) = `200`):


```r
length(mod1@u)
```

```
## [1] 200
```

The relative covariance factor `Lambda`:


```r
dim(mod1@pp$Lambdat)
```

```
## [1] 200 200
```


```r
mod1@pp$Lambdat[1:5, 1:5]
```

```
## 5 x 5 sparse Matrix of class "dgCMatrix"
##                                          
## [1,] 0.9082 -1.1581 .       .      .     
## [2,] .       0.8729 .       .      .     
## [3,] .       .      0.9082 -1.1581 .     
## [4,] .       .      .       0.8729 .     
## [5,] .       .      .       .      0.9082
```


The number of non-zero elements in this matrix:


```r
length(mod1@pp$Lind)
```

```
## [1] 300
```

### Modular call


```r
parsedFormula <- lFormula(formula = respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2),
  data = dat)

devianceFunction <- do.call(mkLmerDevfun, parsedFormula)

optimizerOutput <- optimizeLmer(devianceFunction)

mod2 <- mkMerMod(rho = environment(devianceFunction),
  opt = optimizerOutput,
  reTrms = parsedFormula$reTrms,
  fr = parsedFormula$fr)
```


```r
mod2
```

```
## Linear mixed model fit by REML ['lmerMod']
## REML criterion at convergence: 7908
## Random effects:
##  Groups    Name        Std.Dev. Corr 
##  groupFac1 (Intercept) 0.92          
##            explVar1    1.47     -0.80
##  groupFac2 (Intercept) 1.03          
##            explVar2    1.26     -0.74
##  Residual              1.01          
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##     -0.0487
```

## Extension

### Modular calls

Step 1: parse formula


```r
lfHetero <- lfHomo <- lFormula(respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2), dat, REML = FALSE)
  
if(length(pRE <- unique(sapply(cnms <- lfHomo$reTrms$cnms, length))) > 1L) {
  stop("each random effects term must have the same number\n",
    "of model matrix columns for a homogeneous structure")
}
```


```r
p <- ncol(lfHomo$X)
nth <- choose(pRE + 1, 2)
n_trms <- length(cnms)
```


```r
tab <- ldply(1:3, function(pRE) data.frame(num_var_random_term = pRE, length_theta = choose(pRE + 1, 2)))
tab <- data.frame(tab, theta_names = c("sigma", "2 sigmas, 1 rho", "3 sigmas, 3 rhos"))
tab
```

```
##   num_var_random_term length_theta      theta_names
## 1                   1            1            sigma
## 2                   2            3  2 sigmas, 1 rho
## 3                   3            6 3 sigmas, 3 rhos
```


```r
lfHomo$reTrms <- within(lfHomo$reTrms, {
  theta <- theta[1:nth]
  lower <- lower[1:nth]
  Lind <- rep(1:nth, length = length(lfHomo$reTrms$Lambdat@x))
})
```


```r
tab <- data.frame(length_theta = c(length(lfHetero$reTrms$theta), length(lfHomo$reTrms$theta)), 
  length_lower = c(length(lfHetero$reTrms$lower), length(lfHomo$reTrms$lower)),
  unique_Lind = c(paste(unique(lfHetero$reTrms$Lind), collapse = ", "), 
    paste(unique(lfHomo$reTrms$Lind), collapse = ", ")))
tab
```

```
##   length_theta length_lower      unique_Lind
## 1            6            6 1, 2, 3, 4, 5, 6
## 2            3            3          1, 2, 3
```

Step 2 and 3: deviance function and its optimization.


```r
devf <- do.call(mkLmerDevfun, lfHomo)
opt <- optimizeLmer(devf)
```

Step 4: output object


```r
opt$par
```

```
## [1]  0.9647 -1.0320  0.8753
```

```r
th <- rep(opt$par, n_trms)
devf0 <- do.call(mkLmerDevfun, lfHetero)
devf0(opt$par <- th)
```

```
## [1] 7910
```

```r
opt$par
```

```
## [1]  0.9647 -1.0320  0.8753  0.9647 -1.0320  0.8753
```


```r
setClass("homoLmerMod", representation(thetaUnique = "numeric"), contains = "lmerMod")

logLik.homoLmerMod <- function(object, ...) {
  ll <- lme4:::logLik.merMod(object, ...)
  attr(ll, "df") <- length(object@beta) +
  length(object@thetaUnique) +
  object@devcomp[["dims"]][["useSc"]]
  return(ll)
}

refitML.homoLmerMod <- function(object, newresp, ...) {
  if(!isREML(object) && missing(newresp)) return(object)
  stop("can't refit homoLmerMod objects yet")
}
```

### `homoLmer` function


```r
homoLmer <- function(formula, data, use.mkMerMod = FALSE, test.envir = FALSE) 
{
  mc <- match.call()

  lfHetero <- lfHomo <- lFormula(formula, data = data, REML = FALSE)
  
  if(length(pRE <- unique(sapply(cnms <- lfHomo$reTrms$cnms, length))) > 1L) {
    stop("each random effects term must have the same number\n",
      "of model matrix columns for a homogeneous structure")
  }

  p <- ncol(lfHomo$X)
  nth <- choose(pRE + 1, 2)
  n_trms <- length(cnms)
  
  lfHomo$reTrms <- within(lfHomo$reTrms, {
    theta <- theta[1:nth]
    lower <- lower[1:nth]
    Lind <- rep(1:nth, length = length(lfHomo$reTrms$Lambdat@x))
  })
  
  devf <- do.call(mkLmerDevfun, lfHomo)
  opt <- optimizeLmer(devf)
  
  th <- rep(opt$par, n_trms)
  devf0 <- do.call(mkLmerDevfun, lfHetero)
  if(!test.envir) {
    devf0(opt$par <- th)
  } else {
    opt$par <- th
  }

  mkMerMod(rho = environment(devf0),
    opt = opt,
    reTrms = lfHetero$reTrms,
    fr = lfHetero$fr)
}
```


```r
mod4 <- homoLmer(respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2), dat)
```


```r
mod5 <- as(mod4, "homoLmerMod")
mod5@thetaUnique <- mod4@theta[1:3]
```


```r
c(length(mod4@beta) + length(mod4@theta) + mod4@devcomp[["dims"]][["useSc"]],
  length(mod5@beta) + length(mod5@theta) + mod5@devcomp[["dims"]][["useSc"]],
  length(mod5@beta) + length(mod5@thetaUnique) + mod5@devcomp[["dims"]][["useSc"]])
```

```
## [1] 8 8 5
```

The value of the likelihood is the same, and the difference between the models is in `df`.

```r
logLik(mod4)
```

```
## 'log Lik.' -3955 (df=8)
```

```r
logLik(mod5)
```

```
## 'log Lik.' -3955 (df=5)
```


```r
mod4
```

```
## Linear mixed model fit by maximum likelihood  ['lmerMod']
##      AIC      BIC   logLik deviance df.resid 
##     7926     7973    -3955     7910     2492 
## Random effects:
##  Groups    Name        Std.Dev. Corr 
##  groupFac1 (Intercept) 0.977         
##            explVar1    1.370    -0.76
##  groupFac2 (Intercept) 0.977         
##            explVar2    1.370    -0.76
##  Residual              1.013         
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##     -0.0477
```

```r
mod5
```

```
## Linear mixed model fit by maximum likelihood  ['homoLmerMod']
##      AIC      BIC   logLik deviance df.resid 
##     7920     7949    -3955     7910     2492 
## Random effects:
##  Groups    Name        Std.Dev. Corr 
##  groupFac1 (Intercept) 0.977         
##            explVar1    1.370    -0.76
##  groupFac2 (Intercept) 0.977         
##            explVar2    1.370    -0.76
##  Residual              1.013         
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##     -0.0477
```

### Testing assignment in `devfun0` environment


```r
mod6 <- homoLmer(respVar ~ 1 + (explVar1|groupFac1) + (explVar2|groupFac2), dat, test.envir = TRUE)
mod6 <- as(mod6, "homoLmerMod")
mod6@thetaUnique <- mod6@theta[1:3]
```


```r
mod5
```

```
## Linear mixed model fit by maximum likelihood  ['homoLmerMod']
##      AIC      BIC   logLik deviance df.resid 
##     7920     7949    -3955     7910     2492 
## Random effects:
##  Groups    Name        Std.Dev. Corr 
##  groupFac1 (Intercept) 0.977         
##            explVar1    1.370    -0.76
##  groupFac2 (Intercept) 0.977         
##            explVar2    1.370    -0.76
##  Residual              1.013         
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##     -0.0477
```

```r
mod6
```

```
## Linear mixed model fit by maximum likelihood  ['homoLmerMod']
##      AIC      BIC   logLik deviance df.resid 
##     7920     7949    -3955     7910     2492 
## Random effects:
##  Groups    Name        Std.Dev. Corr
##  groupFac1 (Intercept) 1.03         
##            explVar1    1.03     0.00
##  groupFac2 (Intercept) 1.03         
##            explVar2    1.03     0.00
##  Residual              1.03         
## Number of obs: 2500, groups:  groupFac1, 50; groupFac2, 50
## Fixed Effects:
## (Intercept)  
##      -0.233
```

## Observations and conclusions

* It was shown how to extend the `lmer` call via 4 modular steps.
* In general, methods of the output module (Step 4) must be revised and rewritten with a user-defined control.
  * In the `homoLmerMod` expample, the `logLik` methods does the correct job,
    that can further be used in ANOVA tests.
  * It seems that the degree of freedom is under control, but it is not so.
    For example, look at the print results for two `mod4` and `mod5` models.
    The residual `df` is `2500 - 8 = 2492` for both,
    while this value has to be `2500 - 5 = 2495` for the `mod5` model.
* Unexplained operations.
  * `devf0(opt$par <- th)`: installing the parameters in the environment of the unmodified deviance function,
  as said the `lme4` article.    
    
