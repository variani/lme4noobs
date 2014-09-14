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
## [1] '1.1.7'
```


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





