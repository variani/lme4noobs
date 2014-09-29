% Homogeneous variance over random effects
% Andrey Ziyatdinov
% 29/09/2014







## About

This report explores an example in the article 
[Fitting Linear Mixed-Effects Models using lme4](http://arxiv.org/abs/1406.5823), Appendix A.2.

The model under study:

```
Reaction ~ Days + (Days|Subject)
```

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


## Default model

Modular call:


```r
parsedFormula <- lFormula(formula = Reaction ~ Days + (Days|Subject),
  data = sleepstudy)

devianceFunction <- do.call(mkLmerDevfun, parsedFormula)

optimizerOutput <- optimizeLmer(devianceFunction)

mod1 <- mkMerMod(rho = environment(devianceFunction),
  opt = optimizerOutput,
  reTrms = parsedFormula$reTrms,
  fr = parsedFormula$fr)
```


```r
mod1
```

```
## Linear mixed model fit by REML ['lmerMod']
## REML criterion at convergence: 1744
## Random effects:
##  Groups   Name        Std.Dev. Corr
##  Subject  (Intercept) 24.74        
##           Days         5.92    0.07
##  Residual             25.59        
## Number of obs: 180, groups:  Subject, 18
## Fixed Effects:
## (Intercept)         Days  
##       251.4         10.5
```

