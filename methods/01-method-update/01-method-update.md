# update method of lme4 R package
Andrey Ziyatdinov  
`r Sys.Date()`  


## Include


```r
library(lme4) 
```


```r
library(microbenchmark)
library(plyr)
library(ggplot2)
```

## Benchmarking of update methods


```r
f1 <- Reaction ~ Days + (1 | Subject)
m1 <- lmer(f1, sleepstudy)
start.m1 <- list(theta = getME(m1, "theta"))
start2.m1 <- getME(m1, c("theta", "fixef"))
m2 <- update(m1, Reaction ~ Days + (1|Subject))
m3 <- update(m1, Reaction ~ (1|Subject))
```


```r
tab <- microbenchmark(
  m1 = lmer(f1, sleepstudy),
  m1.opt.none = lmer(f1, sleepstudy, control = lmerControl(optimizer = "none")),
  m1.update = update(m1, f1, sleepstudy),
  m1.start = update(m1, f1, start = start.m1),
  m1.start2 = update(m1, f1, start = start2.m1),
  m1.start.opt.none = lmer(f1, sleepstudy, control = lmerControl(optimizer = "none"), start = start.m1),
  #m2 = lmer(Reaction ~ Days + (1 | Subject), sleepstudy),
  #m3 = lmer(Reaction ~ (1 | Subject), sleepstudy),  
  #m2.update = update(m1, Reaction ~ Days + (1|Subject)),
  #m3.update = update(m1, Reaction ~ (1|Subject)),
  times = 10)
```


```r
tab
```

```
## Unit: milliseconds
##               expr      min       lq     mean   median       uq      max
##                 m1 41.09880 41.31964 41.39687 41.38564 41.50877 41.68065
##        m1.opt.none 33.86017 33.96102 35.25992 34.29755 37.46391 38.40831
##          m1.update 41.52630 41.70893 42.94529 42.13532 45.10790 45.21790
##           m1.start 39.62557 39.71818 40.41815 39.95414 40.03770 44.03654
##          m1.start2 39.47024 39.54944 40.65966 39.75572 39.88887 49.11917
##  m1.start.opt.none 32.84936 33.01097 34.66532 33.34632 36.52874 38.57739
##  neval
##     10
##     10
##     10
##     10
##     10
##     10
```


```r
autoplot(tab)
```

![](01-method-update_files/figure-html/plot_tab-1.png) 

## Check model equvalence


### Atemp 1: main function


```r
f1 <- Reaction ~ Days + (1 | Subject)

mod1 <- lmer(f1, sleepstudy)

start.mod1 <- list(theta = getME(mod1, "theta"))
mod2 <- lmer(f1, sleepstudy, control = lmerControl(optimizer = "none"), start = start.mod1)
#mod2 <- update(mod1, f1, start = start.mod1)

mod3 <- lmer(f1, sleepstudy, control = lmerControl(optimizer = "none"))

# check
fixef(mod1)
```

```
## (Intercept)        Days 
##   251.40510    10.46729
```

```r
fixef(mod2)
```

```
## (Intercept)        Days 
##   251.40510    10.46729
```

```r
getME(mod1, "theta")
```

```
## Subject.(Intercept) 
##            1.197882
```

```r
getME(mod2, "theta")
```

```
## Subject.(Intercept) 
##                   1
```

### Atemp 2: modular call



```r
lmod <- lFormula(f1, sleepstudy)
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
mod1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

start.mod1 <- list(theta = getME(mod1, "theta"))

opt2 <- list(par = as.numeric(start.mod1$theta), fval = NA,conv = 1000, message="start copied")
mod2 <- mkMerMod(environment(devfun), opt2, lmod$reTrms, fr = lmod$fr)

# check
fixef(mod1)
```

```
## (Intercept)        Days 
##   251.40510    10.46729
```

```r
fixef(mod2)
```

```
## (Intercept)        Days 
##   251.40510    10.46729
```

```r
getME(mod1, "theta")
```

```
## Subject.(Intercept) 
##            1.197882
```

```r
getME(mod2, "theta")
```

```
## Subject.(Intercept) 
##            1.197882
```
