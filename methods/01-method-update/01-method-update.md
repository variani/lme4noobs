# update method of lme4 R package
Andrey Ziyatdinov  
`r Sys.Date()`  


## Include


```r
library(lme4)
library(lme4qtl)
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
##                 m1 34.70588 35.04873 35.85028 35.21293 35.43481 41.81900
##        m1.opt.none 29.51423 29.66837 30.61623 29.96862 30.25699 36.83185
##          m1.update 35.89472 36.23624 36.65655 36.42146 36.76047 39.11307
##           m1.start 34.04036 34.48707 35.89125 34.91659 37.44799 38.61308
##          m1.start2 34.01843 34.08569 35.51000 34.48333 35.02715 40.95332
##  m1.start.opt.none 28.74758 28.81002 29.97371 29.09829 31.99967 33.08563
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

![](01-method-update_files/figure-html/plot_tab-1.png)<!-- -->

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

## `relmatLmer` function

### Benchmarking 2


```r
f1 <- Reaction ~ Days + (1 | Subject)
m1 <- lmer(f1, sleepstudy)
m2 <- relmatLmer(f1, sleepstudy)

tab2 <- microbenchmark(
  m1 = lmer(f1, sleepstudy),
  m1.opt.none = lmer(f1, sleepstudy, control = lmerControl(optimizer = "none")),
  m1.start = update(m1, f1, start = start.m1),
  m2 = relmatLmer(f1, sleepstudy),
  m2.update = update(m2, f1),
  m2.start = relmatLmer(f1, sleepstudy, start = start.m1),
  m2.update.start = update(m2, f1, start = start.m1),
  m2.opt.none.start = relmatLmer(f1, sleepstudy, control = lmerControl(optimizer = "none"), start = start.m1),
  m2.update.opt.none.start = update(m2, f1, control = lmerControl(optimizer = "none"), start = start.m1),
  times = 10)
```


```r
autoplot(tab2)
```

![](01-method-update_files/figure-html/plot_tab2-1.png)<!-- -->

### Model comparison


```r
system.time(mod1 <- lmer(f1, sleepstudy))
```

```
##    user  system elapsed 
##   0.032   0.000   0.036
```

```r
system.time(mod2 <- relmatLmer(f1, sleepstudy, control = lmerControl(optimizer = "none"), start = start.m1))
```

```
##    user  system elapsed 
##   0.032   0.000   0.030
```

```r
system.time(mod3 <- update(mod1, f1, control = lmerControl(optimizer = "none"), start = start.m1))
```

```
##    user  system elapsed 
##   0.028   0.000   0.031
```

```r
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
fixef(mod3)
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

```r
getME(mod3, "theta")
```

```
## Subject.(Intercept) 
##                   1
```
