# update method of lme4 R package
Andrey Ziyatdinov  
`r Sys.Date()`  


## Include


```r
library(lme4)
library(solaris)
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
##                 m1 48.58241 49.56911 50.44976 50.00865 50.88206 54.18908
##        m1.opt.none 40.64339 41.11479 43.28585 41.80860 46.32470 49.90008
##          m1.update 49.44773 50.06004 52.78024 52.35340 54.89347 57.94328
##           m1.start 48.39947 48.81600 50.73314 49.68118 53.13262 54.52046
##          m1.start2 47.48920 48.21813 49.26506 49.02211 49.76483 52.64671
##  m1.start.opt.none 39.76549 39.99929 43.10869 42.90461 45.37880 48.79881
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

![](01-method-update_files/figure-html/plot_tab2-1.png) 

### Model comparison


```r
system.time(mod1 <- lmer(f1, sleepstudy))
```

```
##    user  system elapsed 
##   0.049   0.000   0.049
```

```r
system.time(mod2 <- relmatLmer(f1, sleepstudy, control = lmerControl(optimizer = "none"), start = start.m1))
```

```
##    user  system elapsed 
##    0.04    0.00    0.04
```

```r
system.time(mod3 <- update(mod1, f1, control = lmerControl(optimizer = "none"), start = start.m1))
```

```
##    user  system elapsed 
##   0.041   0.000   0.041
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
