---
title: "Bi-variate mixed models with lme4qtl"
author: "Andrey Ziyatdinov"
date: "2019-01-03"
output:
  html_document:
    theme: united
    toc: true
    number_sections: true
    keep_md: true
    includes:
      in_header: header.html
---






```r
library(Matrix)
library(MASS)

library(lme4qtl)
library(sommer)
library(solarius)
```

# Data simulations

## Parameters


```r
h1 <- 0.5
h2 <- 0.5
rg <- 0.5
re <- 0

h12 <- rg * sqrt(h1) * sqrt(h2)
e12 <- re * sqrt(1-h2) * sqrt(1-h2)

h12
```

```
[1] 0.25
```

```r
e12
```

```
[1] 0
```

## Pre-defined kinship matrix


```r
data(dat40)
N <- nrow(dat40)
N
```

```
[1] 234
```

```r
# new ids
ids <- dat40$ID

ids <- paste0("ID", ids)
dat40$ID <- ids
rownames(kin2) <- colnames(kin2) <- ids

mat <- kin2[ids, ids]
stopifnot(all(rownames(mat) == dat40$ID))

dat40 <- within(dat40, RID <- ID)
```

## Function to simulate bi-variate


```r
sim_bivar <- function(h1, h2, rg, re, R = 1, scale = TRUE) 
{
  stopifnot(scale)
  
  h12 <- rg*sqrt(h1)*sqrt(h2)
  Sigma_gen <- matrix(c(h1, h12, h12, h2), 2, 2)
  V_gen <- kronecker(Sigma_gen, mat)

  e12 <- re*sqrt(1-h2)*sqrt(1-h2)
  Sigma_resid <- matrix(c(1-h1, e12, e12, 1-h2), 2, 2)
  V_resid <- kronecker(Sigma_resid, Diagonal(N))

  V <- V_gen + V_resid
  
  lapply(seq(1, R), function(r) {
    y12 <- mvrnorm(1, rep(0, 2*N), V)
    y1 <- y12[seq(1, N)] %>% scale %>% as.numeric
    y2 <- y12[N + seq(1, N)] %>% scale %>% as.numeric

    dat <-  dat40 %>% # dplyr::select(dat40, ID, RID)
      mutate(y1 = y1, y2 = y2, rep = r)
    rownames(dat) <- ids
    
    as_data_frame(dat)
  }) %>% bind_rows
}
```

# Bi-variate model with no repetitions


```r
dat <- sim_bivar(h1, h2, rg, re)

bdat <- gather(dat, trait, y, y1, y2)
```

##  Build model by blocks


```r
#m <- relmatLmer(y ~ (0 + trait|ID) + (0 + dummy(trait)|RID), bdat, relmat = list(ID = kin2))
formula <- y ~ (0 + trait|ID) + (0 + dummy(trait)|RID)
data <- bdat
relmat <- list(ID = kin2)

control <- lmerControl()
control$checkControl$check.nobs.vs.rankZ <- "ignore"
control$checkControl$check.nobs.vs.nlev <- "ignore"
control$checkControl$check.nobs.vs.nRE <- "ignore"
```

### Step 1


```r
lmod <- lFormula(formula, data, control = control)
```

Parse the result of Step 1:


```r
flist <- lmod$reTrms[["flist"]]
fnmns <- names(flist) 

str(flist)
```

```
List of 2
 $ ID : Factor w/ 234 levels "ID101","ID102",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ RID: Factor w/ 234 levels "ID101","ID102",..: 1 2 3 4 5 6 7 8 9 10 ...
 - attr(*, "assign")= int [1:2] 1 2
```

```r
fnmns
```

```
[1] "ID"  "RID"
```

```r
Ztlist <- lmod$reTrms[["Ztlist"]]
```

### Interim between Step 1 and Step 2


```r
i <- 1
fn <- fnmns[i] # "ID"

zn <- lmod$fr[, fn] # ID values: 101 102 103 104
zn.unique <- unique(zn) 

Zt <- Ztlist[[i]]
zlvl <- rownames(Zt) # rownames are (repeated) ID values: ID101 ID102 ID103...
zobs <- colnames(Zt) # values are meaningless: "1" "2" "3", ..., but correspond to `zn`
stopifnot(length(zobs) == length(zn))

rows <- sapply(zn.unique, function(x) which(zlvl == x)%>% head(1))
cols <- sapply(zn.unique, function(x) which(zn == x) %>% head(1))
Zt.unique <- Zt[rows, cols]

reln <- rownames(relmat[[fn]])
stopifnot(!is.null(reln))
stopifnot(all(zn.unique %in% reln))

A <- relmat[[fn]][zn.unique, zn.unique]
A <- Matrix::Matrix(A, sparse = TRUE)
R <- relfac(A, "evd") # A = R' R 

# new `Zt.unique` matrix, named as `Wt.unique`
Wt.unique <- R %*% Zt.unique # Z* = Z L, then Z*' = L' Z' = R Z'

rows <- sapply(zlvl, function(x) which(zn.unique == x))
cols <- sapply(zn, function(x) which(zn.unique == x)) # zobs or colnames(Zt) corresponds to `zn`
Wt <- Wt.unique[rows, cols]

# write updated matrices back to model
Ztlist[[i]] <- Wt
lmod$reTrms[["Ztlist"]] <- Ztlist
```

Dummary check:


```r
A[1:5, 1:5]
```

```
5 x 5 sparse Matrix of class "dsCMatrix"
      ID101 ID102 ID103 ID104 ID105
ID101   1.0   .     0.5   0.5   0.5
ID102   .     1.0   0.5   0.5   0.5
ID103   0.5   0.5   1.0   0.5   0.5
ID104   0.5   0.5   0.5   1.0   0.5
ID105   0.5   0.5   0.5   0.5   1.0
```

```r
crossprod(R)[1:5, 1:5]
```

```
5 x 5 sparse Matrix of class "dsCMatrix"
      ID101 ID102 ID103 ID104 ID105
ID101   1.0   0.0   0.5   0.5   0.5
ID102   0.0   1.0   0.5   0.5   0.5
ID103   0.5   0.5   1.0   0.5   0.5
ID104   0.5   0.5   0.5   1.0   0.5
ID105   0.5   0.5   0.5   0.5   1.0
```

```r
crossprod(Zt.unique)[1:5, 1:5]
```

```
5 x 5 sparse Matrix of class "dsCMatrix"
  1 2 3 4 5
1 1 . . . .
2 . 1 . . .
3 . . 1 . .
4 . . . 1 .
5 . . . . 1
```

```r
crossprod(Wt.unique)[1:5, 1:5]
```

```
5 x 5 sparse Matrix of class "dsCMatrix"
    1   2   3   4   5
1 1.0 0.0 0.5 0.5 0.5
2 0.0 1.0 0.5 0.5 0.5
3 0.5 0.5 1.0 0.5 0.5
4 0.5 0.5 0.5 1.0 0.5
5 0.5 0.5 0.5 0.5 1.0
```

### Step 2


```r
devfun <- do.call(mkLmerDevfun, c(lmod,
    list(control = control)))
```

### Step 3


```r
opt <- optimizeLmer(devfun,
  optimizer = control$optimizer,
  restart_edge = control$restart_edge,
  boundary.tol = control$boundary.tol,
  control = control$optCtrl,
  calc.derivs = control$calc.derivs,
  use.last.params = control$use.last.params)
```

### Step 4


```r
mod <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
```

### Show the model


```r
mod
```

```
Linear mixed model fit by REML ['lmerMod']
REML criterion at convergence: 1320.6
Random effects:
 Groups   Name         Std.Dev. Corr
 ID       traity1      0.7734       
          traity2      0.5794   0.45
 RID      dummy(trait) 0.5123       
 Residual              0.6326       
Number of obs: 468, groups:  ID, 234; RID, 234
Fixed Effects:
(Intercept)  
   4.86e-17  
```

### Compare true vs. observed


```r
vf <- VarCorr(mod) %>% as.data.frame

h1.obs <- vf[1, "vcov"]
h2.obs <- vf[2, "vcov"]
e1.obs <- vf[4, "vcov"] + vf[5, "vcov"]
e2.obs <- vf[5, "vcov"]

rg.obs <- vf[3, "sdcor"]

data_frame(par = c("h1", "h2", "rg"),
    true = c(h1, h2, rg),
    observed = c(h1.obs/(h1.obs+e1.obs), 
      h2.obs/(h2.obs+e2.obs), rg.obs)
  ) %>%
  pander
```


-----------------------
par   true   observed  
----- ------ ----------
h1    0.5    0.4745    

h2    0.5    0.4562    

rg    0.5    0.4524    
-----------------------

# Fit bi-variate model by SOLAR


```r
solarPolygenic(y1 + y2 ~ 1, as.data.frame(dat)) 
```

```

Call: solarPolygenic(formula = y1 + y2 ~ 1, data = as.data.frame(dat))

File polygenic.out:
	Pedigree:    dat.ped 
	Phenotypes:  dat.phe 
	Trait:       y1 y2                 Individuals:  234 
 
			 H2r(y1) is 0.5656683   
	       H2r(y1) Std. Error:  0.1626428 
 
			 H2r(y2) is 0.1604673   
	       H2r(y2) Std. Error:  0.1518465 
 
			 RhoE is 0.0441958   
	       RhoE Std. Error:  0.1854189 
 
			 RhoG is 0.5828165 
	       RhoG Std. Error:  0.3626984 
 
	       Derived Estimate of RhoP is 0.2022802 
 
 
	Loglikelihoods and chi's are in y1.y2/polygenic.logs.out 
	Best model is named poly and null0 
	Final models are named poly, spor 
```
