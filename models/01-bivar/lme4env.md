---
title: "Bi-variate mixed models with lme4qtl"
author: "Andrey Ziyatdinov"
date: "2019-01-04"
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

library(lme4)
library(lme4qtl)
library(sommer)
library(solarius)
```

# Data


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


##  Build model by blocks


```r
#formula <- trait1 ~ AGE + SEX + (1|ID)
formula <- trait1 ~ AGE + SEX + (0 + SEX|ID)

data <- dat40
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
List of 1
 $ ID: Factor w/ 224 levels "ID101","ID102",..: 1 2 3 4 5 6 7 8 9 10 ...
 - attr(*, "assign")= int 1
```

```r
fnmns
```

```
[1] "ID"
```

```r
Ztlist <- lmod$reTrms[["Ztlist"]]
```

### Interim between Step 1 and Step 2


```r
i <- 1
fn <- fnmns[i] # "ID"

zn <- lmod$fr[, fn] # ID values: 101 102 103 104
zn.unique <- levels(zn)
obs.unique <- rownames(lmod$fr)

Zt <- Ztlist[[i]]
rownames.Zt <- rownames(Zt) # rownames are (repeated) ID values: ID101 ID102 ID103...
colnames.Zt <- colnames(Zt) # values are meaningless: "1" "2" "3", ..., but correspond to `zn`
nrow.Zt <- nrow(Zt)
ncol.Zt <- ncol(Zt)
stopifnot(ncol.Zt == length(zn))
stopifnot(all(rownames.Zt %in% zn.unique))
stopifnot(all(colnames.Zt %in% obs.unique))

Zt.unique <- Zt[zn.unique, obs.unique]

reln <- rownames(relmat[[fn]])
stopifnot(!is.null(reln))
stopifnot(all(zn.unique %in% reln))

A <- relmat[[fn]][zn.unique, zn.unique]
A <- Matrix::Matrix(A, sparse = TRUE)

R <- relfac(A) # A = R' R 
dim.R <- nrow(R) # nrow(R) = ncol(R)

# new `Zt.unique` matrix, named as `Wt.unique`
Wt.unique <- R %*% Zt.unique # Z* = Z L, then Z*' = L' Z' = R Z'

if(nrow.Zt == dim.R & ncol.Zt == dim.R) {
  Wt <- Wt.unique
} else {
  stopifnot(!(nrow.Zt %% dim.R))
  stopifnot(!(ncol.Zt %% dim.R))
  nrep.row <- nrow.Zt / dim.R
  nrep.col <- ncol.Zt / dim.R
  
  ind.rows <- lapply(zn.unique, function(x) which(rownames.Zt == x)) 
  ind.cols <- lapply(obs.unique, function(x) which(colnames.Zt == x)) 
  stopifnot(all(sapply(ind.rows, length) == nrep.row))
  stopifnot(all(sapply(ind.cols, length) == nrep.col))
  
  Wt <- Zt
  for(ir in seq(1, nrep.row)) {
    for(ic in seq(1, nrep.col)) {
      rows <- sapply(ind.rows, function(x) x[ir])
      cols <- sapply(ind.cols, function(x) x[ic])
      Wt[rows, cols] <- R %*% Wt[rows, cols]
    }
  }
}

# write updated matrices back to model
Ztlist[[i]] <- Wt
lmod$reTrms[["Ztlist"]] <- Ztlist
lmod$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
```

Dummary check:


```r
#A[1:5, 1:5]
#crossprod(R)[1:5, 1:5]
#crossprod(Zt)[1:5, 1:5]
#crossprod(Wt)[1:5, 1:5]
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
REML criterion at convergence: 987.8249
Random effects:
 Groups   Name Std.Dev. Corr
 ID       SEX1 1.454        
          SEX2 2.068    0.67
 Residual      1.468        
Number of obs: 224, groups:  ID, 224
Fixed Effects:
(Intercept)          AGE         SEX2  
    7.47158      0.01122     -0.50887  
```


```r
(mod0 <- relmatLmer(formula, dat40, relmat = list(ID = kin2)))
```

```
Linear mixed model fit by REML ['lmerMod']
Formula: trait1 ~ AGE + SEX + (0 + SEX | ID)
   Data: dat40
REML criterion at convergence: 987.8249
Random effects:
 Groups   Name Std.Dev. Corr
 ID       SEX1 1.454        
          SEX2 2.068    0.67
 Residual      1.468        
Number of obs: 224, groups:  ID, 224
Fixed Effects:
(Intercept)          AGE         SEX2  
    7.47158      0.01122     -0.50887  
```

```r
stopifnot(all.equal(getME(mod, "Ztlist")[[1]]@x, getME(mod0, "Ztlist")[[1]]@x))
stopifnot(all.equal(getME(mod, "Ztlist")[[1]]@Dimnames[[1]], getME(mod0, "Ztlist")[[1]]@Dimnames[[1]]))
stopifnot(all.equal(getME(mod, "Ztlist")[[1]]@Dimnames[[2]], getME(mod0, "Ztlist")[[1]]@Dimnames[[2]]))
stopifnot(all.equal(logLik(mod0), logLik(mod)))
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

# Fit bi-variate model by SOLAR


```r
solarPolygenic(y1 + y2 ~ 1, as.data.frame(dat)) 
```
