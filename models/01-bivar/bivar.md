---
title: "Bi-variate mixed models"
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
h1 <- 0.25
h2 <- 0.75
rg <- 0.5
re <- 0.1

h12 <- rg * sqrt(h1) * sqrt(h2)
e12 <- re * sqrt(1-h2) * sqrt(1-h2)

h12
```

```
[1] 0.2165064
```

```r
e12
```

```
[1] 0.025
```

## Pre-defined kinship matrix


```r
data(dat40)
N <- nrow(dat40)

ids <- dat40$ID
mat <- kin2[ids, ids]
stopifnot(all(rownames(mat) == dat40$ID))

N
```

```
[1] 234
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

    dat <- dplyr::select(dat40, ID) %>%
      mutate(y1 = y1, y2 = y2, rep = r)
    rownames(dat) <- ids
    
    as_data_frame(dat)
  }) %>% bind_rows
}

# Bi-variate model with no repetitions
```

```r
dat <- sim_bivar(h1, h2, rg, re)

m1 <- mmer(cbind(y1, y2) ~ 1, 
  random =~ vs(ID, Gu = as.matrix(mat)), rcov =~ vs(units), dat, 
  verbose = FALSE)
  
m1$sigma_scaled
```

```
$`u:ID`
          y1        y2
y1 0.4433113 0.0876152
y2 0.0876152 0.6944101

$`u:units`
            y1          y2
y1  0.55381652 -0.01596232
y2 -0.01596232  0.31229209
```


```r
data_frame(par = c("h12", "e12"),
    true = c(h12, e12),
    observed = c(m1$sigma_scaled[["u:ID"]][1, 2],
      m1$sigma_scaled[["u:units"]][1, 2])) %>%
  pander
```


-------------------------
par   true     observed  
----- -------- ----------
h12   0.2165   0.08762   

e12   0.025    -0.01596  
-------------------------

# Bi-variate model with (2) repetitions


```r
dat <- sim_bivar(h1, h2, rg, re, R = 2)

m2 <- mmer(cbind(y1, y2) ~ 1, 
  random =~ vs(ID, Gu = as.matrix(mat)), rcov =~ vs(units), dat, 
  verbose = FALSE)
  
m2$sigma_scaled
```

```
$`u:ID`
           y1         y2
y1 0.02100302 0.02764787
y2 0.02764787 0.17277035

$`u:units`
          y1        y2
y1 0.9784089 0.2753681
y2 0.2753681 0.8014491
```
