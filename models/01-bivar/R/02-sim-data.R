# notes:
# - simulate traits & scale them, so h2 is equal to variances

### inc
library(dplyr)
library(magrittr)

library(Matrix)
library(sparseMVN)
library(MASS)

library(lme4qtl)
library(sommer)
library(solarius)

### par
h1 <- 0.25
h2 <- 0.75
rg <- 0.5
re <- 0.1

### data
data(dat40)
N <- nrow(dat40)

ids <- dat40$ID
mat <- kin2[ids, ids]
stopifnot(all(rownames(mat) == dat40$ID))

### sim.
V1 <- h1 * mat + (1 - h1) * Diagonal(N)
rownames(V1) <- ids
colnames(V1) <- ids

#ch <- Cholesky(V, prec = FALSE)
#y <- rmvn.sparse(1, rep(0, N), ch, FALSE) %>% as.numeric
y1 <- mvrnorm(1, rep(0, N), V1) %>% scale %>% as.numeric

V2 <- h2 * mat + (1 - h2) * Diagonal(N)
rownames(V2) <- ids
colnames(V2) <- ids

y2 <- mvrnorm(1, rep(0, N), V2) %>% scale %>% as.numeric

dat <- mutate(dat40, y1 = y1, y2 = y2)
rownames(dat) <- ids

p1 <- solarPolygenic(y1 ~ 1, dat)
p2 <- solarPolygenic(y2 ~ 1, dat)

m1 <- mmer(y1 ~ 1, ~vs(ID, Gu = as.matrix(mat)), ~vs(units), data = dat)
m2 <- mmer(y2 ~ 1, ~vs(ID, Gu = as.matrix(mat)), ~vs(units), data = dat)

### bivariate
h12 <- rg*sqrt(h1)*sqrt(h2)
Sigma_gen <- matrix(c(h1, h12, h12, h2), 2, 2)
V_gen <- kronecker(Sigma_gen, mat)

e12 <- re*sqrt(1-h2)*sqrt(1-h2)
Sigma_resid <- matrix(c(1-h1, e12, e12, 1-h2), 2, 2)
V_resid <- kronecker(Sigma_resid, Diagonal(N))

V <- V_gen + V_resid

y12 <- mvrnorm(1, rep(0, 2*N), V)
y1 <- y12[seq(1, N)] %>% scale %>% as.numeric
y2 <- y12[N + seq(1, N)] %>% scale %>% as.numeric

bdat <- mutate(dat40, y1 = y1, y2 = y2)
rownames(bdat) <- ids

bp <- solarPolygenic(y1 + y2 ~ 1, bdat)
bm <- mmer(cbind(y1, y2) ~ 1, random =~ vs(ID, Gu = as.matrix(mat)), rcov =~ vs(units), bdat)

