### inc
library(dplyr)
library(magrittr)

library(lme4qtl)
library(sommer)
library(solarius)

### data
data(dat40)

mat <- as.matrix(kin2)
#ids <- dat40$ID
#mat <- mat[ids, ids]
#stopifnot(all(rownames(mat) == dat40$ID))

dat <- within(dat40, {
  trait1 <- as.numeric(scale(trait1)) 
  trait2 <- as.numeric(scale(trait2))
})

### fit by sommer   
#m0 <- mmer(cbind(trait1, trait2) ~ 1, random =~ vs(ID), rcov =~ vs(units), dat)
m1 <- mmer(cbind(trait1, trait2) ~ 1, random =~ vs(ID, Gu = mat), rcov =~ vs(units), dat)

### fit by solarius
m2 <- solarPolygenic(trait1 + trait2 ~ 1, dat)

### print
m1$sigma_scaled

m2$vcf

