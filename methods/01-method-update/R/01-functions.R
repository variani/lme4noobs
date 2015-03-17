library(lme4)

createFormula <- function(resp, fixed, rand) {
  f <- reformulate(c(fixed,rand), response = resp)
  ## use the parent (createModel) environment, not the
  ## environment of this function (which does not contain 'data')
  environment(f) <- parent.frame()
  f
}

createModel <- function(data) {
  mf.final <- createFormula("Reaction", "Days", "(Days|Subject)")
  lmer(mf.final, data=data)
}

#mod <- createModel(data=sleepstudy)
#drop1(mod)

run_assoc <- function(dat, snps)
{
  ### option 1
  #mod <- lmer(Reaction ~ (1 | Subject), dat)
  ### option 2
  sdat <- snps[[1]]
  N <- ncol(sdat)
  cnames <- paste0("X", 1:N)
  colnames(sdat) <- cnames
  dat2 <- cbind(dat, sdat)
  mod2 <- lmer(Reaction ~ (1 | Subject), dat2)
  
  out <- llply(1:length(snps), function(i) {
    sdat <- snps[[i]]

    N <- ncol(sdat)
    cnames <- paste0("X", 1:N)
    colnames(sdat) <- cnames
    dat2[,  ncol(dat) + seq(1, ncol(sdat))] <- sdat
    
    f <- as.formula(paste("~", paste(cnames, collapse = "+")))
    mod2 <- lmer(Reaction ~ (1 | Subject), dat2)
    
    list(batch.index = i, add1 = add1(mod2, f, test = "Chisq"))
  })   
  
  return(out) 
}


N <- 2
set.seed(1)
dat1 <- matrix(rnorm(nrow(sleepstudy)*N), nrow(sleepstudy), N)
dat2 <- matrix(rnorm(nrow(sleepstudy)*N), nrow(sleepstudy), N)
colnames(dat1) <- paste0("X", 1:N)
colnames(dat2) <- paste0("X", 1:N)

out <- run_assoc(sleepstudy, list(dat1, dat2))
