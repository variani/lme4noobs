### inc
library(gait)
library(solaris)

### par
lmerControl <- lmerControl(check.nobs.vs.nlev = "ignore",
  check.nobs.vs.nRE = "ignore", check.nobs.vs.rankZ = "ignore")
  
### data
pdat  <- gait1.phen()
pdat <- mutate(pdat, 
  tr_FXI = FXI_T * 5.1)

k2 <- solarKinship2(pdat)

#ids <- with(pdat, ID[!is.na(AGE)])
#pdat <- subset(pdat, ID %in% ids)
#ind <- which(rownames(k2) %in% ids)
#k2 <- k2[ind, ind]

### polygenic
mod <- solaris(tr_FXI ~ AGE + (1|HHID) + (1|ID), pdat, 
  relmat = list(ID = k2), control = lmerControl)

### assoc function
createFormula <- function(resp, fixed, rand) {
  f <- reformulate(c(fixed,rand), response = resp)
  ## use the parent (createModel) environment, not the
  ## environment of this function (which does not contain 'data')
  environment(f) <- parent.frame()
  f
}

run_assoc <- function(dat, genocov.files)
{
  B <- length(genocov.files)
  
  out <- llply(1:B, function(i) {
    genocov.file <- genocov.files[i]
    
    sdat <- read.table(genocov.file, sep = ",", header = TRUE)
    sdat <- rename(sdat, c(id = "ID"))
    snp.names <- grep("^snp_", names(sdat), value = TRUE)
    
    dat2 <- merge(dat, sdat, by = "ID")

    f <- createFormula("tr_FXI", c("AGE"), c("(1|HHID)", "(1|ID)")) 
    mod2 <- solaris(f, dat2, 
      relmat = list(ID = k2), control = lmerControl)
    
        
    f.snps <- as.formula(paste("~", paste(snp.names, collapse = "+")))
    
    list(batch.index = i, add1 = add1(mod2, f.snps, test = "Chisq"))
  })   
  
  return(out) 
}

### snp files
snp.files <- gait1.snpfiles(num.snpdirs = 2)
genocov.files <- snp.files$genocov.files

pdat2 <- subset(pdat, select = c("tr_FXI", "AGE", "HHID", "ID"))
out <- run_assoc(pdat2, genocov.files)
