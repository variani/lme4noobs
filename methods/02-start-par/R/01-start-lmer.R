### inc
library(lme4)

### models
m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, verbose = 2)

th <- getME(m1, "theta")
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, start = th, verbose = 2)

control <- lmerControl(optCtrl = list(maxfun = 1, maxit = 1))
m3 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, start = th, control = control, verbose = 2)
          
