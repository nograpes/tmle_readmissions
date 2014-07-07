# /usr/bin/R --args data_dump/disease_ami.object data_dump/glmnet_Q_model_ami.object 
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(survival))
registerDoMC(cores=10) # For 10 folds

# arguments<-c('data_dump/disease_ami.object', 'data_dump/glmnet_Q_model_ami.object')
arguments <- commandArgs(trailingOnly=TRUE)

object.file<-arguments[1]
output.file<-arguments[2]
load(object.file)
set.seed(1)


system.time(glmnet.predict.outcome<-cv.glmnet(x=disease.big.matrix[1:1000,], 
                                              y=Surv(disease.df$tte, !disease.df$censor)[1:1000,],
                                              family='cox',
                                              parallel=TRUE))



system.time(glmnet.predict.outcome <- 
              glmnet(x=disease.big.matrix, 
                     y=Surv(disease.df$tte, 
                           !disease.df$censor),
                     family='cox' ))

system.time(glmnet.predict.censor <- 
              glmnet(x=disease.big.matrix, 
                     y=Surv(disease.df$tte, 
                            disease.df$censor),
                     family='cox' ))
# Convergence after 37th lambda not reached. (Too small?)
system.time(glmnet.predict.exposure <- 
              glmnet(x=disease.big.matrix, 
                     y=disease.df$hosp,
                     family='multinomial' ))



# I hacked the survival package and scraped out the 
# survival function estimator, so I can use it with a glmnet
# found Cox model.

get.baseline.surv <- function(y, x, coef) {
  # Clean up x
  x <- as.matrix(x)
  good.indices <- complete.cases(x)
  x <- x[good.indices,,drop=FALSE]
  y <- y[good.indices,,drop=FALSE]
  xcenter <- colMeans(x)
  scaled.x <- scale(x, center = xcenter, scale = FALSE)
  scaled.risk <- c(exp(scaled.x %*% coef))
  list(center=xcenter,coef=coef, fit=survival:::agsurv(y=y, x=scaled.x, 
                    wt=rep(1,nrow(x)), 
                    risk=scaled.risk, 
                    survtype=3, vartype=3))
}

get.surv <- function(y, x, coef, new.x) {
  expand <- function(fit, newrisk) {
    surv <- exp(-fit$cumhaz)
    fit$surv <- surv^newrisk
    fit$cumhaz <- fit$cumhaz * newrisk
    fit
  }

  center.and.fit <- get.baseline.surv(y,x,coef)
  fit <- center.and.fit$fit
  center <- center.and.fit$center
  scaled.new.x <- scale(new.x, center = xcenter, scale = FALSE)
  scaled.new.risk <- c(exp(scaled.new.x %*% coef))
  expand(fit, scaled.new.risk)$surv
}


surv.precomputed.baseline <- function(baseline, new.x) {
  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
  base <- exp(-baseline$fit$cumhaz)
  lapply(c(scaled.new.risk),function(x) base^x)
}


x <- disease.big.matrix
y <- Surv(disease.df$tte, !disease.df$censor)
glmnet.model <- glmnet.predict.outcome
s = 12

betas <- glmnet.model$beta[,s]
betas <- betas[betas!=0]
cut.x <- x[,names(betas)]


debugonce(get.baseline.surv)

baseline <- get.baseline.surv(y, cut.x, coef=betas)
surv.precomputed.baseline(baseline, new.x=cut.x)


get.surv(y=y, x=cut.x, coef=betas, new.x=cut.x[1,,drop=FALSE])




# Estimate hazard

# Estimate cumulative hazard

# surv = exp(-cumhaz)

# surv^newrisk


which(glmnet.predict.censor$beta[,36]!=0)




system.time(glmnet.predict.censor <- 
              glmnet(x=disease.big.matrix, 
                     y=,
                     family='cox' ))

apply(glmnet.predict.outcome$beta[,1:22],2,function(x) names(which(x!=0)))


apply(glmnet.predict.outcome$beta[,1:42],2,function(x) grep('hosp',names(which(x!=0)),value=TRUE))


?grep


glmnet.predict.outcome<-cv.glmnet(x=disease.big.matrix, 
                                  y=Surv(disease.df$tte, !disease.df$censor),
                                  family='cox',
                                  parallel=TRUE)



str(glmnet.predict.outcome)

save(glmnet.predict.outcome,file=output.file)