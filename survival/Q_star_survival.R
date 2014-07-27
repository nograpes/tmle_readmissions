suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(bigrf))
registerDoMC(cores=10) # For 10 folds
options(mc.cores=12)

arguments <- commandArgs(trailingOnly=TRUE)
# For testing purposes.
if (interactive()) arguments<-c('data_dump/disease_heart_failure.object',
                                'data_dump/rf_G_calibrated_model_heart_failure.object',
                                'survival/data_dump/glmnet_g_censor_heart_failure.object',
                                'survival/data_dump/glmnet_Q_heart_failure.object',
                                'survival/Q_star_survival.R',
                                'survival/data_dump/Q_star_survival_heart_failure.object')

object.file <- arguments[1]
rf.G.object <- arguments[2]
survival.censor.object <- arguments[3]
survival.Q.object <- arguments[4]
R.file <- arguments[5] # Because I can't get the Makefile to not pass these.
output.file <- arguments[6]

load(object.file)
load(rf.G.object)
load(survival.Q.object) # glmnet.predict.outcome
load(survival.censor.object) # glmnet.predict.censor

# I hacked the survival package and scraped out the
# survival function estimator, so I can use it with a glmnet
# found Cox model.
get.betas <- function(glmnet.model, s) {
  betas <- glmnet.model$glmnet.fit$beta[,s]
  betas[betas!=0]
}

get.cut.x <- function(x, betas)
  x[,names(betas)]


# Actually, it uses Fleming and Harrington, not Kalbfleisch and Prentice.
get.kalb.prent.surv <- function(y, x, coef) {
  # Clean up x
  x <- as.matrix(x)
  good.indices <- complete.cases(x)
  x <- x[good.indices,,drop=FALSE]
  y <- y[good.indices,,drop=FALSE]
  xcenter <- colMeans(x)
  scaled.x <- scale(x, center = xcenter, scale = FALSE)
  scaled.risk <- c(exp(scaled.x %*% coef))
  fit <- survival:::agsurv(y=y, x=scaled.x,
                    wt=rep(1,nrow(x)),
                    risk=scaled.risk,
                    survtype=3, vartype=3)
  cumhaz <- fit$cumhaz
  unique.times <- fit$time
  cumhaz <- cumhaz[findInterval(seq.int(max(unique.times)), unique.times)]
  list(center=xcenter, max=max(unique.times), coef=coef, surv=exp(-cumhaz)) # Strip the first cumhaz, it is zero.
}

get.surv.at.new.x <- function (baseline, new.x, times) {
  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
  base <- baseline$surv
  surv <- mapply(function(x,y) (base^(x))[1:y], c(scaled.new.risk), times)
  # By definition, you can't have an estimate at the very end of the
  # longest time. So, you'll get NAs for those. Reset those to the last
  # estimated survival time.
  surv <- lapply(surv, function(x) {
    if (length(x)>=baseline$max) {
      last.not.NA <- x[max(which(!is.na(x)))]
	  x[baseline$max:length(x)] <- last.not.NA
    }
	x
  })
  tte <- mapply(function(x,y) sum(base^(x)), c(scaled.new.risk), times)
  list(surv = surv, tte = tte)
}

# I thought I would need this, but I don't.
# I really understand the hazard now, though.
get.hazard.at.new.x <- function (surv)
  lapply(surv, function(x) diff(-log(x)))

inv.cumprod <- function(x) x[-1] / x[-length(x)]

get.conditional.failure.at.new.x <- function(surv)
  lapply(surv, function(s) 1 - inv.cumprod(s))

# The other choice for s is lambda.1se, where you select the biggest
# lambda that is within 1 standard error of lambda min.
drop.all.probs <- function(glmnet.model, y, x, s.type='lambda.min') {
  s <- match(glmnet.model[[s.type]], glmnet.model$lambda)
  betas <- get.betas(glmnet.model, s)
  cut.x <- get.cut.x(x, betas)
  baseline <- get.kalb.prent.surv(y, cut.x, coef=betas)

  hosps <- levels(disease.df$hosp)
  cols <- paste0('hosp', hosps)

  x.mat.list <- mclapply(cols, function(x) {
    if (x %in% colnames(cut.x)) cut.x[,x]<-1
    cut.x
  })

  times <- disease.df$tte
  surv.and.tte <- mclapply(x.mat.list, get.surv.at.new.x,
                           baseline=baseline, times=times+1) # Because first is always 1
  surv <- lapply(surv.and.tte, "[[", 'surv')
  tte <- lapply(surv.and.tte, "[[", 'tte')

  # Also need to calculate the conditional failure for the observed.
  surv.as.observed <- get.surv.at.new.x(cut.x, baseline=baseline, times=times+1)$surv
  cond.fail.as.observed <- get.conditional.failure.at.new.x(surv.as.observed)

  hazard <- mclapply(surv, get.hazard.at.new.x)
  cond.fail <- mclapply(surv, get.conditional.failure.at.new.x)
  surv <- lapply(surv, function(x) lapply(x, function(y) y[-1])) # Now clip the first surv prob.
  names(tte) <- names(cond.fail) <- names(hazard) <- names(surv) <- hosps
  list(surv=surv, cond.fail=cond.fail, hazard=hazard, 
       #s.ratio=s.ratio, 
       tte=tte, cond.fail.as.observed=cond.fail.as.observed)
}

# S0 and Q.
# debugonce(drop.all.probs)
outcome.S0.and.Q.list <- drop.all.probs(
  glmnet.model=glmnet.predict.outcome,
  y=Surv(disease.df$tte, !disease.df$censor),
  x=disease.big.matrix)
hosps=names(outcome.S0.and.Q.list$surv)

# Calculate g2 (probability of censoring by time)
censor.probs <-drop.all.probs(
  glmnet.model=glmnet.predict.censor,
  y=Surv(disease.df$tte, disease.df$censor),
  x=disease.big.matrix)$surv
names(censor.probs) <- hosps

# Calculate g1 (probability of exposure)
# G model - calibrated RF
g.votes <- rf.predict.exposure@oobvotes
g.by.rf.unscaled <- g.votes / rowSums(g.votes)
one.v.all = sapply(colnames(g.by.rf.unscaled), function(x) disease.df$hosp == x)


# Important to scale each column of g
g.by.rf <- mapply(function(y,x) predict(glm(y~x, family=binomial),
                                        type='response'),
                  data.frame(one.v.all, check.names=FALSE),
                  data.frame(g.by.rf.unscaled))

prob.of.exposure.to.exposed <- g.by.rf[cbind(
  seq.int(nrow(disease.df)),
  match(disease.df$hosp,colnames(g.by.rf)))]

# Reorganize, so each hospital has its own data.frame
# with Q, S.ratio, S0, g1, and g2.
hosp.dfs <- lapply(hosps,
			   function(hosp)
				data.frame(
				  # This is the same for all hospitals.
            id      = rep(seq.int(nrow(disease.df)), times=disease.df$tte),
            A       = rep(disease.df$hosp, times=disease.df$tte),
            QAW     = unlist(outcome.S0.and.Q.list$cond.fail.as.observed),
            Q.star  = unlist(outcome.S0.and.Q.list$cond.fail.as.observed),
				    Q       = unlist(outcome.S0.and.Q.list$cond.fail[[hosp]]),
						hazard  = unlist(outcome.S0.and.Q.list$hazard[[hosp]]),
						# S.ratio = unlist(outcome.S0.and.Q.list$s.ratio[[hosp]]),
						S0      = unlist(outcome.S0.and.Q.list$surv[[hosp]]),
# 						g1      = rep(ifelse(disease.df$hosp==hosp,
# 											 prob.of.exposure.to.exposed,
# 											 0), times=disease.df$tte),
            g1      = rep(prob.of.exposure.to.exposed, times=disease.df$tte),
						g2      = unlist(censor.probs[[hosp]]),
						time    = unlist(lapply(disease.df$tte,seq.int)),
						tte     = rep(disease.df$tte, disease.df$tte),
						censor  = rep(disease.df$censor, disease.df$tte)
				))
names(hosp.dfs) <- hosps


# This is a practical consideration to prevent positivity violations.
# There will be very little data at the end of the survival curve.
# By clipping the 1% quantile, we can remove these violations, and prevent
# having enormous weights (1600) over a period of 1000 days, which overwhelms
# anything else in the model.
# It is also necessary to know the max time because it is t_k, and the
# S.ratio depends on t_k.
# But, to be clear, for heart failure this ends up being five years.
max.time <- quantile(unlist(lapply(disease.df$tte,seq.int)),probs=0.975)

# I must now find S(t_k) for everybody.
glmnet.model <- glmnet.predict.outcome
x <- disease.big.matrix
y <- Surv(disease.df$tte, !disease.df$censor)
s.type='lambda.min'
s <- match(glmnet.model[[s.type]], glmnet.model$lambda)
betas <- get.betas(glmnet.model, s)
cut.x <- get.cut.x(x, betas)
baseline <- get.kalb.prent.surv(y, cut.x, coef=betas)
baseline.t_k <- baseline$surv[max.time]
scaled.new.x <- scale(cut.x, center = baseline$center, scale = FALSE)
scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
S.t_k <- baseline.t_k^c(scaled.new.risk)

pinch.low <- function(x, low=0.025)
  ifelse(x < low, low, x)

clip <- function(x){
  second.lowest <- min(x[x!=min(x)])
  x[x==min(x)] <- second.lowest
  x
}

for(i in seq_along(hosp.dfs)) {
  hosp.dfs[[i]]$Q.star <- clip(hosp.dfs[[i]]$Q.star)
  hosp.dfs[[i]]$n1 <- with(hosp.dfs[[i]], !censor & (time==tte))
}

# Calculate S.ratio 
calc.S <- function(hosp.df) {
  split.Q.star <- split(hosp.df$Q.star, hosp.df$id)
  unlist(lapply(split.Q.star, function(x) cumprod(1-x)))
}

# Calculate clever covariate.
calc.clever.covariate <- function(hosp.df) {
  hosp.df$S0 <- calc.S(hosp.df)
  hosp.df$clever.covariate <- hosp.df$fixed.part / hosp.df$S
  # A simple fix
  hosp.df$clever.covariate <- ifelse(hosp.df$clever.covariate > 100, 100, hosp.df$clever.covariate)
  hosp.df
}

# Seriously, I don't need the whole thing.
# Those zeroes will always just add a constant to the likelihood
# in a single variable function.
reduced.hosp.dfs <- mapply(function(hosp.df,hosp) hosp.df[hosp==hosp.df$A,],
                           hosp.dfs, names(hosp.dfs), SIMPLIFY=FALSE)

for(i in seq_along(reduced.hosp.dfs)) {
  hosp <- names(reduced.hosp.dfs)[i]
  S.t_k.hosp <- S.t_k[disease.df$hosp==hosp]
  tte <- disease.df$tte [disease.df$hosp==hosp]
  S.t_k.hosp <- rep(S.t_k.hosp, tte)
  reduced.hosp.dfs[[i]]$fixed.part <- S.t_k.hosp / 
          (reduced.hosp.dfs[[i]]$g1*reduced.hosp.dfs[[i]]$g2)
  cutoff <- quantile(hosp.df$fixed.part,0.975)
  reduced.hosp.dfs[[i]]$fixed.part <-
    ifelse(reduced.hosp.dfs[[i]]$fixed.part > cutoff,
           cutoff, reduced.hosp.dfs[[i]]$fixed.part)
}


fit.epsilon <- function(hosp.df) 
  glm(n1~clever.covariate-1,offset=qlogis(Q.star), data=hosp.df)$coef

update.Q.star <- function(hosp.df,epsilon) {
  hosp.df$Q.star <- with(hosp.df,
                         plogis(qlogis(Q.star) + (epsilon*clever.covariate)))
  hosp.df
}

is.bad<-function(x) is.na(x) | is.infinite(x) | is.nan(x)

iteration <- 1
delta <- 1e-10
while (iteration==1 || any(epsilons>delta)) {
  reduced.hosp.dfs <- lapply(reduced.hosp.dfs, calc.clever.covariate)
  epsilons <- sapply(reduced.hosp.dfs, fit.epsilon) 
  reduced.hosp.dfs <- mapply(update.Q.star,reduced.hosp.dfs,epsilons,
                             SIMPLIFY=FALSE)
  if(iteration!=1) epsilon.mat <- rbind(epsilon.mat,epsilons)
  else epsilon.mat <- t(epsilons)
  iteration <- iteration + 1
}

epsilons <- colSums(epsilon.mat)


# Functions to extract the tte based on the new epsilon
compute.tte.by.baseline.and.epsilon <- function(baseline, x, epsilon) {
  scaled.x <- scale(x, center = baseline$center, scale = FALSE)
  scaled.risk <- exp(scaled.x %*% baseline$coef)
  base <- baseline$surv
  surv <- lapply(scaled.risk, function(x) base^x)
  cond.fail <- get.conditional.failure.at.new.x(surv)
  cond.fail.plus.epsilon <- lapply(cond.fail,function(x)
                                   plogis(qlogis(x) + epsilon))
  sapply(cond.fail.plus.epsilon, function(x) sum(cumprod(1-x)))
}

compute.tte.by.epsilon.and.hosp <- function(glmnet.model, hosp,
                                   s.type='lambda.min', epsilon) {
  y<-disease.df$outcome
  x<-disease.big.matrix
  s <- match(glmnet.model[[s.type]], glmnet.model$lambda)
  betas <- get.betas(glmnet.model, s)
  cut.x <- get.cut.x(x, betas)
  
  if(hosp %in% colnames(cut.x)) cut.x[,hosp] <- 1
  
  baseline <- get.kalb.prent.surv(y, cut.x, coef=betas)
  compute.tte.by.baseline.and.epsilon(baseline, cut.x, epsilon)
}

debugonce(compute.tte.by.baseline.and.epsilon)

computed.tte <- mapply(compute.tte.by.epsilon.and.hosp, 
                       hosp=hosps,epsilon=epsilons,
                       MoreArgs=list(glmnet.model=glmnet.predict.outcome))

tte.mat <- computed.tte
# A couple of renames before the save
save(tte.mat, epsilon.mat, file=output.file)
