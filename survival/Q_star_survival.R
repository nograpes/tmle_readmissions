suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(survival))
registerDoMC(cores=10) # For 10 folds
options(mc.cores=12)

arguments <- commandArgs(trailingOnly=TRUE)
# For testing purposes.
if (interactive()) arguments<-c('data_dump/disease_ami.object', 'data_dump/rf_G_model_ami.object', 'survival/data_dump/glmnet_g_censor_ami.object', 'survival/data_dump/glmnet_Q_ami.object', 'survival/Q_star_survival.R', 'survival/data_dump/Q_star_ami.object')

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
  # Later on, when I need the time-to-event, I'll need the "long" form as well.
  # Wastes memory, but will be super-convenient.
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

S0.to.S.ratio <- function(S0)
  S0[length(S0)]/S0
 
# The other choice for s is lambda.1se, where you select the biggest
# lambda that is within 1 standard error of lambda min.
drop.all.probs <- function(glmnet.model, y, x, s.type='lambda.min') {
  s <- match(glmnet.model[[s.type]], glmnet.model$lambda)
  betas <- get.betas(glmnet.model, s)
  cut.x <- get.cut.x(x, betas)
  baseline <- get.kalb.prent.surv(y, cut.x, coef=betas)

  all.zero <- cut.x
  hosps <- levels(disease.df$hosp)
  cols <- paste0('hosp', hosps)
  
  x.mat.list <- mclapply(cols, function(x) {
    if (x %in% colnames(all.zero)) all.zero[,x]<-1
    all.zero
  })
  
  times <- disease.df$tte
  surv.and.tte <- mclapply(x.mat.list, get.surv.at.new.x, 
                           baseline=baseline, times=times+1) # Because first is always 1
  surv <- lapply(surv.and.tte, "[[", 'surv')
  tte <- lapply(surv.and.tte, "[[", 'tte')
  hazard <- mclapply(surv, get.hazard.at.new.x)
  cond.fail <- mclapply(surv, get.conditional.failure.at.new.x)
  surv <- lapply(surv, function(x) lapply(x, function(y) y[-1])) # Now clip the first surv prob.
  s.ratio <- mclapply(surv, function(x) lapply(x,S0.to.S.ratio))
  names(tte) <- names(cond.fail) <- names(s.ratio) <- names(hazard) <- names(surv) <- hosps
  list(surv=surv, cond.fail=cond.fail, hazard=hazard, s.ratio=s.ratio, tte=tte)
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
				        id      = rep(seq.int(nrow(disease.df)), times=disease.df$tte),
				        Q    = unlist(outcome.S0.and.Q.list$cond.fail[[hosp]]),
						hazard  = unlist(outcome.S0.and.Q.list$hazard[[hosp]]),
						S.ratio = unlist(outcome.S0.and.Q.list$s.ratio[[hosp]]),
						S0      = unlist(outcome.S0.and.Q.list$surv[[hosp]]),
						g1      = rep(ifelse(disease.df$hosp==hosp, 
											 prob.of.exposure.to.exposed, 
											 0), times=disease.df$tte),
						g2      = unlist(censor.probs[[hosp]]),
						time    = unlist(lapply(disease.df$tte,seq.int)),
						tte     = rep(disease.df$tte, disease.df$tte),
						censor  = rep(disease.df$censor, disease.df$tte)
				))
names(hosp.dfs) <- hosps


clip <- function(x){ 
  second.lowest <- min(x[x!=min(x)])
  x[x==min(x)] <- second.lowest
  x
}

pinch.low <- function(x, low=0.025)
  ifelse(x < low, low, x)

# A very useful utility function
is.bad<-function(x) is.na(x) | is.infinite(x) | is.nan(x)

for(i in seq_along(hosp.dfs)) {
  # Some Q will equal exactly zero. This is going to be a problem when 
  # I take the inverse logit.
  hosp.dfs[[i]]$Q <- clip(hosp.dfs[[i]]$Q)
  hosp.dfs[[i]]$g1 <- ifelse(hosp.dfs[[i]]$g1==0,0,
                             1/pinch.low(hosp.dfs[[i]]$g1))
							 
  hosp.dfs[[i]]$g2 <- pinch.low(hosp.dfs[[i]]$g2)
  hosp.dfs[[i]]$clever.covariate <- with(hosp.dfs[[i]],
                                        (g1*(1/g2))*S.ratio)
  hosp.dfs[[i]]$n1 <- with(hosp.dfs[[i]], !censor & (time==tte))
}

# Seriously, I don't need the whole thing.
# Those zeroes will always just add a constant to the likelihood
# in a single variable function.
reduced.hosp.dfs <- lapply(hosp.dfs, function(x) x[x$clever.covariate!=0,])

# This is a practical consideration to prevent positivity violations.
# There will be very little data at the end of the survival curve.
# By clipping the 1% quantile, we can remove these violations, and prevent
# having enormous weights (1600) over a period of 1000 days, which overwhelms
# anything else in the model.
cutoff <- quantile(disease.df$tte,0.99)
reduced.hosp.dfs <- lapply(reduced.hosp.dfs, function(x) x[x$time<cutoff,])

# I need the epsilons for several reasons.
Q.star.by.hosp.df <- function(hosp.df) {
	model <- glm(n1~clever.covariate-1, 
              family=binomial,
              offset=qlogis(Q), 
              data=hosp.df)
    list(epsilon=coef(model), Q.star=predict(model, type='response'))
}
Q.star.and.epsilons <- lapply(reduced.hosp.dfs, Q.star.by.hosp.df)

# Setup the "initial" iteration.
Q.iteration <- lapply(Q.star.and.epsilons, function(x) matrix(x$Q.star,ncol=1))
names(Q.iteration) <- names(Q.star.and.epsilons)
epsilon.iteration <- t(sapply(Q.star.and.epsilons, "[[", "epsilon"))
colnames(epsilon.iteration) <- gsub('.clever.covariate','',colnames(epsilon.iteration),fixed=TRUE)

# Now setup for the iterations.
update.hosp.df <- function(hosp, Q.star.and.epsilons, reduced.hosp.dfs){
  Q.star <- Q.star.and.epsilons[[hosp]]$Q.star
  hosp.df <- reduced.hosp.dfs[[hosp]]
  hosp.df$Q <- Q.star
  tte.set <- disease.df$tte[disease.df$hosp==hosp]
  split.Q.star <- split(Q.star, hosp.df$id)
  new.surv <- lapply(split.Q.star, function(x) cumprod(1-x))
  hosp.df$S0 <- unlist(new.surv)
  S.ratio <- lapply(new.surv, S0.to.S.ratio)
  hosp.df$S.ratio <- unlist(S.ratio)
  hosp.df$clever.covariate <- with(hosp.df,
                                        (g1*(1/g2))*S.ratio)
  hosp.df
}

updated.hosp.dfs <- reduced.hosp.dfs
updated.Q.star.and.epsilons <- Q.star.and.epsilons
iteration <- 2
while (any((epsilon.iteration[iteration-1,]) > 1e-10)) { # Will change to delta
	updated.hosp.dfs <- lapply(hosps, update.hosp.df, 
							   updated.Q.star.and.epsilons, updated.hosp.dfs)
    names(updated.hosp.dfs) <- hosps							   
	updated.Q.star.and.epsilons <- lapply(updated.hosp.dfs, 
	                                      Q.star.by.hosp.df)
	epsilon.iteration <- rbind(epsilon.iteration,
	                           sapply(updated.Q.star.and.epsilons, "[[", "epsilon"))
    for (i in seq_along(Q.iteration))
	   Q.iteration[[i]] <- cbind(Q.iteration[[i]],
	                             updated.Q.star.and.epsilons[[i]]$Q.star)
	iteration <- iteration + 1
}
final.iteration <- iteration - 1

# Make a tte matrix.
base.tte <- do.call(cbind,outcome.S0.and.Q.list$tte)

# Make a matrix of the update order.
row.col <- as.matrix(do.call(rbind, 
                    lapply(hosps, function(hosp) 
                          data.frame(row=which(disease.df$hosp==hosp),
				                     col=match(hosp,colnames(base.tte))))))

# Extract a final Q.star
ids <- sapply(reduced.hosp.dfs, function(x) x$id)
final.Q.star <- mapply(function(x,y) split(x[,final.iteration],y), Q.iteration, ids)
# Calculate the updated tte from the final Q.star
final.updated.tte <- unlist(lapply(final.Q.star, function(x) sapply(x, 
                                                  function(y) sum(cumprod(1-y)) )))
updated.tte.matrix <- base.tte
updated.tte.matrix [row.col] <- final.updated.tte

# A couple of renames before the save
tte.mat <- updated.tte.matrix
Q.star.list <- final.Q.star
epsilon.mat <- epsilon.iteration
save(tte.mat, Q.star.list, epsilon.mat, file=output.file)
