# /usr/bin/R --args data_dump/disease_ami.object data_dump/glmnet_Q_model_ami.object 
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(survival))
registerDoMC(cores=10) # For 10 folds
options(mc.cores=12)

# arguments<-c('data_dump/disease_ami.object', 'data_dump/glmnet_Q_model_ami.object')
arguments <- commandArgs(trailingOnly=TRUE)

object.file<-arguments[1]
output.file<-arguments[2]
load(object.file)
set.seed(1)

load('data_dump/rf_G_model_ami.object')
load('survival/glmnet_Q.object')
load('survival/glmnet_g_censor.object')

# arguments <- c('data_dump/rf_G_model_ami.object','data_dump/rf_Q_model_ami.object','data_dump/glmnet_Q_model_ami.object','data_dump/rf_G_calibrated_model_ami.object','data_dump/rf_Q_calibrated_model_ami.object','data_dump/disease_ami.object', 'build_rf_Q_star_model.R','data_dump/rf_Q_star_model_ami.object','./matrix_cache')

# system.time(glmnet.predict.outcome<-cv.glmnet(x=disease.big.matrix[1:1000,], 
#                                               y=Surv(disease.df$tte, !disease.df$censor)[1:1000,],
#                                               family='cox',
#                                               parallel=TRUE))



system.time(glmnet.predict.outcome <- 
              glmnet(x=disease.big.matrix, 
                     y=Surv(disease.df$tte, 
                            !disease.df$censor),
                     family='cox' ))
save(glmnet.predict.outcome,file='survival/glmnet_Q.object')

system.time(glmnet.predict.censor <- 
              glmnet(x=disease.big.matrix, 
                     y=Surv(disease.df$tte, 
                            disease.df$censor),
                     family='cox' ))
save(glmnet.predict.censor,file='survival/glmnet_g_censor.object')



# Convergence after 37th lambda not reached. (Too small?)
# disease.big.matrix.wo.hosp <- disease.big.matrix[,!grepl('^hosp',
#                                                          colnames(disease.big.matrix))]
# 
# # Takes too long... probably just use random forest.
# system.time(glmnet.predict.exposure <- 
#               glmnet(x=disease.big.matrix.wo.hosp, 
#                      y=disease.df$hosp,
#                      family='multinomial' ))

# I hacked the survival package and scraped out the 
# survival function estimator, so I can use it with a glmnet
# found Cox model.
# Gold
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
  list(center=xcenter, coef=coef, surv=exp(-cumhaz)) # Strip the first cumhaz, it is zero.
}

get.surv.at.new.x <- function (baseline, new.x, times) {
  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
  base <- baseline$surv
  mapply(function(x,y) (base^(x))[1:y], c(scaled.new.risk), times)
}

get.hazard.at.new.x <- function (baseline, new.x, times) {
  s <- get.surv.at.new.x(baseline, new.x)
  lapply(x, function(x) diff(-log(x)))
}

S0.to.S.ratio <- function(S0)
  S0[length(S0)]/S0
 
drop.all.probs <- function(glmnet.model, y, x, s) {
  betas <- get.betas(glmnet.model, s)
  cut.x <- get.cut.x(x, betas)
  baseline <- get.kalb.prent.surv(y, x, coef=betas)

  all.zero <- cut.x
  hosps <- levels(disease.df$hosp)
  cols <- paste0('hosp', hosps)
  
  x.mat.list <- lapply(cols, function(x) {
    if (x %in% colnames(all.zero)) all.zero[,x]<-1
    all.zero
  })
  
  surv <- mclapply(x.mat.list, get.surv.at.new.x, 
                baseline=baseline, times=times)
  hazard <- mclapply(x.mat.list, get.hazard.at.new.x, 
                     baseline=baseline, times=times)
  s.ratio <- mclapply(surv, S0.to.S.ratio)
  names(s.ratio) <- names(hazard) <- names(surv) <- hosps

  # Unwrap them into data.frames.
  surv.unlist <- lapply(surv, unlist)
  hazard.unlist <- lapply(hazard, unlist)
  s.ratio.unlist <- lapply(s.ratio, unlist)
  
  lapply(hosp,
	  function(hosp) data.frame(S=surv.unlist[[hosp]], 
								Q=hazard.unlist[[hosp]], 
								S.ratio=s.ratio.unlist[[hosp]])
  )
}


# End Gold



get.haz.from.hazard.object <- function(hazard.object) {
  hazard <- hazard.object$hazard
  unique.times <- hazard.object$time
  filled.hazard <- c()
  filled.hazard[unique.times] <- hazard
  ifelse(is.na(filled.hazard),0,filled.hazard) [-1]
}

get.surv.from.hazard.object <- function(hazard.object, new.risk, times) {
  baseline <- get.cumhaz.from.hazard.object(hazard.object)
  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
  base <- exp(-baseline$cumhaz)
  mapply(function(x,y) (base^(x))[1:y], c(scaled.new.risk), times)
}


get.cumhaz.from.hazard.object <- function(hazard.object) {
  cumhaz <- hazard.object$cumhaz
  unique.times <- hazard.object$time
  cumhaz <- cumhaz[findInterval(seq.int(max(unique.times)), unique.times)]
  list(center=xcenter, coef=coef, cumhaz=cumhaz[-1]) 
}



#surv.precomputed.baseline <- function(baseline, new.x, times) {
#  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
#  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
#  base <- exp(-baseline$cumhaz)
#  mapply(function(x,y) (base^(x))[1:y], c(scaled.new.risk), times)
#}

get.hazard.object <- function(y, x, coef) {

}


survival:::agsurv

get.kalb.prent.cum.haz <- function(y, x, coef) {
  # Clean up x
  x <- as.matrix(x)
  good.indices <- complete.cases(x)
  x <- x[good.indices,,drop=FALSE]
  y <- y[good.indices,,drop=FALSE]
  xcenter <- colMeans(x)
  scaled.x <- scale(x, center = xcenter, scale = FALSE)
  scaled.risk <- c(exp(scaled.x %*% coef))
  survival:::agsurv(y=y, x=scaled.x, 
                    wt=rep(1,nrow(x)), 
                    risk=scaled.risk, 
                    survtype=3, vartype=3) # Kalbfleisch-Prentice. (Breslow)
}


  glmnet.model=glmnet.predict.outcome
  y=Surv(disease.df$tte, !disease.df$censor)
  x=disease.big.matrix
  s=72

 





Q.precomputed.baseline <- function(baseline, new.x, times) {
  scaled.new.x <- scale(new.x, center = baseline$center, scale = FALSE)
  scaled.new.risk <- exp(scaled.new.x %*% baseline$coef)
  base <- exp(-baseline$cumhaz)
  mapply(function(x,y) inv.cumprod((base^(x))[1:(y+1)]), c(scaled.new.risk), times)
}

get.betas <- function(glmnet.model, s) {
  betas <- glmnet.model$beta[,s]
  betas[betas!=0]
}

get.cut.x <- function(x, betas) {
  x[,names(betas)]
}

drop.probs <- function(glmnet.model, y, x, s, unlist=TRUE) {
  betas <- get.betas(glmnet.model, s)
  cut.x <- get.cut.x(x, betas)
  
  times <- as.matrix(y)[,1]
  baseline <- get.baseline.surv(y, cut.x, coef=betas)
  surv.all.times <- surv.precomputed.baseline(baseline, 
                                              new.x=cut.x,
                                              times=times) 
  #   m <- mapply(function(x,y) x[seq.int(y)],surv.all.times, times)
  if (unlist) return(unlist(surv.all.times))
  return(surv.all.times)
}

inv.cumprod <- function(x) x[-1] / x[-length(x)]
S0.to.Q <- function (x) 1-inv.cumprod

# Calculate g2 (probability of censoring by time)
censor.probs <-drop.probs(
  glmnet.model=glmnet.predict.censor,
  y=Surv(disease.df$tte, disease.df$censor),
  x=disease.big.matrix,
  s=12
)

# Calculate g1 (probability of exposure)
# Just going to reuse the random forest model.
# G model - calibrated RF 
g.votes <- rf.predict.exposure@oobvotes
g.by.rf.unscaled <- g.votes / rowSums(g.votes)

# Important to scale each column of g
one.v.all = sapply(colnames(g.by.rf.unscaled), function(x) disease.df$hosp == x)

g.by.rf <- mapply(function(y,x) predict(glm(y~x, family=binomial),
                                        type='response'),  
                  data.frame(one.v.all, check.names=FALSE), 
                  data.frame(g.by.rf.unscaled))

prob.of.exposure.to.exposed <- g.by.rf[cbind(
  seq.int(nrow(disease.df)),
  match(disease.df$hosp,colnames(g.by.rf)))]

# Now, get my first iteration of Q.
# Grab all the Q stuff and surv probs.
# debugonce(drop.all.probs)
outcome.S0.and.Q.list <- drop.all.probs(
  glmnet.model=glmnet.predict.outcome,
  y=Surv(disease.df$tte, !disease.df$censor),
  x=disease.big.matrix,
  s=72)

outcome.S0.list <- outcome.S0.and.Q.list$surv  
Q.init.list <- outcome.S0.and.Q.list$Q
hosps <- names(outcome.S0.list)

# outcome.S0.list[[1]][[1]])


#Q.init.list <- mclapply(outcome.S0.list,
#                        function(outcome.S0)
#                          lapply(outcome.S0, S0.to.Q))
# Why in god's name is there no USE.NAMES argument for mclapply?
#names(Q.init.list) <- hosps

S.ratio.init <- mclapply(outcome.S0.list, 
                         function(outcome.S0) 
                           lapply(outcome.S0, S0.to.S.ratio))
names(S.ratio.init) <- hosps


# sapply is super slow for reasons I can't understand.
Q.by.hosp <- do.call(cbind, lapply(Q.init.list, unlist))
S.ratio.by.hosp <- do.call(cbind, lapply(S.ratio.init, unlist))
S0.by.hosp <- do.call(cbind, lapply(outcome.S0.list, unlist))

colnames(Q.by.hosp) <- 
  colnames(S.ratio.by.hosp) <- colnames(S0.by.hosp) <- hosps

# Reorganize, so each hospital has its own data.frame
# with Q, S.ratio, S0, g1, and g2.
hosp.dfs <- lapply(hosps,
                   function(hosp) 
                     data.frame(Q       = Q.by.hosp[,hosp],
                                S.ratio = S.ratio.by.hosp[,hosp],
                                S0      = S0.by.hosp[,hosp],
                                g1      = rep(ifelse(disease.df$hosp==hosp, 
                                                     1/prob.of.exposure.to.exposed, 
                                                     0), times=disease.df$tte),
                                g2      = censor.probs,
                                time    = unlist(lapply(disease.df$tte,seq.int)),
                                tte     = rep(disease.df$tte, disease.df$tte),
                                censor  = rep(disease.df$censor, disease.df$tte)
                     ))
names(hosp.dfs) <- hosps

for(i in seq_along(hosp.dfs)) {
  hosp.dfs[[i]]$clever.covariate <- with(hosp.dfs[[i]],(g1*(1/g2))*S.ratio)
  hosp.dfs[[i]]$n1 <- with(hosp.dfs[[i]], !censor & (time==tte))
  # Some Q will equal exactly one. This is going to be a problem when 
  # I take the inverse logit.
  second.highest <- max(hosp.dfs[[i]]$Q [hosp.dfs[[i]]$Q!=1])
  hosp.dfs[[i]]$Q [hosp.dfs[[i]]$Q==1] <- second.highest
}

# I actually don't even need the epsilons, I just need the new Q*
#epsilons <- sapply(hosp.dfs,
#                   function(hosp.df) glm(n1~clever.covariate-1, offset=qlogis(Q), 
#                                         data=hosp.df)$coef['clever.covariate']
#)

# I need the epsilons for several reasons.
Q.star.by.hosp.df <- function(hosp.df) {
	model <- glm(n1~clever.covariate-1, 
              family=binomial,
              offset=qlogis(Q), 
              data=hosp.df)
    list(epsilon=coef(model), Q=predict(model, type='response'))
}

hosp.df <- hosp.dfs[[1]]
system.time(
	model <- glm(n1~clever.covariate-1, 
              family=binomial,
              offset=qlogis(Q), 
              data=hosp.df)
)

Q.star.and.epsilons <- lapply(hosp.dfs, Q.star.by.hosp.df)

Q.star <- lapply(Q.star.and.epsilons,"[[",'Q')
epsilons <- sapply(Q.star.and.epsilons,"[[",'epsilon')

match(hosps[1],disease.df$hosp)[1]
disease.df$hosp[1:2]

colnames(Q.star)

disease.df$tte[1]

hosp.df <- hosp.dfs[[1]]
hosp.df[789:(789+54-1), c('Q','S0','clever.covariate','time','censor','n1')]

predict(model, type='response') [789:(789+54-1)]

model <- glm(n1~clever.covariate-1,family=binomial, offset=qlogis(Q),data=hosp.df)

identical(min(Q.star),0)
head(which(Q.star[,6]==max(Q.star)))

Q.star[549:558,6]

head(which(Q.star[],0))

Q.star <- do.call(cbind,Q.star)
# save(Q.star, file='survival/Q.star.object')


# For some reason it is way faster to lapply and then combine.
Q.star <- do.call(cbind, Q.star)

colnames(Q.star)

Q.star.to.S.ratio <- function(Q.star) {
  starts <- cumsum(c(1,disease.df$tte[-length(disease.df$tte)]))
  ends <- c(starts[-1]-1, sum(disease.df$tte))
  
  cut.S <- mapply(function(x,y) {
                    a <- apply(Q.star[x:y,,drop=FALSE],2,cumprod)
                    dim(a) <- c(y-x+1,ncol(Q.star))
                    a
                  },
                  starts, ends, SIMPLIFY=FALSE)
  
  S.ratio <- lapply(cut.S, function(x) apply(x,2,S0.to.S.ratio))
  ret <- do.call(rbind, S.ratio)
  colnames(ret) <- colnames(Q.star)
  ret
}

hosp.df.by.S.ratio <- function(hosp.df, S.ratio, Q.star){
  hosp.df$clever.covariate <- with(hosp.df, (g1*(1/g2))*S.ratio)
  hosp.df$Q <- Q.star
  hosp.df
}

new.S.ratio <- Q.star.to.S.ratio(Q.star)

head(new.S.ratio[,'Hôpital Lanaudière'],100)

head(hosp.dfs[['Hôpital Lanaudière']])
head(new.hosp.dfs[['Hôpital Lanaudière']])


new.hosp.dfs <- mapply(hosp.df.by.S.ratio, 
                       hosp.dfs, 
                       data.frame(new.S.ratio),
                       data.frame(Q.star), SIMPLIFY=FALSE)



disease.df$hosp[1]


new.Q.star <- lapply(new.hosp.dfs, Q.star.by.hosp.df)





disease.df[145,]

cut.S[[145]]

which(sapply(sapply(cut.S,dim),is.null))

dim(cut.S[[1]])

S0.to.S.ratio <- function(S0)
  S0[length(S0)]/S0


S <- do.call(rbind,cut.S)



cut.S[[1]]

which(starts==76550)

which(disease.df$tte==0)

starts[145]
ends[145]

disease.df[143:147,]

?mapply

head(starts)
head(ends)

ends <- c(0,disease.df$tte[-nrow(nrow(disease.df))])

coef.val <- glm(n1~clever.covariate-1, offset=qlogis(Q), data=hosp.dfs[[1]])$coef['clever.covariate']



# Q.init <- mclapply(outcome.S0, S0.to.Q)
# S.ratio.init <- mclapply(outcome.S0, S0.to.S.ratio)

init.mat <- data.frame(Q=unlist(Q.init), 
                       S.ratio=unlist(S.ratio.init), S0=unlist(outcome.S0), 
                       g1=rep(prob.of.exposure.to.exposed, times=disease.df$tte), 
                       g2=censor.probs)

init.mat$clever.covariate <- with(init.mat, S.ratio / (g1 * g2))
init.mat$time <- unlist(lapply(disease.df$tte,seq.int))
init.mat$tte <- rep(disease.df$tte, disease.df$tte)
init.mat$censor <- rep(disease.df$censor, disease.df$tte)
init.mat$n1 <- with(init.mat,!censor & (time==tte))

debugonce(glm.fit)

efficient.influence.func <- glm(n1~Q+clever.covariate-1,
                                data=init.mat,
                                family=binomial)



p=(predict(glmnet.predict.exposure, newx=disease.big.matrix, type='response'))
s=9
exposure.mat <- p[,,s]

extract.betas <- function(model,outcome,s)
  which(model$beta[[outcome]][,s]!=0)


which(glmnet.predict.exposure$beta[['Hôpital Pierre-Le Gardeur']][,s]!=0)

sapply(unique(disease.df$hosp),
       extract.betas, model=glmnet.predict.exposure,s=9)

table(exposure.mat[,1], big.hosp.mat[,1])

platt.scale <- function(y,x) predict(glm(y~x,family=binomial), type='response')

# Fix up an dummy variable matrix of the exposures.
grep('^hosp',colnames(disease.big.matrix), value=TRUE)
big.hosp.mat <- model.matrix(~hosp-1,disease.df['hosp'])
colnames(big.hosp.mat) <- gsub('^hosp','',colnames(big.hosp.mat))
big.hosp.mat <- big.hosp.mat[,colnames(exposure.mat)]

# Platt scale all the predictions against the dummy variables.
scaled.exposure.prob <- mapply(platt.scale, data.frame(big.hosp.mat), data.frame(exposure.mat))

colnames(big.hosp.mat)
testo=as.numeric(disease.df$hosp=='Hôpital Pierre-Le Gardeur')


which(glmnet.predict.exposure$beta$'Hôpital Pierre-Le Gardeur'[,23]!=0)

(which(glmnet.predict.exposure$beta[,s]!=0))

glm(big.hosp.mat[,1] ~ exposure.mat[,1], family=binomial)



dim(exposure.mat)

str(glmnet.predict.exposure)

censor.probs <-drop.probs(
  y=Surv(disease.df$tte, disease.df$censor),
  x=disease.big.matrix,
  s=12,
  times=disease.df$tte)

length(censor.probs)

disease.df$tte

glmnet.predict.censor

x=disease.big.matrix 
y=Surv(disease.df$tte, disease.df$censor)
s=12



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