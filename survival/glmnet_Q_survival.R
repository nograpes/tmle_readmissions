suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(survival))
registerDoMC(cores=10) # For 10 folds
options(mc.cores=12)

arguments<-c('data_dump/disease_ami.object','data_dump/rf_G_model_ami.object', 'survival/glmnet_Q.object', 'survival/glmnet_g_censor.object')
# arguments <- commandArgs(trailingOnly=TRUE)

object.file <- arguments[1]
rf.G.object <- arguments[2]
survival.Q.object <- arguments[3]
survival.censor.object <- arguments[4]

load(object.file)
load(rf.G.object)
load(survival.Q.object)
load(survival.censor.object)

# I hacked the survival package and scraped out the 
# survival function estimator, so I can use it with a glmnet
# found Cox model.
get.betas <- function(glmnet.model, s) {
  betas <- glmnet.model$beta[,s]
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
  ret <- mapply(function(x,y) (base^(x))[1:y], c(scaled.new.risk), times)
  # By definition, you can't have an estimate at the very end of the 
  # longest time. So, you'll get NAs for those. Reset those to the last
  # estimated survival time. 
  lapply(ret, function(x) {
    if (length(x)>=baseline$max) {
      last.not.NA <- x[max(which(!is.na(x)))]
	  x[baseline$max:length(x)] <- last.not.NA
    }
	x
  })
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
 
drop.all.probs <- function(glmnet.model, y, x, s) {
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
  surv <- mclapply(x.mat.list, get.surv.at.new.x, 
                baseline=baseline, times=times+1) # Because first is always 1
  hazard <- mclapply(surv, get.hazard.at.new.x)
  cond.fail <- mclapply(surv, get.conditional.failure.at.new.x)
  surv <- lapply(surv, function(x) lapply(x, function(y) y[-1])) # Now clip the first surv prob.
  s.ratio <- mclapply(surv, function(x) lapply(x,S0.to.S.ratio))
  names(cond.fail) <- names(s.ratio) <- names(hazard) <- names(surv) <- hosps
  list(surv=surv, cond.fail=cond.fail, hazard=hazard, s.ratio=s.ratio)
}

# S0 and Q.
# debugonce(drop.all.probs)
outcome.S0.and.Q.list <- drop.all.probs(
  glmnet.model=glmnet.predict.outcome,
  y=Surv(disease.df$tte, !disease.df$censor),
  x=disease.big.matrix,
  s=72)
hosps=names(outcome.S0.and.Q.list$surv)

# Calculate g2 (probability of censoring by time)
censor.probs <-drop.all.probs(
  glmnet.model=glmnet.predict.censor,
  y=Surv(disease.df$tte, disease.df$censor),
  x=disease.big.matrix,
  s=12
)$surv
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
				data.frame(Q    = unlist(outcome.S0.and.Q.list$cond.fail[[hosp]]),
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

# I need the epsilons for several reasons.
Q.star.by.hosp.df <- function(hosp.df) {
	model <- glm(n1~clever.covariate-1, 
              family=binomial,
              offset=qlogis(Q), 
              data=hosp.df)
    list(epsilon=coef(model), Q.star=predict(model, type='response'))
}
Q.star.and.epsilons <- lapply(reduced.hosp.dfs, Q.star.by.hosp.df)


mean(Q.star.and.epsilons[[1]]$Q / Q.star.and.epsilons[[2]]$Q)

length(Q.star.and.epsilons[[1]]$Q)
length(Q.star.and.epsilons[[2]]$Q)

sapply(Q.star.and.epsilons, "[[", "epsilon")

hosp.df = hosp.dfs[["Hôpital Saint-Luc du CHUM"]]
hosp.df[which(hosp.df$clever.covariate==max(hosp.df$clever.covariate)),]


m = Q.star.by.hosp.df(hosp.df)
m$epsilon

hosp.df = hosp.df[hosp.df$clever.covariate!=0,]
m2 = Q.star.by.hosp.df.test.offset(hosp.df)

summary(n1 - hosp.df$Q)

guy.54 <- hosp.dfs[[1]][789:(789+54-1),]

cumprod(1 - (guy.54$Q + Q.star.and.epsilons[[1]]$epsilon))
cumprod(1 - guy.54$Q)


head(hosp.df)
summary(hosp.df$clever.covariate)

sapply(hosp.dfs, function(x) max(x$clever.covariate))

1/50079
1/0.025

system.time(model <- glm(n1~clever.covariate-1, 
              family=binomial,
              offset=qlogis(Q), 
              data=hosp.dfs[[1]]))

# End Gold





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