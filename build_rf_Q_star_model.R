suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.
options(mc.cores=12)

arguments<-commandArgs(trailingOnly=TRUE)
if(interactive()) arguments<-c('data_dump/rf_G_calibrated_model_ami.object','data_dump/rf_Q_calibrated_model_ami.object','data_dump/disease_ami.object','build_rf_Q_star_model.R','data_dump/rf_Q_star_model_ami.object','./matrix_cache')

calibrated.G.model.file <- arguments[1]
calibrated.rf.Q.model.file <- arguments[2]
object.file  <- arguments[3]
R.file  <- arguments[4] # Ugly hack -- I can't get the Makefile to work without passing the R file as a param.
output.file  <- arguments[5]
matrix.cache <- arguments[6]

load(object.file)
load(calibrated.G.model.file)
load(calibrated.rf.Q.model.file)

# For a given tree, get predicted class for data
predict.by.tree<-function(tree.num, forest, data, cachepath=matrix.cache, n){
  xtype <- as.integer(.Call("CGetType", data@address, PACKAGE = "bigmemory"))
  tree<-forest[[tree.num]]
  treepredict.result <- .Call('treepredictC',
                              data@address, xtype, n, forest, tree, 
                              PACKAGE = "bigrf")
  sapply(seq_along(levels(forest@y)),function(c){
    w <- treepredict.result$testpredclass == c
    ifelse(w,tree@nodewt[treepredict.result$testprednode],0)
  })  
}

# This will calculate vote proportions for any rf model
rf.probs <- function(rf.model, data, oob.only=TRUE){
  data.pointer <- bigrf:::makex(data, "xtest", cachepath=matrix.cache)
  # Use foreach instead of 
  prediction.by.tree<-foreach(tree=seq_along(rf.model)) %dopar% 
    predict.by.tree(tree.num=tree, forest=rf.model, data=data.pointer, n=nrow(data))

  # Now set all in-bag to NA
  if(oob.only){
    # Assume provided data is the first n of training sample.
    in.bag.mat<-sapply(rf.model, function(tree) tree@insamp!=0)[1:nrow(data),]
    for(col in 1:ncol(in.bag.mat)){
      in.bag <- in.bag.mat[,col]
      prediction.by.tree[[col]][in.bag,] <- 0
    }
    # prediction.by.tree[in.bag.mat]<-NA
  }
  votes <- Reduce("+", prediction.by.tree)
  classes <- names(rf.model@ytable)
  colnames(votes) <- classes
  votes / rowSums(votes)
}

# This returns a disease matrix with the exposure fixed to a certain hospital.
fixed.hosp.data <- function(hosp) {
  set.to.hosp <- disease.big.matrix # Fix A to the hospital.
  set.to.hosp[,grep('^hosp',colnames(disease.big.matrix))] <- 0
  var.name <- paste0('hosp',hosp)
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name] <- 1
  set.to.hosp
}

rf.prob.by.hosp <- function(rf.model, hosp, oob.only=TRUE)
  rf.probs(rf.model=rf.model, data=fixed.hosp.data(hosp), oob.only=oob.only)

# G model - calibrated RF 
g.votes <- rf.predict.exposure@oobvotes
g.by.rf.unscaled <- g.votes / rowSums(g.votes)

# Important to scale each column of g
one.v.all = sapply(colnames(g.by.rf.unscaled), function(x) disease.df$hosp == x)

g.by.rf <- mapply(function(y,x) predict(glm(y~x, family=binomial),
                                type='response'),  
                  data.frame(one.v.all, check.names=FALSE), 
                  data.frame(g.by.rf.unscaled))

# Q model - calibrated RF - as observed (A=a)
votes <- rf.predict.outcome@oobvotes
prob <- votes/rowSums(votes)
vote.prop <- prob[,"TRUE"]

# Very important to scale the vote proportions.
platt.scaler <- 
  glm(disease.df$day_30_readmit ~ vote.prop, family=binomial(link='logit'))

Q.as.observed.by.rf <-  predict(platt.scaler, type='response')

# Q model - calibrated RF - manipulating exposure to each of twenty levels.
# (A=1, A=2..,A=20)
unscaled.all.rf.Q.by.hosp<-sapply(levels(disease.df$hosp),
                                  function(x,rf.model) 
                                    rf.prob.by.hosp(rf.model,x)[,"TRUE"],
                                  rf.model=rf.predict.outcome)

all.rf.Q.by.hosp <- apply(unscaled.all.rf.Q.by.hosp,2, 
                           function(x) predict(platt.scaler,
                                               newdata=data.frame(prob=x),
                                               type='response'))

# Create a matrix of indicator variables.
ff<-~hosp-1
exposure.mat<-(model.matrix(ff,model.frame(ff, disease.df)))
colnames(exposure.mat)<-gsub('^hosp','',colnames(exposure.mat))

# Temporary fix:
# For seven people, the probability of exposure (going to the hospital the went to) is exactly zero.
# I'm going to bump it up slightly as a quick fix.
# But this should be subject to a lot of analysis.
bump.zeroes <- function(x){
  min.baseline<-min(x[x!=0])
  ifelse(x==0, min.baseline, x)
}

# Exposure
modified.rf.prob.of.exposure <- bump.zeroes(g.by.rf)

# Outcome
# modified.rf.prob.of.outcome <- bump.zeroes(rf.prob.of.outcome)

# I want to find one epsilon for each hospital.
# Each column becomes (A==a)/g
iptw <- exposure.mat  / modified.rf.prob.of.exposure 

# You don't have to run a model for each level, since all of the iptw vars are mutually exclusive (when one is nonzero, all rest are zero) putting them all in the same model is essentially the same thing.

# Standard GLM fitter.
# rf.epsilon.model<-glm(Y ~ .-1, 
#                       offset=qlogis(modified.rf.prob.of.outcome), 
#                       data=data.frame(Y=as.factor(disease.df$day_30_readmit),iptw),
#                       family=binomial(link=logit))
# 
# 
# rf.epsilons<-coef(rf.epsilon.model)

# Normally, I would use the standard GLM fitter for these data.
# It appears that sometimes the fitter, which uses Iterative Reweighted Least Squares (IRLS)
# to fit the betas doesn't always fit well.
# Indeed, for pneumonia, the one variable kept flipping back and forth, and would
# never converge, despite having an extremely smooth, convex likelihood function.
# Using BFGS, I have much better chance of convergence, and slightly better fits.
# It is a little slow but who cares?
# Plus the names don't get mucked up because of rowname character restrictions.
glm.BFGS<-function(x,y,offset=rep(plogis(0),length(Y))) {
  likelihood<-function(betas) {
    preds<-plogis(((x %*% t(betas)) + qlogis(offset)))
    sum((y * log(preds)) + ((1-y)*log(1-preds)))
  }
  setNames(optim(rep(0,ncol(as.matrix(x))), 
                 function(betas) abs(likelihood(betas)), method='BFGS')$par,
           colnames(x))
}

epsilons<-function(offset, iptw) {
  glm.BFGS(x=iptw,
           y=as.numeric(disease.df$day_30_readmit),
           offset=offset)
}

# Offset needs to change every time.
# IPTW needs to change every time.

# I have chosen to use only the calibrated IPTW, because it is more theoretically sound. (I really don't care about accuracy here.)
rf.epsilons <- mapply(epsilons, offset=data.frame(all.rf.Q.by.hosp) ,iptw=data.frame(iptw))

Q.star<-function(Q, ptw, epsilons)
				plogis(qlogis(Q) + ((ptw) %*% t(epsilons)))

rf.Q.star <- mapply(Q.star, Q=data.frame(all.rf.Q.by.hosp), ptw=data.frame(modified.rf.prob.of.exposure), epsilons=rf.epsilons)	

save(rf.Q.star, rf.epsilons, all.rf.Q.by.hosp, file=output.file)

