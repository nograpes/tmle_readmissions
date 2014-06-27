suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.
options(mc.cores=12)

# /usr/bin/R ./build_rf_Q_star_model.R --args data_dump/rf_G_model_heart_failure.object data_dump/rf_Q_model_heart_failure.object data_dump/glmnet_Q_model_heart_failure.object data_dump/rf_G_calibrated_model_heart_failure.object data_dump/rf_Q_calibrated_model_heart_failure.object data_dump/disease_heart_failure.object  data_dump/build_rf_Q_star_model.R data_dump/rf_Q_star_model_heart_failure.object  ./matrix_cache

# arguments <- c('data_dump/rf_G_model_pneumonia.object','data_dump/rf_Q_model_pneumonia.object','data_dump/glmnet_Q_model_pneumonia.object','data_dump/rf_G_calibrated_model_pneumonia.object','data_dump/rf_Q_calibrated_model_pneumonia.object','data_dump/disease_pneumonia.object', 'build_rf_Q_star_model.R','data_dump/rf_Q_star_model_pneumonia.object','./matrix_cache')

# /usr/bin/R ./build_rf_Q_star_model.R  --args data_dump/rf_G_model_ami.object data_dump/rf_Q_model_ami.object data_dump/glmnet_Q_model_ami.object data_dump/rf_G_calibrated_model_ami.object data_dump/rf_Q_calibrated_model_ami.object data_dump/disease_ami.object data_dump/rf_Q_star_model_ami.object ./matrix_cache

arguments<-commandArgs(trailingOnly=TRUE)

# G.model.file <- arguments[1]
# rf.Q.model.file <- arguments[2]
glmnet.Q.model.file <- arguments[3]
calibrated.G.model.file <- arguments[4]
calibrated.rf.Q.model.file <- arguments[5]
object.file  <- arguments[6]
R.file  <- arguments[7] # Ugly hack -- I can't get the Makefile to work without passing the R file as a param.
output.file  <- arguments[8]
matrix.cache <- arguments[9]

load(object.file)
load(glmnet.Q.model.file) # glmnet.predict.outcome
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
fixed.hosp.data <- function(hosp){
  set.to.hosp <- disease.big.matrix # Fix A to the hospital.
  set.to.hosp[,grep('^hosp',colnames(disease.big.matrix))] <- 0
  var.name <- paste0('hosp',hosp)
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name] <- 1
  set.to.hosp
}

rf.prob.by.hosp <- function(rf.model, hosp, oob.only=TRUE)
  rf.probs(rf.model=rf.model, data=fixed.hosp.data(hosp), oob.only=oob.only)

glmnet.prob.by.hosp<-function(hosp){
  set.to.hosp<-disease.big.matrix
  set.to.hosp[,grep('^hosp',colnames(disease.big.matrix))]<-0
  var.name<-paste0('hosp',hosp)
  # Because one hospital was left one -- just by tradition for interpretation.
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name]<-1
  c(plogis(predict(glmnet.predict.outcome, 
                   newx=set.to.hosp, 
              		    s=glmnet.predict.outcome$lambda.1se)))  
}

# G model - calibrated RF 
votes <- rf.predict.exposure@oobvotes
rf.prob.of.exposure <- votes/rowSums(votes)
# rf.prob.of.exposure <- prob[cbind(1:length(disease.df$hosp),
#                             match(disease.df$hosp,colnames(prob)))]

# Q model - calibrated RF - as observed (A=a)
votes <- rf.predict.outcome@oobvotes
prob <- votes/rowSums(votes)
vote.prop <- prob[,"TRUE"]
# Very important to scale the vote proportions.
platt.scaler <- 
  glm(disease.df$day_30_readmit ~ vote.prop, family=binomial(link='logit'))
rf.prob.of.outcome <-  predict(platt.scaler, type='response')

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

# Q model - glmnet - as observed (A=a)
glmnet.prob.of.outcome <-c(plogis(predict(glmnet.predict.outcome, 
                                          newx=disease.big.matrix, 
                                          s=glmnet.predict.outcome$lambda.1se)))  

# Q model - glmnet - manipulating exposure to each of twenty levels
# There is almost certainly a more clever way to do this:
# Calculate for baseline, and then just add in.
all.glmnet.Q.by.hosp <- sapply(levels(disease.df$hosp),glmnet.prob.by.hosp)

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
modified.rf.prob.of.exposure <- bump.zeroes(rf.prob.of.exposure)

# Outcome
modified.rf.prob.of.outcome <- bump.zeroes(rf.prob.of.outcome)

# I want to find one epsilon for each hospital.
# Each column becomes (A==a)/g
# iptw <- exposure.mat  / modified.rf.prob.of.exposure 
iptw <- exposure.mat  / modified.rf.prob.of.exposure 

# You don't have to run a model for each level, since all of the iptw vars are mutually exclusive (when one is nonzero, all rest are zero) putting them all in the same model is essentially the same thing.

# Standard GLM fitter.
# rf.epsilon.model<-glm(Y ~ .-1, 
#                       offset=qlogis(modified.rf.prob.of.outcome), 
#                       data=data.frame(Y=as.factor(disease.df$day_30_readmit),iptw),
#                       family=binomial(link=logit))
# 
# 
# glmnet.epsilon.model<-glm(Y ~ .-1, 
#                       offset=qlogis(glmnet.prob.of.outcome), 
#                       data=data.frame(Y=as.factor(disease.df$day_30_readmit),iptw),
#                       family=binomial(link=logit))
# rf.epsilons<-coef(rf.epsilon.model)
# glmnet.epsilons<-coef(glmnet.epsilon.model)

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
# rf.epsilons <- epsilons(offset=modified.rf.prob.of.outcome, iptw=iptw)
# glmnet.epsilons <- epsilons(offset=glmnet.prob.of.outcome, iptw=iptw)
rf.epsilons <- mapply(epsilons, offset=data.frame(all.rf.Q.by.hosp) ,iptw=data.frame(iptw))
glmnet.epsilons <- mapply(epsilons, offset=data.frame(all.glmnet.Q.by.hosp) ,iptw=data.frame(iptw))

Q.star<-function(Q, iptw, epsilons)
				plogis(qlogis(Q) + ((1/iptw) %*% epsilons))

mapply(Q, iptw, )	
	
# rf.Q.star <- Q.star(all.rf.Q.by.hosp, 
#                     iptw=iptw, 
# 	 				  rf.epsilons)

# glmnet.Q.star <- Q.star(all.glmnet.Q.by.hosp, 
#                         iptw=iptw, 
#                         glmnet.epsilons)

save(rf.Q.star, glmnet.Q.star,
     rf.epsilons, glmnet.epsilons,
     all.rf.Q.by.hosp, all.glmnet.Q.by.hosp,
  	 file=output.file)

