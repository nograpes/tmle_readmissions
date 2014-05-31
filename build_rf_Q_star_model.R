suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.

# /usr/bin/R ./build_rf_Q_star_model.R --args data_dump/rf_G_model_pneumonia.object data_dump/rf_Q_model_pneumonia.object data_dump/glmnet_Q_model_pneumonia.object data_dump/rf_G_calibrated_model_pneumonia.object data_dump/rf_Q_calibrated_model_pneumonia.object data_dump/disease_pneumonia.object data_dump/rf_Q_star_model_pneumonia.object ./matrix_cache

# arguments <- c('data_dump/rf_G_model_pneumonia.object','data_dump/rf_Q_model_pneumonia.object','data_dump/glmnet_Q_model_pneumonia.object','data_dump/rf_G_calibrated_model_pneumonia.object','data_dump/rf_Q_calibrated_model_pneumonia.object','data_dump/disease_pneumonia.object','data_dump/rf_Q_star_model_pneumonia.object','./matrix_cache')

# /usr/bin/R ./build_rf_Q_star_model.R  --args data_dump/rf_G_model_ami.object data_dump/rf_Q_model_ami.object data_dump/glmnet_Q_model_ami.object data_dump/rf_G_calibrated_model_ami.object data_dump/rf_Q_calibrated_model_ami.object data_dump/disease_ami.object data_dump/rf_Q_star_model_ami.object ./matrix_cache

arguments<-commandArgs(trailingOnly=TRUE)

G.model.file <- arguments[1]
rf.Q.model.file <- arguments[2]
glmnet.Q.model.file <- arguments[3]
calibrated.G.model.file <- arguments[4]
calibrated.rf.Q.model.file <- arguments[5]
object.file  <- arguments[6]
output.file  <- arguments[7]
R.file  <- arguments[8] # Ugly hack -- I can't get the Makefile to work without passing the R file as a param.
matrix.cache <- arguments[9]

load(object.file)
load(G.model.file) # rf.predict.exposure
load(rf.Q.model.file) # rf.predict.outcome
load(glmnet.Q.model.file) # glmnet.predict.outcome

# Since I stupidly named the variable the same in the calibrated versions, I'll have to load into an environment to avoid overwriting (masking).
e=new.env()
load(calibrated.G.model.file,e)
load(calibrated.rf.Q.model.file,e)
calibrated.rf.predict.exposure<-e[["rf.predict.exposure"]]
calibrated.rf.predict.outcome<-e[["rf.predict.outcome"]]
rm(e)


rf.prob.by.hosp<-function(hosp, rf.predict.outcome){
  set.to.hosp<-disease.big.matrix
  # Fix A to the hospital.
  set.to.hosp[,grep('^hosp',colnames(disease.big.matrix))]<-0
  var.name<-paste0('hosp',hosp)
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name]<-1
  
  predict.by.tree<-function(tree.num,forest,data,cachepath=matrix.cache,n){
    xtype <- as.integer(.Call("CGetType", data@address, PACKAGE = "bigmemory"))
    tree<-forest[[tree.num]]
    .Call('treepredictC',data@address,xtype,n,forest,tree,PACKAGE = "bigrf")$testpredclass 
  }
  data <- bigrf:::makex(set.to.hosp, "xtest", cachepath=matrix.cache)
  # Use foreach instead of 
  prediction.by.tree<-foreach(i=seq_along(rf.predict.outcome),.combine=cbind) %dopar% predict.by.tree(tree=i,forest=rf.predict.outcome,data=data,n=nrow(set.to.hosp))
  prediction.by.tree<-prediction.by.tree == 2 # Because class 2 is true
  # Now set all in-bag to NA
  in.bag.mat<-sapply(rf.predict.outcome,function(tree)tree@insamp!=0)
  prediction.by.tree[in.bag.mat]<-NA
  (rowSums(prediction.by.tree,na.rm=TRUE) / rowSums(!in.bag.mat)[1])
}

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

# G model - RF 
votes<-rf.predict.exposure@oobvotes
prob<-votes/rowSums(votes)
rf.prob.of.exposure <- prob[cbind(1:length(disease.df$hosp),
                                  match(disease.df$hosp,colnames(prob)))]

# G model - Calibrated RF								  
votes<-calibrated.rf.predict.exposure@oobvotes
prob<-votes/rowSums(votes)
calibrated.rf.prob.of.exposure <- prob[cbind(1:length(disease.df$hosp),
                                       match(disease.df$hosp,colnames(prob)))]
								  
								  
# Q model - RF - as observed (A=a)
votes<-rf.predict.outcome@oobvotes
prob<-votes/rowSums(votes)
rf.prob.of.outcome <- prob[,"TRUE"]

# Q model - calibrated RF - as observed (A=a)
votes<-calibrated.rf.predict.outcome@oobvotes
prob<-votes/rowSums(votes)
calibrated.rf.prob.of.outcome <- prob[,"TRUE"]


# Q model - glmnet - as observed (A=a)
glmnet.prob.of.outcome <-c(plogis(predict(glmnet.predict.outcome, 
                                          newx=disease.big.matrix, 
										  s=glmnet.predict.outcome$lambda.1se)))  

# Q model - RF - manipulating exposure to each of twenty levels. (A=1, A=2..,A=20)
all.rf.Q.by.hosp<-sapply(levels(disease.df$hosp),rf.prob.by.hosp,rf.predict.outcome=rf.predict.outcome)

# Q model - calibrated RF - manipulating exposure to each of twenty levels. (A=1, A=2..,A=20)
all.calibrated.rf.Q.by.hosp<-sapply(levels(disease.df$hosp),rf.prob.by.hosp,rf.predict.outcome=calibrated.rf.predict.outcome)


# Q model - glmnet - manipulating exposure to each of twenty levels
# There is almost certainly a more clever way to do this:
# Calculate for baseline, and then just add in.
all.glmnet.Q.by.hosp<-sapply(levels(disease.df$hosp),glmnet.prob.by.hosp)

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
modified.calibrated.rf.prob.of.exposure <- bump.zeroes(calibrated.rf.prob.of.exposure)

# Outcome
modified.rf.prob.of.outcome <- bump.zeroes(rf.prob.of.outcome)
modified.calibrated.rf.prob.of.outcome <- bump.zeroes(calibrated.rf.prob.of.outcome)

# I want to find one epsilon for each hospital.
# Each column becomes (A==a)/g
iptw <- exposure.mat  / modified.rf.prob.of.exposure 
calibrated.iptw <- exposure.mat  / modified.calibrated.rf.prob.of.exposure 

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
# Indeed, for pneumonia, the Charles-Lemoyne variable kept flipping back and forth, and would
# never converge, despite having an extremely smooth, convex likelihood function.
# Using BFGS, I have much better chance of convergence, and slightly better fits.
# It is a little slow but who cares?
# Plus the names don't get mucked up because of rowname character restrictions.
glm.BFGS<-function(x,y,offset=rep(plogis(0),length(Y))) {
  likelihood<-function(betas) {
    preds<-plogis(((x %*% betas) + qlogis(offset)))
    sum((y * log(preds)) + ((1-y)*log(1-preds)))
  }
  setNames(optim(rep(0,ncol(x)), function(betas) abs(likelihood(betas)), method='BFGS')$par,
           colnames(x))
}

epsilons<-function(offset, iptw) {
  glm.BFGS(x=iptw,
           y=as.numeric(disease.df$day_30_readmit),
           offset=offset)
}

# I have chosen to use only the calibrated IPTW, because it is more theoretically sound. (I really don't care about accuracy here.)
rf.epsilons <- epsilons(modified.rf.prob.of.outcome, calibrated.iptw)
calibrated.rf.epsilons <- epsilons(modified.calibrated.rf.prob.of.outcome, calibrated.iptw)
glmnet.epsilons <- epsilons(glmnet.prob.of.outcome, calibrated.iptw)

Q.star<-function(Q, iptw, epsilons)
				plogis(qlogis(Q) + ((1/iptw) %*% t(epsilons)))

rf.Q.star <- Q.star(all.rf.Q.by.hosp, 
                    modified.rf.prob.of.outcome, 
					rf.epsilons)

calibrated.rf.Q.star <- Q.star(all.calibrated.rf.Q.by.hosp, 
                               modified.calibrated.rf.prob.of.outcome, 
					           calibrated.rf.epsilons)

glmnet.Q.star <- Q.star(all.glmnet.Q.by.hosp, 
                        glmnet.prob.of.outcome, 
                        glmnet.epsilons)

						
# sapply(list(rf.Q.star, 
#             calibrated.rf.Q.star, 
# 			 calibrated.glmnet.Q.star), colMeans)
					
						
# rf.Q.star<-plogis(qlogis(all.rf.Q.by.hosp) + 
#                    ((1/modified.rf.prob.of.outcome) %*% t(rf.epsilons)))

# calibrated.rf.Q.star<-plogis(qlogis(all.rf.Q.by.hosp) + 
#                    ((1/modified.calibrated.rf.prob.of.outcome) %*%  t(rf.epsilons)))

# glmnet.Q.star<-plogis(qlogis(all.glmnet.Q.by.hosp) + 
#                    ((1/glmnet.prob.of.outcome) %*% t(glmnet.epsilons)))
				   
# save(rf.epsilons, glmnet.epsilons, rf.Q.star, glmnet.Q.star, file=output.file)

save(rf.Q.star, calibrated.rf.Q.star, glmnet.Q.star,
     rf.epsilons, calibrated.rf.epsilons, glmnet.epsilons,
     all.rf.Q.by.hosp, all.calibrated.rf.Q.by.hosp, all.glmnet.Q.by.hosp,
	 file=output.file)

