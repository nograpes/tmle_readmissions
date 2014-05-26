suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.

# /usr/bin/R --args data_dump/rf_G_model_heart_failure.object data_dump/rf_Q_model_heart_failure.object data_dump/glmnet_Q_model_heart_failure.object data_dump/disease_heart_failure.object data_dump/rf_Q_star_model_heart_failure.object matrix_cache
# /usr/bin/R --args data_dump/rf_G_model_ami.object data_dump/rf_Q_model_ami.object data_dump/disease_ami.object data_dump/rf_Q_star_model_ami.object
# setwd('~/repo/thesis/code/tmle')
# arguments<-paste('data_dump',c('rf_G_model_pneumonia.object','rf_Q_model_pneumonia.object','disease_pneumonia.object','rf_Q_star_model_pneumonia.object'),sep='/')
# arguments<-c(paste('data_dump',c('rf_G_model_ami.object','rf_Q_model_ami.object','disease_ami.object','rf_Q_star_model_ami.object'),sep='/'),'matrix_cache')
arguments<-commandArgs(trailingOnly=TRUE)

G.model.file <- arguments[1]
rf.Q.model.file <- arguments[2]
glmnet.Q.model.file <- arguments[3]
object.file  <- arguments[4]
output.file  <- arguments[5]
matrix.cache <- arguments[6]
load(object.file)
load(G.model.file)
load(rf.Q.model.file) # rf.predict.outcome
load(glmnet.Q.model.file) # glmnet.predict.outcome

rf.prob.by.hosp<-function(hosp){
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

# Q model - RF - as observed (A=a)
votes<-rf.predict.outcome@oobvotes
prob<-votes/rowSums(votes)
rf.prob.of.outcome <- prob[,"TRUE"]
# Q model - glmnet - as observed (A=a)
glmnet.prob.of.outcome <-c(plogis(predict(glmnet.predict.outcome, 
                                          newx=disease.big.matrix, 
										  s=glmnet.predict.outcome$lambda.1se)))  

# Q model - RF - manipulating exposure to each of twenty levels. (A=1, A=2..,A=20)
all.rf.Q.by.hosp<-sapply(levels(disease.df$hosp),rf.prob.by.hosp)
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
min.baseline<-min(rf.prob.of.exposure[rf.prob.of.exposure!=0])
modified.rf.prob.of.exposure<-ifelse(rf.prob.of.exposure==0,min.baseline,rf.prob.of.exposure)

# w=which(glmnet.predict.outcome$lambda.1se == glmnet.predict.outcome$lambda)
# non.zero=glmnet.predict.outcome$glmnet.fit$beta[(which(glmnet.predict.outcome$glmnet.fit$beta[ ,w]!=0)),w]
# t(t(exp(sort(non.zero))))
# Totally different, but not the same because it includes size of change.
# tail(sort(fastimp(rf.predict.outcome)),100)

# The same for the outcome
min.baseline<-min(rf.prob.of.outcome[rf.prob.of.outcome!=0])
modified.rf.prob.of.outcome<-
  ifelse(rf.prob.of.outcome==0,min.baseline,rf.prob.of.outcome)

# I want to find one epsilon for each hospital.
iptw <- exposure.mat  / modified.rf.prob.of.exposure # Each column is now the (A==a)/g

# You don't have to run a model for each level, since all of the iptw vars are mutually exclusive (when one is nonzero, all rest are zero) putting them all in the same model is essentially the same thing.

data.frame(Y=as.factor(disease.df$day_30_readmit),iptw)

rf.epsilon.model<-glm(Y ~ .-1, 
                      offset=qlogis(modified.rf.prob.of.outcome), 
                      data=data.frame(Y=as.factor(disease.df$day_30_readmit),iptw),
                      family=binomial(link=logit))

glmnet.epsilon.model<-glm(Y ~ .-1, 
                      offset=qlogis(glmnet.prob.of.outcome), 
                      data=data.frame(Y=as.factor(disease.df$day_30_readmit),iptw),
                      family=binomial(link=logit))

rf.epsilons<-coef(rf.epsilon.model)
glmnet.epsilons<-coef(glmnet.epsilon.model)
					  
# epsilon<-function(hospital.covariate) 
#   glm(as.factor(disease.df$day_30_readmit) ~ 
#         hospital.covariate-1, 
#       offset=qlogis(modified.prob.of.outcome), 
#       family=binomial(link=logit))$coef

	  
# epsilons<-apply(iptw,2,epsilon) # The fitter complains that fitted probs of 0 or 1 occurred.

Q.star=plogis(qlogis(all.rf.Q.by.hosp[,hosp]) + (rf.epsilons[hosp] * (1/modified.rf.prob.of.outcome)))
Q=all.rf.Q.by.hosp[,hosp] 


rf.Q.star<-plogis(qlogis(all.rf.Q.by.hosp) + 
                   ((1/modified.rf.prob.of.outcome) %*% t(rf.epsilons)))

glmnet.Q.star<-plogis(qlogis(all.glmnet.Q.by.hosp) + 
                     ((1/glmnet.prob.of.outcome) %*% t(glmnet.epsilons)))
				   
data.frame(colMeans(glmnet.Q.star), colMeans(rf.Q.star))

mean(Q)
mean(Q.star)


Q.star.by.epsilon<-sapply(1:length(epsilons),
                          function(x) plogis(qlogis(all.Q.by.hosp[,x]) + 
                                               (epsilons[x]* (1/modified.prob.of.exposure))))


											   
											   
save(epsilons, Q.star.by.epsilon, all.Q.by.hosp, file=output.file)

# Simple way to get crude:
# aggregate(day_30_readmit~hosp,disease.df,mean)
# head(disease.df)
