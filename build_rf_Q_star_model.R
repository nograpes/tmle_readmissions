suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a serial backend to avoid warning.

# /usr/bin/R --args data_dump/rf_G_model_heart_failure.object data_dump/rf_Q_model_heart_failure.object data_dump/disease_heart_failure.object data_dump/rf_Q_star_model_heart_failure.object
# /usr/bin/R --args data_dump/rf_G_model_ami.object data_dump/rf_Q_model_ami.object data_dump/disease_ami.object data_dump/rf_Q_star_model_ami.object
# setwd('~/repo/thesis/code/tmle')
# arguments<-paste('data_dump',c('rf_G_model_pneumonia.object','rf_Q_model_pneumonia.object','disease_pneumonia.object','rf_Q_star_model_pneumonia.object'),sep='/')
arguments<-commandArgs(trailingOnly=TRUE)

G.model.file<-arguments[1]
Q.model.file<-arguments[2]
object.file<-arguments[3]
output.file<-arguments[4]
load(object.file)
load(G.model.file)
load(Q.model.file)

votes<-rf.predict.exposure@oobvotes
prob<-votes/rowSums(votes)
prob.of.exposure <- prob[cbind(1:length(disease.df$hosp),match(disease.df$hosp,colnames(prob)))]

votes<-rf.predict.outcome@oobvotes
prob<-votes/rowSums(votes)
prob.of.outcome <- prob[,"TRUE"]

prob.by.hosp<-function(hosp){
  set.to.hosp<-disease.big.matrix
  # Fix A to the hospital.
  set.to.hosp[,grep('^hosp',colnames(disease.big.matrix))]<-0
  var.name<-paste0('hosp',hosp)
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name]<-1
  
  predict.by.tree<-function(tree.num,forest,data,cachepath=getwd(),n){
    xtype <- as.integer(.Call("CGetType", data@address, PACKAGE = "bigmemory"))
    tree<-forest[[tree.num]]
    .Call('treepredictC',data@address,xtype,n,forest,tree,PACKAGE = "bigrf")$testpredclass 
  }
  data <- bigrf:::makex(set.to.hosp, "xtest", cachepath=getwd())
  # Use foreach instead of 
  prediction.by.tree<-foreach(i=seq_along(rf.predict.outcome),.combine=cbind) %dopar% predict.by.tree(tree=i,forest=rf.predict.outcome,data=data,n=nrow(set.to.hosp))
  prediction.by.tree<-prediction.by.tree == 2 # Because class 2 is true
  # Now set all in-bag to NA
  in.bag.mat<-sapply(rf.predict.outcome,function(tree)tree@insamp!=0)
  prediction.by.tree[in.bag.mat]<-NA
  (rowSums(prediction.by.tree,na.rm=TRUE) / rowSums(!in.bag.mat)[1])
}

# Predict by all trees.
all.Q.by.hosp<-sapply(levels(disease.df$hosp),prob.by.hosp)

# Create a matrix of indicator variables.
ff<-~hosp-1
# contrasts(disease.df$hosp)
exposure.mat<-(model.matrix(ff,model.frame(ff, disease.df)))
colnames(exposure.mat)<-gsub('^hosp','',colnames(exposure.mat))

# Temporary fix:
# For seven people, the probability of exposure (going to the hospital the went to) is exactly zero.
# I'm going to bump it up slightly as a quick fix.
# But this should be subject to a lot of analysis.
min.baseline<-min(prob.of.exposure[prob.of.exposure!=0])
modified.prob.of.exposure<-ifelse(prob.of.exposure==0,min.baseline,prob.of.exposure)

# The same for the outcome
min.baseline<-min(prob.of.outcome[prob.of.outcome!=0])
modified.prob.of.outcome<-ifelse(prob.of.outcome==0,min.baseline,prob.of.outcome)

# I want to find one epsilon for each hospital.
iptw <- exposure.mat  / modified.prob.of.exposure # Each column is now the (A==a)/g

epsilon<-function(hospital.covariate)
  glm(as.factor(disease.df$day_30_readmit) ~ 
        hospital.covariate-1, 
      offset=qlogis(modified.prob.of.outcome), 
      family=binomial(link=logit))$coef

 
epsilons<-apply(iptw,2,epsilon) # The fitter complains that fitted probs of 0 or 1 occurred.

Q.star.by.epsilon<-sapply(1:length(epsilons),
                          function(x) plogis(qlogis(all.Q.by.hosp[,x]) + 
                                               (epsilons[x]* (1/modified.prob.of.exposure))))

save(epsilons, Q.star.by.epsilon, all.Q.by.hosp, file=output.file)



