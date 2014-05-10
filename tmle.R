suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12)
options(mc.cores=12)

dump.dir<-'data_dump'
setwd('~/repo/thesis/code/tmle')
object.file<-'heart.failure.object'
prefix<-gsub('object$','',object.file)
load(paste(dump.dir,object.file,sep='/'))

# I'm not sure why this isn't saved. I tested it, and it should have been.
hosps<-c(levels(disease.df$hosp))

# Time for some quick data reduction!
# If I don't do this, then the cross-validation for lambda becomes too much.
# 5000 variables are exactly zero!
disease.big.matrix<-disease.big.matrix[,colSums(disease.big.matrix)>30]
testo<-as.data.frame(disease.big.matrix)
testo<-testo[-match(hosps[-1],names(testo))]
testo<-testo[-grep('csd_ct_uid',names(testo))]
testo$csd_ct_uid<-as.factor(paste(disease.df$csd,disease.df$cma,disease.df$ct,sep='_'))

ntrees=1200

set.seed(1)
system.time(
  rf.predict.exposure<-bigrfc(x=testo, 
                              y=disease.df$hosp, 
                              ntrees = ntrees, 
                              cachepath='~',
                              trace = 1)
)

save(rf.predict.exposure,file=paste(dump.dir,'rf.predict.exposure.object',sep='/'))

set.seed(1)
system.time(
  rf.predict.outcome<-bigrfc(x=cbind(testo,disease.df$hosp), 
                             y=as.factor(disease.df$day_30_readmit), 
                             ntrees = ntrees, 
                             cachepath='~',
                             trace = 1)
)

save(rf.predict.outcome,file=paste(dump.dir,'rf.predict.outcome.object',sep='/'))

load(paste(dump.dir,'rf.predict.outcome.object',sep='/'))
load(paste(dump.dir,'rf.predict.exposure.object',sep='/'))

votes<-rf.predict.exposure@oobvotes
prob<-votes/rowSums(votes)
prob.of.exposure <- prob[cbind(1:length(disease.df$hosp),match(disease.df$hosp,colnames(prob)))]

votes<-rf.predict.outcome@oobvotes
prob<-votes/rowSums(votes)
prob.of.outcome <- prob[,"TRUE"]

# Just a couple of checks that the models aren't crazy.
# data.frame(prop.table(table(disease.df$hosp)),colMeans(prob))
sum(diag(rf.predict.exposure@trainconfusion)) / nrow(disease.df) # A 58% accuracy.


# Find majority vote by number of trees.
votes<-rf.predict.exposure@oobvotes

# A function that gives all oob votes by tree.
# A zero vote for in-bag.
get.oob.pred.by.tree<-function(tree) {
  in.bag<-tree@insamp!=0
  pred.class<-tree@trainpredclass
  pred.class[in.bag]<-0
  pred.class
}

oob.pred.mat<-sapply(rf.predict.exposure,get.oob.pred.by.tree)
#get.pred.by.row<-function(x) as.numeric(names(sort(table(x[x!=0]),decreasing=TRUE))[1])

# All right, so this one is complicated.
# What this generates is a list of 20, one for each hospital.
# Each element of the list is a matrix which is num.discharges by num.trees
# So, for element 1 of the list of 20 that is for hospital 1.
# So, the index 1,1 of that matrix is the number of votes for hospital 1 for the first patient when the first 1 tree(s) are included.
# Don't think about it too hard.
l<-lapply(1:20, function(x) t(apply(oob.pred.mat == x, 1, cumsum)))
extract.by.tree.num<-function(tree.num) sapply(l,function(x,tree.num) x[,tree.num],tree.num=tree.num)

get.accuracy<-function(num.trees){
  # s<-sapply(1:20,function(x) rowSums(oob.pred.mat[,1:num.trees] == x))
  s<-extract.by.tree.num(num.trees)
  set.seed(1) # Because tie breakage is random. This is how randomForest does it too.
  prop.table(table(max.col(s,ties.method='random')==as.numeric(disease.df$hosp)))['TRUE']
}

acc.by.trees<-sapply(1:ncol(l[[1]]),get.accuracy)
plot(1-acc.by.trees,type='l') # Second degree bend around 200.

# Now, let's take a look at the Q estimates.
# get.oob.pred.by.tree(rf.predict.outcome[[1]])

# get.oob.outcome.pred.by.tree<-function(tree) {
#   in.bag<-tree@insamp!=0
#   pred.class<-tree@trainpredclass
#   pred.class[in.bag]<-NA
#   pred.class
# }
# 
# oob.outcome.pred.mat<-sapply(rf.predict.outcome,get.oob.outcome.pred.by.tree)
# 

prob.by.hosp<-function(hosp){
  set.to.hosp<-cbind(testo,disease.df$hosp)
  set.to.hosp[,"disease.df$hosp"]<-
    factor(hosp,levels=levels(disease.df$hosp))
  p=predict(rf.predict.outcome,set.to.hosp)
  prop.table(p@testvotes,margin=1)[,"TRUE"]
}
all.Q.by.hosp<-sapply(as.character(levels(disease.df$hosp)),prob.by.hosp)
# The Q
colMeans(all.Q.by.hosp)


# Create a matrix of indicator variables.
ff<-~hosp-1
# contrasts(disease.df$hosp)
exposure.mat<-(model.matrix(ff,model.frame(ff, disease.df)))
colnames(exposure.mat)<-gsub('^hosp','',colnames(exposure.mat))

# Q is the prob.of.outcome
# G is the prob.of.exposure

# Temporary fix:
# For seven people, the probability of exposure (going to the hospital the went to) is exactly zero.
# I'm going to bump it up slightly as a quick fix.
# But this should be subject to a lot of analysis.
min.baseline<-min(prob.of.exposure[prob.of.exposure!=0])
modified.prob.of.exposure<-ifelse(prob.of.exposure==0,min.baseline,prob.of.exposure)

# I want to find one eta for each hospital.
iptw <- exposure.mat  / modified.prob.of.exposure # Each column is now the (A==a)/g

eta<-function(hospital.covariate)
  glm(disease.df$day_30_readmit ~ 
        hospital.covariate-1, 
      offset=qlogis(prob.of.outcome), 
      family=binomial(link=logit))$coef

etas<-apply(iptw,2,eta) # The fitter complains that fitted probs of 0 or 1 occurred.
Q.star.by.eta<-sapply(1:length(etas),function(x) plogis(qlogis(prob.of.outcome) + (etas[x]*iptw[,x])))


res<-colMeans(Q.star.by.eta)
names(res)<-names(etas)


# The crude probability of readmission.
crude.prob<-prop.table(table(disease.df[c('hosp','day_30_readmit')]),margin=1)

# The crude probability of readmission.
data.frame(crude=crude.prob[,2], Q=colMeans(all.Q.by.hosp), `Q.star`=res)