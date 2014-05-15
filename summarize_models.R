# /usr/bin/R --args ami heart_failure pneumonia
suppressPackageStartupMessages(library(bigrf))

# Just a couple of checks that the models aren't crazy.
# data.frame(prop.table(table(disease.df$hosp)),colMeans(prob))
# sum(diag(rf.predict.exposure@trainconfusion)) / nrow(disease.df) # A 58% accuracy.

# Read this in programatically.
prefix<-'heart_failure'
setwd('~/repo/thesis/code/tmle')
load(paste0('data_dump/rf_Q_star_model_',prefix,'.object')) # epsilons, Q.star.by.epsilon
load(paste0('data_dump/disease_',prefix,'.object')) # epsilons, Q.star.by.epsilon

# A little plot for the random forest:
get.accuracy<-function(num.trees){
  # s<-sapply(1:20,function(x) rowSums(oob.pred.mat[,1:num.trees] == x))
  s<-extract.by.tree.num(num.trees)
  set.seed(1) # Because tie breakage is random. This is how randomForest does it too.
  prop.table(table(max.col(s,ties.method='random')==as.numeric(disease.df$hosp)))['TRUE']
}

acc.by.trees<-sapply(1:ncol(l[[1]]),get.accuracy)
# plot(1-acc.by.trees,type='l') # Second degree bend around 200.

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



crude.prob<-prop.table(table(disease.df[c('hosp','day_30_readmit')]),margin=1)
comparison<-data.frame(crude=crude.prob[,2], Q=colMeans(all.Q.by.hosp), `Q.star`=colMeans(Q.star.by.epsilon))

round(comparison,3)

epsilons
Q.star.by.epsilon

dim(Q.star.by.epsilon)

colMeans(Q.star.by.epsilon)
epsilons

res<-colMeans(Q.star.by.eta)
names(res)<-names(etas)

# The crude probability of readmission.


# What about logistic regression?
# Time for some quick data reduction!
hosps<-c(levels(disease.df$hosp))

disease.big.matrix<-disease.big.matrix[,colSums(disease.big.matrix)>30]
testo<-disease.big.matrix
csd_ct_uid.df<-data.frame(csd_ct_uid=as.factor(paste(disease.df$csd,disease.df$cma,disease.df$ct,sep='_')))

ff<-~csd_ct_uid-1
# contrasts(disease.df$hosp)
ct.mat<-(model.matrix(ff,model.frame(ff, csd_ct_uid.df)))
testo<-cbind(testo,ct.mat)
library(glmnet)

system.time(g<-glmnet(x=testo,
       y=disease.df[,'day_30_readmit'],
       family='binomial'))


testo<-testo[-match(hosps[-1],names(testo))] # This removes the hospitals?
testo<-testo[-grep('csd_ct_uid',names(testo))]
testo$csd_ct_uid<-as.factor(paste(disease.df$csd,disease.df$cma,disease.df$ct,sep='_'))


round(prop.table(table(data$csd_ct_uid,data$day_30_readmit),margin=1),2)

data=data.frame(day_30_readmit=as.factor(disease.df$day_30_readmit),testo,hosp=disease.df$hosp)


glm.predict.outcome<-
  glm(day_30_readmit~.,
      family=binomial(link=logit),
      data=data
  )


confint(glm.predict.outcome)













# Accuracy by number of trees.
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
# End plot accuracy by number of trees.


# res<-colMeans(Q.star.by.epsilon)
# names(res)<-names(epsilons)


# The crude probability of readmission.
# crude.prob<-prop.table(table(disease.df[c('hosp','day_30_readmit')]),margin=1)

# The crude probability of readmission.
# data.frame(crude=crude.prob[,2], Q=colMeans(all.Q.by.hosp), `Q.star`=res)