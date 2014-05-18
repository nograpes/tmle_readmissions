# /usr/bin/R --args ami heart_failure pneumonia
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))

# Read this in programatically.
prefixes<-c('ami','heart_failure')
pretty.names<-c('Acute myocardial infarction','Heart failure')
names(pretty.names)<-prefixes
setwd('~/repo/thesis/code/tmle')

models<-list()
for (prefix in prefixes){
  e<-new.env()
  load(paste0('data_dump/rf_G_model_',prefix,'.object'),envir=e) 
  load(paste0('data_dump/rf_Q_star_model_',prefix,'.object'),envir=e)
  load(paste0('data_dump/disease_',prefix,'.object'),envir=e) 
  models[[prefix]]<-mget(ls(e),envir=e)
}
rm(e)

# A zero vote for in-bag.
get.oob.pred.by.tree<-function(tree) {
  in.bag<-tree@insamp!=0
  pred.class<-tree@trainpredclass
  pred.class[in.bag]<-0
  pred.class
}

get.accuracy<-function(rf.predict.exposure, truth){
  class.num<-length(rf.predict.exposure@ytable)
  oob.pred.mat<-sapply(rf.predict.exposure, get.oob.pred.by.tree)
  l<-lapply(1:class.num, function(x) t(apply(oob.pred.mat == x, 1, cumsum)))
  extract.by.tree.num<-function(tree.num) {
    s<-sapply(l,function(x,tree.num) x[,tree.num],tree.num=tree.num)
    prop.table(table(max.col(s,ties.method='random')==as.numeric(truth)))['TRUE']
  }
  unname(sapply(seq_along(rf.predict.exposure),extract.by.tree.num))
}

rf.exposure.accuracy.by.disease<-sapply(models,function(x) with(x,get.accuracy(rf.predict.exposure, disease.df$hosp)))

colnames(rf.exposure.accuracy.by.disease)<-pretty.names[colnames(rf.exposure.accuracy.by.disease)]

accuracy.df<-melt(rf.exposure.accuracy.by.disease,varnames = c('trees','disease'),value.name='accuracy')

theme_set(theme_gray(base_size = 16))
ggplot(accuracy.df,
       aes(x=trees,y=1-accuracy,col=disease)) + 
  geom_line(size=1.5) +
  labs(title='Error rate for random forest model of hospital choice',
       x='Number of trees',
       y='Error rate (out-of-bag)') +
  scale_colour_discrete(name = 'Admission diagnosis') +
  theme(legend.position = 'bottom')

# Also, let's make a table of the variable importance.
help(paste0(class(rf.predict.exposure),'-class'))

# Now, let's take a look at the Q estimates.
rf.predict.exposure<-models[['ami']][['rf.predict.exposure']]
disease.df<-models[['ami']][['disease.df']]
disease.big.matrix<-models[['ami']][['disease.big.matrix']]
t(t(head(rev(sort(fastimp(rf.predict.exposure))),25)))
table(disease.df$day_30_readmit,disease.big.matrix[,'1893'])

# 189.3 - Malignant neoplasm of urethra
addmargins(table(disease.df$hosp,disease.big.matrix[,'1893']))





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


#get.pred.by.row<-function(x) as.numeric(names(sort(table(x[x!=0]),decreasing=TRUE))[1])

# All right, so this one is complicated.
# What this generates is a list of 20, one for each hospital.
# Each element of the list is a matrix which is num.discharges by num.trees
# So, for element 1 of the list of 20 that is for hospital 1.
# So, the index 1,1 of that matrix is the number of votes for hospital 1 for the first patient when the first 1 tree(s) are included.
# Don't think about it too hard.

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