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
