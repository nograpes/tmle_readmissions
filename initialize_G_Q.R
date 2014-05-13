suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12)

object.file<-commandArgs(trailingOnly=TRUE)[1]
output.file<-commandArgs(trailingOnly=TRUE)[2]
load(object.file)

file.name.sub<-function(target,pattern,sub) 
   file.path(dirname(target),sub(pattern,sub,basename(target)))

# I'm not sure why this isn't saved. I tested it, and it should have been.
hosps<-c(levels(disease.df$hosp))

# Time for some quick data reduction!
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
                              trace = FALSE)
)


save(rf.predict.exposure,file=file.name.sub(object.file,'disease_','rf_predict_exposure_'))

set.seed(1)
system.time(
  rf.predict.outcome<-bigrfc(x=cbind(testo,disease.df$hosp), 
                             y=as.factor(disease.df$day_30_readmit), 
                             ntrees = ntrees, 
                             cachepath='~',
                             trace = FALSE)
)

save(rf.predict.exposure,file=file.name.sub(object.file,'disease_','rf_predict_outcome_'))
