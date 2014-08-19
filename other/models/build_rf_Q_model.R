# /usr/bin/Rscript --vanilla ./build_rf_Q_model.R data_dump/disease_ami.object data_dump/rf_Q_model_ami.object ./matrix_cache
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12)

object.file<-commandArgs(trailingOnly=TRUE)[1]
output.file<-commandArgs(trailingOnly=TRUE)[2]
matrix.cache<-commandArgs(trailingOnly=TRUE)[3]
load(object.file)

# Time for some quick data reduction!
ntrees=1200

set.seed(1)
rf.predict.outcome<-bigrfc(x=as.data.frame(disease.big.matrix), 
						 y=as.factor(disease.df$day_30_readmit), 
						 ntrees = ntrees, 
						 cachepath=matrix.cache,
						 trace = 0,
						 printclserr = FALSE)

save(rf.predict.outcome,file=output.file)