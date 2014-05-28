# /usr/bin/R --args data_dump/disease_ami.object data_dump/rf_G_model_ami.object ./matrix_cache
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12)

object.file<-commandArgs(trailingOnly=TRUE)[1]
output.file<-commandArgs(trailingOnly=TRUE)[2]
matrix.cache<-commandArgs(trailingOnly=TRUE)[3]
load(object.file)

disease.big.df<-as.data.frame(disease.big.matrix)
disease.big.df<-disease.big.df[!grepl('^hosp',colnames(disease.big.df))]
ntrees=1200

set.seed(1)
rf.predict.exposure<-bigrfc(x=disease.big.df, 
						  y=as.factor(disease.df$hosp), 
						  ntrees = ntrees, 
						  cachepath=matrix.cache,
						  trace = 0,
						  printclserr = FALSE)
						  
save(rf.predict.exposure,file=output.file)
