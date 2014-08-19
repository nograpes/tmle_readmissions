# /usr/bin/R --args data_dump/disease_ami.object data_dump/glmnet_Q_model_ami.object 
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
registerDoMC(cores=10) # For 10 folds

object.file<-commandArgs(trailingOnly=TRUE)[1]
output.file<-commandArgs(trailingOnly=TRUE)[2]
load(object.file)

set.seed(1)
glmnet.predict.outcome<-cv.glmnet(x=disease.big.matrix, 
                             y=disease.df$day_30_readmit,
                             family='binomial',
                             type.logistic='modified.Newton',
                             parallel=TRUE
                            )

save(glmnet.predict.outcome,file=output.file)