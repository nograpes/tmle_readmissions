suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doMC))
registerDoMC(cores=10) # For 10 folds

# arguments<-c('data_dump/disease_ami.object', 'survival/data_dump/glmnet_Q_ami.object')
arguments<-commandArgs(trailingOnly=TRUE)
object.file<-arguments[1]
output.file<-arguments[2]
load(object.file)

set.seed(1)
glmnet.predict.censor <- 
              cv.glmnet(x=disease.big.matrix, 
                     y=Surv(disease.df$tte, 
                            disease.df$censor),
					 nfolds=10,
                     family='cox',
 					 parallel=TRUE)
save(glmnet.predict.censor,file=output.file)