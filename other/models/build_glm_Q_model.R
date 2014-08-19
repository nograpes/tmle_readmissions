# /usr/bin/R --args data_dump/disease_ami.object data_dump/glm_Q_model_ami.object 
object.file<-commandArgs(trailingOnly=TRUE)[1]
output.file<-commandArgs(trailingOnly=TRUE)[2]
load(object.file)

glm.predict.outcome<-glm.fit(x=disease.big.matrix, 
                             y=disease.df$day_30_readmit,
                             family=binomial(link=logit))

save(glm.predict.outcome,file=output.file)