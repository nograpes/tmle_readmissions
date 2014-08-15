suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gam))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.

prefixes<-c('ami','heart_failure','pneumonia')
pretty.names<-c('Acute myocardial infarction','Heart failure','Pneumonia')
names(pretty.names)<-prefixes

get.data<-function(prefix) {
  dump.dir<-'data_dump'
  data.files<-c('rf_Q_star_model_', 'disease_', 'glmnet_Q_model_')
  files<-paste0('data_dump/', data.files, prefix, '.object')
  e<-new.env()
  for(file in files) load(file,envir=e)
  # The calibrated models have the same names, unfortunately.
  calibrated.data.files<-c('rf_G_calibrated_model_', 'rf_Q_calibrated_model_')
  calibrated.files<-paste0('data_dump/', calibrated.data.files, prefix, '.object')
  e2<-new.env()
  for(calibrated.file in calibrated.files)
    load(calibrated.file, envir=e2)
  for(item in ls(e2)) assign(paste0('calibrated.',item), e2[[item]], envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.

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

# Exposure
rf.exposure.accuracy.by.disease <- sapply(models,function(x)
  get.accuracy(x[['calibrated.rf.predict.exposure']], x$disease.df$hosp))
colnames(rf.exposure.accuracy.by.disease) <-pretty.names[names(models)]

# Outcome
rf.outcome.accuracy.by.disease <- sapply(models,function(x)
  get.accuracy(x[['calibrated.rf.predict.outcome']], 1))
colnames(rf.outcome.accuracy.by.disease) <-  pretty.names[names(models)]

exposure.accuracy.df<-melt(rf.exposure.accuracy.by.disease,
                           varnames = c('trees', 'disease'),
                           value.name = 'accuracy')

outcome.accuracy.df<-melt(rf.outcome.accuracy.by.disease,
                          varnames = c('trees', 'disease'),
                          value.name = 'accuracy')




both.accuracy.df<-rbind(data.frame(exposure.accuracy.df,type='g model\n(Hospital choice)'),
                        data.frame(outcome.accuracy.df,type='Q model\n(Readmission)'))

p<-ggplot(both.accuracy.df,
          aes(x=trees,
              y=1-accuracy,
              col=disease)) +
  facet_grid(type~.) +
  geom_line(size=1.5) +
  labs(x='Number of trees',
       y='Error rate (out-of-bag)') +
  scale_x_continuous(breaks=seq(0,1200,by=300)) +
  scale_y_continuous(labels=c(0,1), breaks=c(0,1), limits=c(0,1)) +
  scale_colour_discrete(name = 'Admission diagnosis') +
  theme(legend.position = 'right',
        text=element_text(family="Cambria"))
ggsave(filename="figures/error_rate_g_and_Q.png",
       plot=p,width=(8.5 - (0.5*2)),height=11/3,dpi=300)
