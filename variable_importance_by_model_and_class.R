suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(ggplot2))

prefixes <- list.files('disease_subsets')
pretty.names<-c('AMI','Heart failure','Pneumonia')
names(pretty.names)<-prefixes

get.data<-function(prefix) {
  data.files<-c('rf_G_calibrated_model_', 'rf_Q_calibrated_model_','disease_')
  files<-paste0('data_dump/', data.files, prefix, '.object')
  e<-new.env()
  for(file in files) load(file,envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.

# Now let's make a nice density of the gini coefficient by variable type.
get.gini.types<-function(disease,model=c('G','Q')) {
  if(model=='G') rf.predict<-models[[disease]][['rf.predict.exposure']]
  if(model=='Q') rf.predict<-models[[disease]][['rf.predict.outcome']]
  
  disease.df<-models[[disease]][['disease.df']]
  disease.big.matrix<-models[[disease]][['disease.big.matrix']]
  var.names<-models[[disease]][['var.names']]
  
  gini.importance<-fastimp(rf.predict)
  gini.df<-data.frame(var=names(gini.importance),gini=gini.importance)
  
  types<-unname(rep(names(var.names),times=sapply(var.names,length)))
  vars<-unname(unlist(var.names))
  
  cts<-as.character(gini.df$var[grep('^csd_ct_uid',gini.df$var)])
  types<-c(types,rep('ct',length(cts)))
  vars<-c(vars,cts)
  var.types<-data.frame(var=vars,type=types)
  cbind(disease=unname(pretty.names[disease]),merge(gini.df,var.types))
}

G.gini.types<-do.call(rbind,lapply(prefixes,get.gini.types,model='G'))
Q.gini.types<-do.call(rbind,lapply(prefixes,get.gini.types,model='Q'))

gini.types<-rbind(cbind(model='g (Outcome - hospital choice)',G.gini.types),cbind(model='Q (Outcome - hospital readmission)',Q.gini.types))

p=ggplot(gini.types, aes(x=log(gini+exp(-12)),fill=type, order=type)) +
  facet_grid(disease~model) +
  geom_density(alpha=0.5) +
  theme(legend.position = 'bottom',text=element_text(family="Cambria")) +
  ylim(c(0,0.3)) +
  scale_fill_discrete('Variable class',labels=c(ct='Census tract',diagnosis='Diagnosis',drug='Drug',procedure='Procedure')) +
  labs(x=expression(paste("Variable importance (",log(gini + e^-12),")",sep=""))  ,
       y='Density of variables (within class)') 
ggsave(filename="figures/variable_importance_by_model_and_class.png", 
       plot=p,width=8.5 - (0.5*2),height=11*2/3,dpi=600)