suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(gam))
registerDoParallel(cores=12) # Register a parallel backend -- prediction is slow.

# Read this in programatically.
prefixes<-c('ami','heart_failure','pneumonia')
pretty.names<-c('Acute myocardial infarction','Heart failure','Pneumonia')
names(pretty.names)<-prefixes
setwd('~/repo/thesis/code/tmle')
get.data<-function(prefix) {
  dump.dir<-'data_dump'
  data.files<-c('rf_Q_star_model_', 'disease_','crude_readmissions_risk_')
  files<-paste0('data_dump/', data.files, prefix, '.object')
  files<-c(files,paste0('survival/data_dump/Q_star_survival_',prefix,'.object'))
  e<-new.env()
  for(file in files) load(file,envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.

# Now, make a table for Q*?
dump.base.stats<-function(disease,los=FALSE) {
  disease.df <- models[[disease]][['disease.df']]
  crude.readmitted.n<-table(disease.df$hosp,disease.df$day_30_readmit)[,'TRUE']
  n<-c(table(disease.df$hosp))
  crude.props<-prop.table(table(disease.df$hosp,disease.df$day_30_readmit),margin=1)[,'TRUE']
  # Died during stay.
  died.during.stay <- models[[disease]] [['died.during.stay']]
  n.died <- c(table(died.during.stay$hosp))
  # Length-of-stay
  
  dash.NAs <- function(x){
    m<-merge(data.frame(hosp=levels(disease.df$hosp)),x,all.x=TRUE)
    # m[is.na(m$los),'los']<-'-'
    m
  }

  if (los) {
  mean.survived.los <- dash.NAs(aggregate(los~hosp,disease.df,mean))[,'los']
  mean.died.los <- dash.NAs(aggregate(los~hosp,died.during.stay,mean))[,'los']
  mean.both.los <- dash.NAs(aggregate(los~hosp,
                                      rbind(disease.df[c('hosp','los')], 
                                            died.during.stay[c('hosp','los')]),
                                      mean)) [,'los']
  }
  
  d = data.frame(admitted=n+n.died, died=n.died, died.prop=n.died/(n+n.died), live.discharge=n, readmitted=crude.readmitted.n, prop=crude.props)
  if(los) return(cbind(d,data.frame(overall.los=mean.both.los,died.los=mean.died.los,survived.los=mean.survived.los)))
  d
}

dump.results.df<-function(disease, model) {
  # all.Q.by.hosp <- models[[disease]][[paste0('all.',model,'.Q.by.hosp')]]
  # epsilon <- models[[disease]][[paste0(model,'.epsilons')]]
  Q.star.by.epsilon <- models[[disease]][[paste0(model,'.Q.star')]]
  tte <- models[[disease]][['tte.mat']]
  
  hazards.mat <- models[[disease]][['hazards.mat']]
  logistic.mat <- models[[disease]][['logistic.mat']]
  base <-  models[[disease]]$base
  
  add.base <- function (x) 
    rbind(x[1:base-1,],c(NA,NA,NA),x[base:nrow(x),])
  
  odds <- add.base(logistic.mat)
  colnames(odds) <- c('odds.ratio','odds.ratio.ci.low','odds.ratio.ci.high')

  hazards <- add.base(hazards.mat)
  colnames(hazards) <- c('hazard.ratio','hazard.ratio.ci.low','hazard.ratio.ci.high')
  
  data.frame(odds,Q.star=colMeans(Q.star.by.epsilon),
             hazards,mean.tte=colMeans(tte),median.tte=apply(tte,2,median))
}

disease.table<-function(disease, models='rf') {
  model.tables <- lapply(models, dump.results.df, disease=disease)
  cbind(dump.base.stats(disease), Reduce(cbind, model.tables))
}

disease.results.table<-sapply(prefixes, disease.table, simplify=FALSE)
save(disease.results.table,file='tables/disease.results.table.object')