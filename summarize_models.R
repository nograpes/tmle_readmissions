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

colnames(rf.exposure.accuracy.by.disease) <-
  pretty.names[names(models)]

# Outcome
rf.outcome.accuracy.by.disease <- sapply(models,function(x)
  get.accuracy(x[['calibrated.rf.predict.outcome']], 1))

colnames(rf.outcome.accuracy.by.disease) <-
  pretty.names[names(models)]

exposure.accuracy.df<-melt(rf.exposure.accuracy.by.disease,
                  varnames = c('trees', 'disease'),
                  value.name = 'accuracy')

outcome.accuracy.df<-melt(rf.outcome.accuracy.by.disease,
                  varnames = c('trees', 'disease'),
                  value.name = 'accuracy')


both.accuracy.df<-rbind(data.frame(exposure.accuracy.df,type='g model (Exposure - Hospital choice)'),
                        data.frame(outcome.accuracy.df,type='Q model (Outcome - Readmission)'))

p<-ggplot(both.accuracy.df,
          aes(x=trees,
              y=1-accuracy,
              col=disease)) +
  facet_grid(type~.) +
  geom_line(size=1.5) +
  labs(x='Number of trees',
       y='Error rate (out-of-bag)') +
  scale_colour_discrete(name = 'Admission diagnosis') +
  theme(legend.position = 'right',
          text=element_text(family="Cambria"))
ggsave(filename="figures/error_rate_for_hospital_choice.png",
       plot=p,width=(8.5 - (0.5*2)),height=11/3,dpi=300)

# What is the best sensitivty specificity?

# Variable importance.
get.top.gini<-function(disease,top=10) {
  rf.predict.exposure<-models[[disease]][['rf.predict.exposure']]
  disease.df<-models[[disease]][['disease.df']]
  disease.big.matrix<-models[[disease]][['disease.big.matrix']]

  top.10<-t(t(head(rev(sort(fastimp(rf.predict.exposure))),top)))
  names.top.10<-rownames(top.10)

  convert<-c(proc='Procedure',diag='Diagnosis',drug='Drug')
  pretty.names.func<-function(x){
    strip<-function(y) {
      name<-sub('^[^_]*_','',y)
      # paste0(name,' (',type,')')
      name
    }
    if (is.factor(x)) factor(strip(x),levels=strip(levels(x)))
    else strip(x)
  }

  heat.map<-t(sapply(names.top.10, function(x) prop.table(table(disease.df$hosp,disease.big.matrix[,x])[,'1'])))

  melted.heat.map<-melt(heat.map,varnames=c('raw_variable','Hospital'))
  melted.heat.map$Hospital<-factor(as.character(as.numeric(melted.heat.map$Hospital)),levels=as.character(1:20))
  melted.heat.map$raw_variable<-factor(melted.heat.map$raw_variable,levels=rev(names.top.10))
  melted.heat.map$gini<-top.10[melted.heat.map$raw_variable]
  melted.heat.map$Variable<-pretty.names.func(melted.heat.map$raw_variable)
  melted.heat.map$disease<-pretty.names[disease]
  melted.heat.map$type<-convert[sapply(lapply(strsplit(as.character(melted.heat.map$raw_variable),''),`[`,1:4),paste0,collapse='')]
  melted.heat.map
}

top.gini.list<-lapply(prefixes,get.top.gini)
top.gini<-do.call(rbind,top.gini.list)

top.gini$Variable<-as.character(top.gini$Variable)
top.gini$var.dis<-reorder(paste(top.gini$Variable,top.gini$disease),top.gini$gini)

name.list<-as.character(top.gini$Variable)
names(name.list)<-as.character(top.gini$var.dis)

name.list[grep('Scinti',name.list)]<-'Scintigraphie et étude de fonctionnement cardiovasculaire'

line_breaker<-function(x,after=20)
  sub(paste0('(.{',after,'}[^ ]*) '),'\\1\n',x)

p <- ggplot(top.gini,
            aes(x=Hospital, y=var.dis, fill=value, group=disease, label=Variable)) +
  geom_tile(colour = "white") +
  facet_grid(disease~., scale="free_y") +
  scale_y_discrete('Variable',labels=line_breaker(name.list)) +
  scale_fill_gradient(name='Proportion choosing hospital',
                      low = "lightblue", high = "steelblue") +
  theme(legend.position = 'bottom',
        text=element_text(family="Cambria"))
# ,axis.text.y = element_text(colour = top.gini$type)
ggsave(filename="figures/top_10_variable_importance_and_hospital.png",
       plot=p,width=8.5 - (0.5*2),height=11 - (0.5 * 2) -0.5,dpi=300)



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

gini.types<-rbind(cbind(model='G (Outcome - hospital choice)',G.gini.types),cbind(model='Q (Outcome - hospital readmission)',Q.gini.types))

p=ggplot(gini.types, aes(x=log(gini+exp(-12)),fill=type, order=type)) +
  facet_grid(disease~model) +
  geom_density(alpha=0.5) +
  theme(legend.position = 'bottom',text=element_text(family="Cambria")) +
  ylim(c(0,0.4)) +
  scale_fill_discrete('Variable class',labels=c(ct='Census tract',diagnosis='Diagnosis',drug='Drug',procedure='Procedure')) +
  labs(x=expression(paste("Variable importance (",log(gini + e^-12),")",sep=""))  ,
       y='Density of variables (within class)')
ggsave(filename="figures/variable_importance_by_model_and_class.png",
       plot=p,width=8.5 - (0.5*2),height=11*2/3,dpi=300)



# Hospital choropleth.
# Counts by CT?
suppressPackageStartupMessages(library(RPostgreSQL))
suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library(maptools))
suppressPackageStartupMessages(library(scales))

# Pull in Montreal shapefile
montreal<-readShapePoly('../maps/ct/Census_Tract_2006')
names(montreal@data) [names(montreal@data) == 'CTUID'] <- 'ct'
test<-fortify(montreal,region='ct')
test$x<-test$long # necessary
test$y<-test$lat # necessary
# Calculate bounding box.
xlim<-montreal@bbox[1,] + (c(0.25,-0.25)*diff(montreal@bbox[1,]))
ylim<-montreal@bbox[2,] + (c(0.25,-0.25)*diff(montreal@bbox[2,]))

# Get counts by hosp for each disease
melted.ct.by.hosp.by.disease<-function(disease) {
  disease.df<-models[[disease]][['disease.df']]
  ct.by.hosp<-with(disease.df[,c('cma','ct','hosp')],table(paste0(cma,ct),hosp))
  melt(ct.by.hosp, varnames=c('ct','hosp'), value.name='count', as.is=TRUE)
}
melted.ct.by.hosp<-sapply(prefixes, melted.ct.by.hosp.by.disease, simplify=FALSE)

# Grab hospital location data from etude 2
drv <- dbDriver("PostgreSQL")
con.etude2 <- dbConnect(drv, dbname="new_etude2")
hosps<-dbGetQuery(con.etude2,'select y.hosp,x.* from hosp_locs x left join hosp_map y on x.name=y.name',encoding='UTF-8')

# Merge on the population:
# Connect to geo to add in the population information.
con.geo <- dbConnect(drv, dbname="geo")
ct.pop<-dbGetQuery(con.geo,'select geocode,col2 as pop2006 from census.census2006_first1200')
montreal@data$geocode<-gsub('\\.','',montreal@data$ct)
montreal@data<-merge(montreal@data,ct.pop,all.x=TRUE)

normalized.counts<-function(melted.ct.by.hosp,hosp){
  counts.by.ct<-melted.ct.by.hosp[melted.ct.by.hosp$hosp==hosp,c('ct','count')]
  test.values<-merge(montreal@data, counts.by.ct, all.x=TRUE)
  test.values$normalized.by.pop<-test.values$count/test.values$pop2006
  test.values$count[is.na(test.values$count)]<-0
  test.values$normalized.by.pop[is.na(test.values$normalized.by.pop)]<-0
  test.values
}

normalized.counts.disease<-lapply(melted.ct.by.hosp, normalized.counts, hosp='Hôpital Charles Lemoyne')
normalized.counts.disease<-lapply(names(normalized.counts.disease),function(x) data.frame(disease=unname(pretty.names[x]),normalized.counts.disease[[x]],stringsAsFactors=FALSE))


map.theme<-theme(axis.text.y = element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(),axis.ticks=element_blank(),legend.position='left')


normalized.counts.disease.bound<-do.call(rbind,normalized.counts.disease)

p=ggplot(data=normalized.counts.disease.bound, aes(fill=normalized.by.pop*1000)) +
  geom_map(aes(map_id=ct), map=test, size=0.0375, colour='white') +
  expand_limits(test) + xlim(xlim) + ylim(ylim) +
  facet_grid(~disease) +
  scale_fill_continuous(name='Rate of seeking care at a selected hospital\nby disease per thousand residents',high = "#132B43",low = "#56B1F7", guide = guide_colorbar(barwidth=12)) +
  map.theme +
  theme(legend.position = 'bottom',text=element_text(family="Cambria"))

ggsave(filename="figures/hosp_choro.png", plot=p,dpi=700,
       width=8.5 - (0.5*2),
       height=((8.5 - (0.5*2))/3) + 1)


# Just to be clear... do the Q models suck?
# We can measure this by calibration.
get.accuracy.outcome<-function(rf.predict.exposure, truth){
  class.num<-length(rf.predict.exposure@ytable)
  oob.pred.mat<-sapply(rf.predict.exposure, get.oob.pred.by.tree)
  l<-lapply(1:class.num, function(x) t(apply(oob.pred.mat == x, 1, cumsum)))
  extract.by.tree.num<-function(tree.num) {
    s<-sapply(l,function(x,tree.num) x[,tree.num],tree.num=tree.num)
    prop.table(table(max.col(s,ties.method='random')==as.numeric(as.factor(truth))))['TRUE']
  }
  unname(sapply(seq_along(rf.predict.exposure),extract.by.tree.num))
}

extract.by.tree.num<-function(tree.num) {
  s<-sapply(l,function(x,tree.num) x[,tree.num],tree.num=tree.num)
  prop.table(table(max.col(s,ties.method='random')==as.numeric(as.factor(truth))))['TRUE']
}


get.oob.predictions<-function(disease){
  e=new.env()
  load(paste0('data_dump/rf_Q_calibrated_model_',disease,'.object'), envir=e)
  calibrated.rf<-get('rf.predict.outcome', e)
  disease.df<-models[[disease]][['disease.df']]

  class.num<-length(calibrated.rf@ytable)
  oob.pred.mat<-sapply(calibrated.rf, get.oob.pred.by.tree)
  l<-lapply(1:class.num, function(x) t(apply(oob.pred.mat == x, 1, cumsum)))
  tree.num=1200
  s<-sapply(l,function(x,tree.num) x[,tree.num],tree.num=tree.num)
  prob.of.outcome<-prop.table(s,margin=1)[,2]

  data.frame(disease=disease, prob=prob.of.outcome, truth=as.numeric(disease.df$day_30_readmit))
}

# Perhaps a calibration fix?
ami.preds<-get.oob.predictions('ami')
g.logit<-glm(truth~prob,data=ami.preds,family=binomial(link=logit))
g.linear<-glm(truth~prob,data=ami.preds)
g.spline.logit<-glm(truth~bs(prob),data=ami.preds,family=binomial(link=logit))

ami.preds$linear.smooth<-predict(g.linear)
ami.preds$logit.smooth<-predict(g.logit,type='response') # Basically, multiply by 2/3
ami.preds$spline.logit<-predict(g.spline.logit,type='response')

sum((g.linear$residuals)^2)
sum((plogis(g.logit$residuals))^2)

x=seq(0,1,by=0.001)
plot(type='l',x,predict(g.logit, newdata=data.frame(prob=seq(0,1,by=0.001)),type='response'),col='red',lwd=2)
lines(x,predict(g.spline.logit, newdata=data.frame(prob=x),type='response'),type='l',col='blue',lwd=2)

range(ami.preds$logit.smooth)
range(ami.preds$spline.logit)

par(mfrow=c(2,1))
hist(ami.preds$spline.logit)
hist(ami.preds$logit.smooth)
par(mfrow=c(1,1))

?geom_density

ggplot(m[m$variable %in% c('linear.smooth','prob','logit.smooth','spline.logit'),]) + geom_density(aes(x=value,fill=variable),alpha=0.5)

range(ami.preds$spline.logit)

m<-melt(ami.preds,id.vars=c('disease', 'truth'))
head(m)

ggplot(m[m$variable %in% c('logit.smooth','spline.logit'),], aes(x=value,y=truth,group=variable,colour=variable)) +
  geom_smooth(method=gam,size=1.25) +
  geom_abline(intercept = 0, slope = 1, colour='black')

head(ami.preds,10)

predictions.by.disease<-do.call(rbind,lapply(prefixes,get.oob.predictions))

p<-ggplot(predictions.by.disease,aes(x=prob, y=truth, group=disease, colour=disease)) +
  geom_smooth(method=gam,size=1.25) +
  geom_abline(intercept = 0, slope = 1, colour='black') +
  scale_x_continuous('Out-of-bag predicted probability of 30-day readmission',expand=c(0,0),limits=c(0,1)) +
  scale_y_continuous('GAM smoothed actual outcome',expand=c(0,0),limits=c(0,1)) +
  scale_colour_discrete('Admission diagnosis', labels=pretty.names) +
  theme(legend.position = 'bottom',text=element_text(family="Cambria"))
ggsave(filename="figures/rf_calibration.png", plot=p,dpi=700,
       width=(8.5) - (0.5*2),
       height=((8.5 - (0.5*2))/3) + 1)




# What about the GLM?
get.glmnet.predictions<-function(disease){
  e=new.env()
  load(paste0('data_dump/glmnet_Q_model_',disease,'.object'), envir=e)
  glm.predict.outcome<-get('glmnet.predict.outcome', e)
  disease.df<-models[[disease]][['disease.df']]
  disease.big.matrix<-models[[disease]][['disease.big.matrix']]

  prob.of.outcome<-c(plogis(predict(glm.predict.outcome, newx=disease.big.matrix, s=glm.predict.outcome$lambda.1se)))

  data.frame(disease=disease, prob=prob.of.outcome, truth=as.numeric(disease.df$day_30_readmit))
}

glmnet.predictions.by.disease<-do.call(rbind,lapply(prefixes,get.glmnet.predictions))

p<-ggplot(glmnet.predictions.by.disease, aes(x=prob,y=truth,colour=disease)) +
  geom_smooth(method='gam') +
  geom_abline(intercept = 0, slope = 1, colour='red') +
  scale_x_continuous('GLM predicted probability of readmission',expand=c(0,0),limits=c(0,1)) + scale_y_continuous('Proportion readmitted (GAM Smoothed)')+
  theme(legend.position = 'bottom',text=element_text(family="Cambria"))
ggsave(filename="figures/glmnet_calibration.png", plot=p,dpi=700,
       width=(8.5*0.75) - (0.5*2),
       height=((8.5 - (0.5*2))/3) + 1)

# A summary function for the top covariates.
dump.glmnet.coefs<-function(disease) {
  glmnet.model<-models[[disease]][['glmnet.predict.outcome']]
  match(glmnet.model$lambda.1se,glmnet.model$lambda)
  rownames(as.matrix(coef(glmnet.model)))[c(as.matrix(coef(glmnet.model)!=0))]
}

# dump.glmnet.coefs('pneumonia')
# dump.glmnet.coefs('ami')
# dump.glmnet.coefs('heart_failure')
# And what does the random forest outcome look like?

# Now, make a table for Q*?
# rf.Q.star, calibrated.rf.Q.star, glmnet.Q.star,
# rf.epsilons, calibrated.rf.epsilons, glmnet.epsilons,
# all.rf.Q.by.hosp, all.calibrated.rf.Q.by.hosp, all.glmnet.Q.by.hosp,
# ls(models[[disease]])
dump.base.stats<-function(disease) {
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

  mean.survived.los <- dash.NAs(aggregate(los~hosp,disease.df,mean))[,'los']
  mean.died.los <- dash.NAs(aggregate(los~hosp,died.during.stay,mean))[,'los']
  mean.both.los <- dash.NAs(aggregate(los~hosp,
                                      rbind(disease.df[c('hosp','los')],
                                            died.during.stay[c('hosp','los')]),
                                      mean)) [,'los']

  data.frame(admitted=n+n.died, died=n.died, died.prop=n.died/(n+n.died), live.discharge=n,
             overall.los=mean.both.los,died.los=mean.died.los,survived.los=mean.survived.los,
             readmitted=crude.readmitted.n, prop=crude.props)
}

dump.results.df<-function(disease, model) {
  all.Q.by.hosp <- models[[disease]][[paste0('all.',model,'.Q.by.hosp')]]
  epsilon <- models[[disease]][[paste0(model,'.epsilons')]]
  Q.star.by.epsilon <- models[[disease]][[paste0(model,'.Q.star')]]
  data.frame(Q=colMeans(all.Q.by.hosp),epsilon, Q.star=colMeans(Q.star.by.epsilon))
}

disease.table<-function(disease, models=c('rf','glmnet')) {
  model.tables <- lapply(models, dump.results.df, disease=disease)
  cbind(dump.base.stats(disease), Reduce(cbind, model.tables))
}

disease.results.table<-sapply(prefixes, disease.table, simplify=FALSE)
save(disease.results.table,file='tables/disease.results.table.object')
