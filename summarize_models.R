# /usr/bin/R --args ami heart_failure pneumonia
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))

# Read this in programatically.
prefixes<-c('ami','heart_failure','pneumonia')
pretty.names<-c('Acute myocardial infarction','Heart failure','Pneumonia')
names(pretty.names)<-prefixes
setwd('~/repo/thesis/code/tmle')

get.data<-function(prefix) {
  dump.dir<-'data_dump'
  data.files<-c('rf_G_model_', 'rf_Q_model_', 'rf_Q_star_model_', 'disease_')
  files<-paste0('data_dump/', data.files, prefix, '.object')
  e<-new.env()
  for(file in files) load(file,envir=e)
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

# sapply only assigns names if X is a character. Annoying, but maybe for the best.
rf.exposure.accuracy.by.disease<-sapply(models,function(x) with(x,get.accuracy(rf.predict.exposure, disease.df$hosp)))

colnames(rf.exposure.accuracy.by.disease)<-
  pretty.names[names(models)]

accuracy.df<-melt(rf.exposure.accuracy.by.disease,
                  varnames = c('trees','disease'),
                  value.name='accuracy')

p<-ggplot(accuracy.df,
          aes(x=trees,
              y=1-accuracy,
              col=disease)) + 
  geom_line(size=1.5) +
  labs(x='Number of trees',
       y='Error rate (out-of-bag)') +
  scale_colour_discrete(name = 'Admission diagnosis') +
  theme(legend.position = 'right',
        text=element_text(family="Cambria"))
ggsave(filename="figures/error_rate_for_hospital_choice.png", 
       plot=p,width=(8.5 - (0.5*2)),height=11/3,dpi=300)

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








