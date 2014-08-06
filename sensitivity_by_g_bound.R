suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(ggplot2))

prefixes <- list.files('disease_subsets')
pretty.names<-c('AMI','Heart failure','Pneumonia')
names(pretty.names)<-prefixes

get.data<-function(prefix) {
  data.files<-'rf_Q_star_model_'
  files<-paste0('data_dump/', data.files, prefix, '.object')
  e<-new.env()
  for(file in files) load(file,envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.

hosps <- colnames(models[[1]]$g.by.rf)
hosp.name <- paste('Hsp',seq_along(hosps))
names(hosp.name) <- hosps

bounds <- models[[1]]$bounds
l = lapply(models, function(x) melt(sapply(x$Q.star.by.bound, colMeans), c('Hospital','idx')))
melted.s <- do.call(rbind, l)
melted.s$disease <- pretty.names[rep(names(l),sapply(l,nrow))]
melted.s$delta <- bounds[melted.s$idx]
colnames(melted.s)[colnames(melted.s)=="value"]<-'Q.star'
melted.s$Hospital <- factor(hosp.name[melted.s$Hospital], levels=hosp.name)


# Load the table file to figure out where we set the bounds.
table.file<-'disease.results.table'
file<-paste0('tables/', table.file, '.object')
load(file)
bounds <- sapply(bound.expression,eval)


breaks=c(0.1,0.075,0.05,0.025,0)
as.character(breaks)

# Plot sensitivity analysis on bounds.
p=ggplot(melted.s, aes(x=delta, y=Q.star, colour=Hospital)) + 
  scale_x_reverse(name='δ', breaks=breaks, labels=as.character(breaks)) +
  geom_line() + facet_wrap(~disease) +
  scale_y_continuous(name='Marginal risk (Q*)') + geom_vline(xintercept = bounds,
                                             colour='red',
                                             lwd=0.5,lty=2) + 
  theme(text = element_text(family="Cambria"))
ggsave(filename='figures/effect_of_gbound.png',p,
       width=(8.5-1), height=11-2,units="in",dpi=600)  

# A histogram of the distribution of the g(A|W)
l = lapply(models,function(x) melt(x$g.by.rf))
melted.g <- do.call(rbind,l)
melted.g$disease <- pretty.names[rep(names(l),sapply(l,nrow))]
colnames(melted.g) <- c('ID','Hospital','g','disease')
melted.g$Hospital <- factor(hosp.name[melted.g$Hospital], levels=hosp.name)

binwidth = c(zoomed=0.0005, full=0.0005)

# Full x: (Display from 0 to 1)
# To difficult to display because of the differences in scale at each end.
# ggplot(melted.g, aes(x=g)) + 
#   geom_bar(binwidth=binwidth['full']) +
#   facet_grid(Hospital~disease) + 
#   # scale_y_continuous(name='Count') +
#   scale_y_log10(name='Count') +
#   scale_x_continuous(name='Probability of exposure g=Pr(A=a|W)') + 
#   geom_vline(xintercept = bounds, colour='red', lwd=0.5,lty=2)
# ggsave(filename='figures/dist_g_full.png',p,
#        width=(8.5-1), height=11-2,units="in",dpi=600)  

# Zoomed 
p=ggplot(melted.g, aes(x=g)) + 
  geom_bar(binwidth=binwidth['zoomed']) +
  facet_grid(Hospital~disease) + 
  scale_y_continuous(limits = c(0,400), name='Count') +
  scale_x_continuous(limits = c(0,0.05), name='Probability of exposure g=Pr(A=a|W)')  + 
  geom_vline(xintercept = bounds, colour='red', lwd=0.5,lty=2)
ggsave(filename='figures/dist_g_zoomed.png',p,
       width=(8.5-1), height=11-2,units="in",dpi=600)  


# I should move this into the document for reproducibility:
#prop.table(table(melted.g$g<10^(-2))) ["TRUE"] # 39 %
#prop.table(table(melted.g$g<10^(-2.5))) ["TRUE"] # 4 %
