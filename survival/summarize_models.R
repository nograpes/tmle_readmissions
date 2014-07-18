suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

# I'll have to load all three for the summary.
# Read this in programatically.
prefixes<-c('ami','heart_failure','pneumonia')
# prefixes<-c('ami','pneumonia')
pretty.names<-c('AMI','Heart failure','Pneumonia')
# pretty.names<-c('Acute myocardial infarction','Pneumonia')
names(pretty.names)<-prefixes

get.data<-function(prefix) {
  data.files<-c('Q_star_survival_')
  files<-paste0('survival/data_dump/', data.files, prefix, '.object')
  e<-new.env()
  for(file in files) load(file,envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.
  
get.melted.df <- function(tte.mat, pretty.name) {
  melted.tte <- melt(data.frame(tte.mat, check.names=FALSE), 
                     id.vars = NULL,
                     variable.name = "hospital",
                     value.name="tte")
  data.frame(melted.tte, disease=unname(pretty.name))
}

melted.df <- do.call(rbind,mapply(get.melted.df, 
                    lapply(models,"[[","tte.mat"), pretty.names,
                    SIMPLIFY=FALSE))
rownames(melted.df) <- NULL

# Add in the median
median.tte <- aggregate(data.frame(median.tte=melted.df$tte),by=melted.df[c('hospital','disease')],median)
melted.df <- merge(melted.df,median.tte)

melted.df$hosp_num <- match(melted.df$hospital,colnames(models[[1]]$epsilon.mat))

melted.df$hospital_name <- paste('Hsp',melted.df$hosp_num,sep=' ')

melted.df <- melted.df[order(as.numeric(melted.df$hosp_num)),]

# Otherwise the hospitals will appear in character order.
melted.df$hospital_name <- factor(melted.df$hospital_name,levels=unique(melted.df$hospital_name))

p<-ggplot(data=melted.df[melted.df$tte <=60,] , aes(x=tte)) + 
   facet_grid(hospital_name ~ disease) +
   scale_fill_continuous(name = 'Median\ndays-to-readmission') +
   scale_x_continuous(breaks=c(0,30,60),
                      name = 'Marginal days-to-readmission') +
   scale_y_continuous(breaks=c(0,0.02),
                     name = 'Density') +
   geom_density(aes(fill=median.tte), na.rm=TRUE, trim=TRUE) +
   theme(legend.position = 'bottom',
         text = element_text(family="Cambria"),
         strip.text.x = element_text(size=8),
         strip.text.y = element_text(size=8))
ggsave(filename='figures/tte_distribution.png',p,
       width=(8.5-1)/2, height=11-2,units="in",dpi=600)

