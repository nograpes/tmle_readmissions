# Dump directory
# /usr/bin/R --args ./../data_clean/data_source.R ./data_dump ./disease_subsets/ami ./disease_subsets/heart_failure ./disease_subsets/pneumonia
# setwd('~/repo/thesis/code/tmle')
# data.source.file='../data_clean/data_source.R'
arguments <- commandArgs(trailingOnly=TRUE)
data.source.file<-arguments[1]
dump.dir<-arguments[2]
paths<-arguments[-(1:2)]

# Extremely important to set this before sourcing.. data.source.file uses this
cohort.table <- 'readmissions_top20_one_year'
source(data.source.file)

fixed.vars<-c('admission_type','dob',
              'sex','age','prev_readmissions',
              'insurance_plan','dow','year','month','admit.diag.mdc','csd_ct_uid','hosp')

ff <- as.formula(paste('~',paste(fixed.vars,collapse='+'),'-1'))
# for (var in fixed.vars) df[,var]<-as.character(df[,var])
for (var in c('age','prev_readmissions','transfers')) df[,var]<-as.numeric(df[,var])
for (var in c('dob')) df[,var]<-as.numeric(as.Date(df[,var],'%Y-%m-%d'))

set.to.none<-c('insurance_plan','admit.diag.mdc')
for (var in set.to.none) {
  df[,var]<-as.character(df[,var])
  df[is.na(df[,var]),var]<-'None'
  df[,var]<-as.factor(df[,var])
}

m <- model.frame(ff, df,na.action=NULL)
fixed.vars.mat<-model.matrix(ff,m)

make.matrix<-function(table_name,item_name,threshold=30){
  id.discharge.item <- dbReadTable(con, table_name)
  names(id.discharge.item) <- c('id', 'discharge', item_name)  
  id.discharge.item$id_dat <- with(id.discharge.item, paste(id, discharge))
  counts <- table(id.discharge.item[, item_name])
  above.thresh <- names(which(counts>=threshold))
  items.above.thresh <- id.discharge.item[id.discharge.item[,item_name] %in% above.thresh,]
  removed.ids <- unique(id.discharge.item$id_dat) [!(unique(id.discharge.item$id_dat) %in% unique(items.above.thresh$id_dat))]

  # Some ids have *no* entries.  
  never.ids <- df$id_dat[!(df$id_dat %in% id.discharge.item$id_dat)]

  x <- merge(c(removed.ids, never.ids), NA)
  names(x) <- c('id_dat',item_name)
  items.above.thresh <- rbind(items.above.thresh[,c('id_dat',item_name)], x)

  tab <- table(items.above.thresh[,'id_dat'],items.above.thresh[,item_name])
  item.names <- colnames(tab)
  x <- as.data.frame.matrix(tab)

  n <- do.call(rbind, strsplit(row.names(x), ' '))
  id.discharge.item <- data.frame(n[,1], as.Date(n[,2]),x)
  names(id.discharge.item) <- c('id', 'discharge', item.names)
  row.names(id.discharge.item) <- NULL
  id.discharge.item <- id.discharge.item[order(id.discharge.item$id,id.discharge.item$discharge),]
  as.matrix(id.discharge.item[,-(1:2)])
}
tables <- c('procs_one_year', 'diags_one_year', 'drugs_one_year')
item.names <- c('procedure', 'diagnosis', 'drug')
item.matrices <- mapply(make.matrix,tables,item.names,SIMPLIFY=FALSE) # 378s
names(item.matrices) <- item.names

# Cut the data into several groups.
fixed.var.names<-colnames(fixed.vars.mat)
var.names<-sapply(item.names,function(x)colnames(item.matrices[[x]]))
big.matrix<-do.call(cbind,c(list(fixed.vars.mat),item.matrices))
rm(fixed.vars.mat,item.matrices)
invisible(gc())

# Now, read in all the variable names. The only thing that should be passed is the file names in disease subsets.
# paths<-c('./disease_subsets/ami', './disease_subsets/heart_failure', './disease_subsets/pneumonia')
all.died.during.stay<-died.during.stay

for (path in paths){
  regex<-scan(path,what='character',quiet=TRUE)
  rows<-which(!is.na(df$admit_diag) & grepl(regex,df$admit_diag))
  disease.df<-df[rows,]
  disease.big.matrix<-cbind(big.matrix[rows,])
  # Died during stay
  rows2<-which(!is.na(all.died.during.stay$admit_diag) & grepl(regex,all.died.during.stay$admit_diag))
  died.during.stay<-all.died.during.stay[rows2,]
  save(died.during.stay,var.names,disease.df,disease.big.matrix,file=paste0(dump.dir,'/disease_',basename(path),'.object'))
}
