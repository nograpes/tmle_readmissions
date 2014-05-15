# Dump directory
# /usr/bin/R --args ./../data_clean/data_source.R ./data_dump/fixed.vars.mat.object
# setwd('~/repo/thesis/code/tmle')
# data.source.file='../data_clean/data_source.R'
data.source.file<-commandArgs(trailingOnly=TRUE)[1]
dump.dir<-commandArgs(trailingOnly=TRUE)[2]
paths<-commandArgs(trailingOnly=TRUE)[-(1:2)]
source(data.source.file)

# This code would work but somehow it causes a segmentation fault.
fixed.vars<-c('admission_type','transfers','to_hosp_type','los','dob','sex','age','prev_readmissions',
              'insurance_plan','dow','year','month','admit.diag.mdc','csd_ct_uid','hosp')

ff <- as.formula(paste('~',paste(fixed.vars,collapse='+'),'-1'))
for (var in fixed.vars) df[,var]<-as.character(df[,var])
for (var in c('age','prev_readmissions','los','transfers')) df[,var]<-as.numeric(df[,var])
for (var in c('dob')) df[,var]<-as.numeric(as.Date(df[,var],'%Y-%m-%d'))

df$to_hosp_type[is.na(df$to_hosp_type)] <- 'None'
df$insurance_plan[is.na(df$insurance_plan)] <- 'None'
df$admit.diag.mdc[is.na(df$admit.diag.mdc)] <- 'None'

m <- model.frame(ff, df,na.action=NULL)
fixed.vars.mat<-model.matrix(ff,m)

make.matrix<-function(table_name,item_name,threshold=30){
  id.discharge.item<-dbReadTable(con,table_name)
  names(id.discharge.item)<-c('id','discharge',item_name)  
  id.discharge.item$id_dat<-with(id.discharge.item,paste(id,discharge))
  counts<-table(id.discharge.item[,item_name])
  above.thresh<-names(which(counts>=threshold))
  items.above.thresh<-id.discharge.item[id.discharge.item[,item_name] %in% above.thresh,]
  removed.ids<-unique(id.discharge.item$id_dat) [!(unique(id.discharge.item$id_dat) %in% unique(items.above.thresh$id_dat))]
  x<-merge(removed.ids,NA)
  names(x)<-c('id_dat',item_name)
  items.above.thresh<-rbind(items.above.thresh[,c('id_dat',item_name)],x)
  tab<-table(items.above.thresh[,'id_dat'],items.above.thresh[,item_name])
  item.names<-colnames(tab)
  x<-as.data.frame.matrix(tab)
  n<-do.call(rbind,strsplit(row.names(x),' '))
  id.discharge.item<-data.frame(n[,1],as.Date(n[,2]),x)
  names(id.discharge.item)<-c('id','discharge',item.names)
  row.names(id.discharge.item)<-NULL
  id.discharge.item<-id.discharge.item[order(id.discharge.item$id,id.discharge.item$discharge),]
  as.matrix(id.discharge.item[,-(1:2)])
}

tables<-c('chandan_procedures','chandan_diagnoses','chandan_drugs')
item.names<-c('procedure','diagnosis','drug')
item.matrices<-mapply(make.matrix,tables,item.names,SIMPLIFY=FALSE) # 378s
names(item.matrices)<-item.names


# Cut the data into several groups.
fixed.var.names<-colnames(fixed.vars.mat)
var.names<-sapply(item.names,function(x)colnames(item.matrices[[x]]))

# save(fixed.var.names,drug.names,procedure.names,hosp.names,diagnosis.names,admit.diag.names,file=paste(dump.dir,'names.object',sep='/'))

big.matrix<-c(fixed.vars.mat,
              item.matrices[['drug']],
              item.matrices[['procedure']],
              item.matrices[['diagnosis']])

rm(fixed.vars.mat,item.matrices)
invisible(gc())

l<-length(big.matrix)
dim(big.matrix)<-c(nrow(df),l/nrow(df))

# Very important that you don't use lapply(c(fixed.vars,drugs),colnames)... it will copy the entire data set just for the function. Again, stupidly.
column.names<-unlist(list(fixed.var.names,var.names))
# Very important that you don't do colnames(big.matrix).. that will copy the whole matrix, stupidly.
attr(big.matrix,'dimnames')[[2]]<-column.names # The "two" isn't a magic number, it is the (unnamed) dimension, 1 would be rows.


# Now, read in all the variable names. The only thing that should be passed is the file names in disease subsets.

for (path in paths){  
  regex<-scan(path,what='character',quiet=TRUE)
  rows<-which(!is.na(df$admit_diag) & grepl(regex,df$admit_diag))
  disease.df<-df[rows,]
  disease.big.matrix<-cbind(big.matrix[rows,])
  save(var.names,disease.df,disease.big.matrix,file=paste0(dump.dir,'/disease_',basename(path),'.object'))
}
