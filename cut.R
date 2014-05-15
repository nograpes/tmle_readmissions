# Cut the data into several groups.
dump.dir<-'data_dump'
setwd('~/repo/thesis/code/chandan')
load('df.object')
files=c('drugs.ahfs.object','drugs.object','diagnoses.object','procedures.object','just.hosps.object','admit.diags.object')
for (file in files) load(paste(dump.dir,file,sep='/'))
# load(paste(dump.dir,'admit.diags.object',sep='/'))

setwd('~/repo/thesis/code/tmle')
load(paste(dump.dir,'fixed.vars.mat.object',sep='/'))

fixed.var.names<-colnames(fixed.vars.mat)
drug.names<-colnames(drugs)
procedure.names<-colnames(procedures)
hosp.names<-colnames(just.hosps)
diagnosis.names<-colnames(diagnoses)
admit.diag.names<-colnames(admit.diags)
save(fixed.var.names,drug.names,procedure.names,hosp.names,diagnosis.names,admit.diag.names,file=paste(dump.dir,'names.object',sep='/'))


big.matrix<-c(fixed.vars.mat,
              drugs,
              procedures,
              just.hosps,
              diagnoses,
              admit.diags)

l<-length(big.matrix)
dim(big.matrix)<-c(nrow(df),l/nrow(df))

# Very important that you don't use lapply(c(fixed.vars,drugs),colnames)... it will copy the entire data set just for the function. Again, stupidly.
column.names<-unlist(lapply(c('fixed.vars.mat','drugs','procedures','just.hosps','diagnoses','admit.diags'),function(x) colnames(get(x))))
# Very important that you don't do colnames(big.matrix).. that will copy the whole matrix, stupidly.
attr(big.matrix,'dimnames')[[2]]<-column.names # The "two" isn't a magic number, it is the (unnamed) dimension, 1 would be rows.

rm(fixed.vars.mat,
   drugs,
   procedures,
   just.hosps,
   diagnoses,
   admit.diags)
invisible(gc())

# Now, read in all the variable names. The only thing that should be passed is the file names in disease subsets.
paths<-commandArgs(trailingOnly=TRUE)

for (path in paths){  
  regex<-scan(path,what='character',quiet=TRUE)
  rows<-which(!is.na(df$admit_diag) & grepl(regex,df$admit_diag))
  disease.df<-df[rows,]
  disease.big.matrix<-cbind(big.matrix[rows,])
  save(disease.df,disease.big.matrix,file=paste0(dump.dir,'/disease_',basename(path),'.object'))
}
