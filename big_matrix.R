# Unfortunately, cbind doesn't support long vectors yet.
# So strip the matrices of a dimension, use c, and then assign the dimension.
dump.dir<-'data_dump'
setwd('~/repo/thesis/code/chandan')
source('~/repo/thesis/code/notebook/data_source.R')
files=c('drugs.ahfs.object','drugs.object','diagnoses.object','procedures.object','just.hosps.object','admit.diags.object')
for (file in files) load(paste(dump.dir,file,sep='/'))
setwd('~/repo/thesis/code/tmle')
load(paste(dump.dir,'fixed.vars.mat.object',sep='/'))

dim(fixed.vars.mat)<-NULL
dim(drugs)<-NULL
dim(procedures)<-NULL
dim(just.hosps)<-NULL
dim(diagnoses)<-NULL
dim(admit.diags)<-NULL

nrow.df<-nrow(df)

big.matrix<-c(fixed.vars.mat,
              drugs,
              procedures,
              just.hosps,
              diagnoses,
              admit.diags)

l<-length(big.matrix)
dim(big.matrix)<-c(nrow.df,l/nrow.df)
rm(fixed.vars.mat,
   drugs,
   procedures,
   just.hosps,
   diagnoses,
   admit.diags)
gc()
save(big.matrix,df,file=paste(dump.dir,'big.matrix.object',sep='/'))