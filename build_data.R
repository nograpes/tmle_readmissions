# TODO:
# Fix the tables so they don't have to be joined before pulling them in.
# Give the procedure codes English names.
# Give the icd9 codes names.

# Set mirror and ncpus:
options(mc.cores=12)

# Make a simple model of readmissions:
setwd('~/repo/thesis/code/tmle')
# print('Data source.')
source('~/repo/thesis/code/notebook/data_source.R')

# Dump directory
dump.dir<-'data_dump'

# Sort both
df<-df[order(df$id,df$discharge),]
df$csd_ct_uid<-paste(df$csd,df$cma,df$ct,sep='_')

print(paste('Fixed vars: ',Sys.time()))
# This code would work but somehow it causes a segmentation fault.
fixed.vars<-c('admission_type','transfers','to_hosp_type','los','dob','sex','age','prev_readmissions','insurance_plan','dow','year','month','admit.diag.mdc','csd_ct_uid')
ff <- as.formula(paste('~',paste(fixed.vars,collapse='+'),'-1'))
df$to_hosp_type[is.na(df$to_hosp_type)] <- 'None'
df$insurance_plan[is.na(df$insurance_plan)] <- 'None'
df$admit.diag.mdc<-as.character(df$admit.diag.mdc)
df$admit.diag.mdc[is.na(df$admit.diag.mdc)] <- 'None'

df$to_hosp_type<-as.factor(df$to_hosp_type)
df$insurance_plan<-as.factor(df$insurance_plan)
df$admit.diag.mdc<-as.factor(df$admit.diag.mdc)
df$admit_diag<-as.factor(df$admit_diag)

contrasts(df$to_hosp_type) <- contr.treatment(levels(df$to_hosp_type),base=match('None',levels(df$to_hosp_type)))
contrasts(df$insurance_plan) <- contr.treatment(levels(df$insurance_plan),base=match('None',levels(df$insurance_plan)))
contrasts(df$admit.diag.mdc) <- contr.treatment(levels(df$admit.diag.mdc),base=match('None',levels(df$admit.diag.mdc)))

df$id_dat<-paste(df$id,df$discharge)

m <- model.frame(ff, df,na.action=NULL)
fixed.vars.mat<-model.matrix(ff,m)
fixed.vars.mat<-fixed.vars.mat[,colnames(fixed.vars.mat)!='(Intercept)']
# Save and dump.
save(fixed.vars.mat,file=paste(dump.dir,'fixed.vars.mat.object',sep='/'))

# Get the rest from the chandan directory.
