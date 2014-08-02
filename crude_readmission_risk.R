suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(RPostgreSQL))
arguments<-commandArgs(trailingOnly=TRUE)
if(interactive()) arguments<-c('data_dump/disease_heart_failure.object','data_dump/crude_readmissions_risk.object')
object.file<-arguments[1]
output.file<-arguments[2]
load(object.file)

# Grab the Charlson Comorbidities?
# Setup database driver.
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="IrisQuebec")
cc <- dbReadTable(con,"charlson_comorbidities_one_year")
cc <- merge(cc, disease.df) [colnames(cc)]
cc$id_dat <- paste(cc$id,cc$discharge)
tab <- table(cc$id_dat,cc$cc)
item.names <- colnames(tab)
x <- as.data.frame.matrix(tab)
n <- do.call(rbind, strsplit(row.names(x), ' '))
colnames(n) <- c('id','discharge')
cc.id.discharge <- cbind(n,x)
rownames(cc.id.discharge) <- NULL
cc.id.discharge$discharge <- as.Date(cc.id.discharge$discharge)
disease.df.cc <- merge(cc.id.discharge,disease.df)

x.vars <- c('age','sex','hosp','prev_readmissions',
            grep('^cc_',colnames(disease.df.cc),value=TRUE))
base=14
contrasts(disease.df.cc$hosp) <- contr.treatment(levels(disease.df$hosp), base=base)


logistic.model <- glm(disease.df.cc$day_30_readmit~.,data=disease.df.cc[x.vars])
coef.mat <- cbind(Estimate=coef(logistic.model),suppressMessages(confint(logistic.model)))
logistic.mat<-exp(coef.mat[match(paste0('hosp',levels(disease.df$hosp))[-base],rownames(coef.mat)),])


manip.hosp <- function(hosp,df) {
  df$hosp <- hosp
  df
}

marginal.risk <- function(hosp)
  mean(predict(logistic.model, 
               newdata=manip.hosp(hosp, disease.df.cc[x.vars])))

marginal.risk.by.hosp <- sapply(levels(disease.df.cc$hosp), marginal.risk)

save(logistic.mat, marginal.risk.by.hosp,
     file=output.file, base=base)
