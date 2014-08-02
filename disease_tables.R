arguments<-commandArgs(trailingOnly=TRUE)
if(interactive()) arguments<-'tables/disease.results.table.object'
output.file<-arguments[1]

prefixes<-list.files('disease_subsets')

pretty.names<-c('Acute myocardial infarction','Heart failure','Pneumonia')
names(pretty.names)<-prefixes

get.data<-function(prefix) {
  dump.dir<-'data_dump'
  data.files<-c('rf_Q_star_model_', 'disease_','crude_readmissions_risk_')
  files<-paste0('data_dump/', data.files, prefix, '.object')
  # files<-c(files,paste0('survival/data_dump/Q_star_survival_',prefix,'.object'))
  e<-new.env()
  for(file in files) load(file,envir=e)
  mget(ls(e),envir=e)
}
models<-sapply(prefixes, get.data, simplify=FALSE) # Using sapply here to get the names.

# Reconstruct QAW.
for (i in seq_along(models)){
  model=models[[i]]
  Q <- model$all.rf.Q.by.hosp
  cols <- match(model$disease.df$hosp,colnames(Q))
  models[[i]]$QAW <- Q[cbind(seq(length(cols)), cols)]
}

# Now, make a table for Q*?
dump.base.stats<-function(disease,los=FALSE) {
  disease.df <- models[[disease]][['disease.df']]
  crude.readmitted.n<-table(disease.df$hosp,disease.df$day_30_readmit)[,'TRUE']
  n<-c(table(disease.df$hosp))
  crude.props<-prop.table(table(disease.df$hosp,disease.df$day_30_readmit),margin=1)[,'TRUE']
  # Died during stay.
  died.during.stay <- models[[disease]] [['died.during.stay']]
  n.died <- c(table(died.during.stay$hosp))
  # Length-of-stay
  dash.NAs <- function(x)
    merge(data.frame(hosp=levels(disease.df$hosp)),x,all.x=TRUE)

  # Logistic model
  logistic.mat <- models[[disease]][['logistic.mat']]
  base <-  models[[disease]]$base
  
  add.base <- function (x) 
    rbind(x[1:base-1,],c(NA,NA,NA),x[base:nrow(x),])
  
  odds <- add.base(logistic.mat)
  colnames(odds) <- c('odds.ratio','odds.ratio.ci.low','odds.ratio.ci.high')
  
  data.frame(admitted=n+n.died, 
             died=n.died, 
             died.prop=n.died/(n+n.died), 
             live.discharge=n, 
             readmitted=crude.readmitted.n, 
             prop=crude.props, 
             odds,
             marginal.risk=models[[disease]]$marginal.risk.by.hosp
             )
}

# The traditional way is to bound:
bound <- function (x, bounds=c(0.025,0.975)) 
{
  x <- ifelse(x>max(bounds),max(bounds),x)
  ifelse(x<min(bounds),min(bounds),x)
}

dump.results.df<-function(disease, bound) {
  Q.star.by.bound <- models[[disease]]$Q.star.by.bound
  bounds <- models[[disease]]$bounds
  Q.star <- Q.star.by.bound[[match(bound, bounds)]]
  average.risk <- colMeans(Q.star)
  base <-  models[[disease]]$base

  # Calculate OR for TMLE
  QAW <- models[[disease]]$QAW
  mu0 <- average.risk[base]
  psi <- (average.risk/(1-average.risk)) / (mu0/(1-mu0))
  g.by.rf <- bound(models[[disease]]$g.by.rf, 
                   bounds=c(bound,1-bound))
  Q.by.rf <- models[[disease]]$all.rf.Q.by.hosp
  hosps <- names(average.risk)
  disease.df <- models[[disease]]$disease.df
  
  base.var.calc <- function(hosp) {
    g1W <- g.by.rf[,hosp]
    Q1W <- Q.by.rf[,hosp]
    Y <- as.numeric(disease.df$day_30_readmit)
    A<-disease.df$hosp
    mu = average.risk[hosp]
    (1/(mu*(1-mu))) * ((((A==hosp) / g1W) * (Y-QAW)) + Q1W)
  }

  base.var <- sapply(hosps, base.var.calc)
  compared.var <- base.var - base.var[,base]
  var.psi <- apply(compared.var,2,var) / nrow(disease.df)
  
  odds.func <- function(x) x/(1-x)
  psi <- odds.func(average.risk) / odds.func(mu0)
  psi.ci <- exp(log(psi) + t((qnorm(0.975) * c(-1,1)) %*% t(var.psi)))
  
  tmle.odds <- cbind(psi,psi.ci)
  tmle.odds [base,] <- NA
  colnames(tmle.odds) <- c('or','low.ci','high.ci')

  data.frame(tmle.odds,
             Q.star=average.risk)
}

disease.table<-function(disease, bound) {
  model.tables <- lapply(bound, dump.results.df, disease=disease)
  cbind(dump.base.stats(disease), Reduce(cbind, model.tables))
}

disease.results.table<-sapply(prefixes, 
                              disease.table, 
                              bound=c(10^-1.6,10^-2), 
                              simplify=FALSE)
save(disease.results.table,file=output.file)
