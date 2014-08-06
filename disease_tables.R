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

bound.expression <- list(expression(10^-2),expression(10^-2.5))
disease.results.table <- sapply(prefixes, 
                                disease.table, 
                                bound=sapply(bound.expression,eval), 
                                simplify=FALSE)

# Also want to calculate a linear correlation between the in-hospital deaths and the proportion readmitted.
# I'm not sure how to account for the variance here.
# 
# df=with(disease.results.table[[2]], data.frame(died, 
#                                             admitted, 
#                                             readmitted,
#                                             discharged=live.discharge))
# 
# # Data
# df = structure(list(died = c(141L, 166L, 134L, 122L, 181L, 107L, 386L, 197L, 111L, 149L, 153L, 162L, 102L, 234L, 190L, 94L, 139L, 212L, 99L, 167L), admitted = c(1229L, 2071L, 1243L, 1076L, 1550L, 827L, 2917L, 1456L, 881L, 1410L, 1297L, 1323L, 1231L, 2110L, 1389L, 681L, 1438L, 1984L, 932L, 1048L), readmitted = c(248L, 441L, 285L, 214L, 288L, 128L, 666L, 232L, 157L, 311L, 258L, 192L, 262L, 424L, 203L, 111L, 328L, 438L, 163L, 171L), discharged = c(1088L, 1905L, 1109L, 954L, 1369L, 720L, 2531L, 1259L, 770L, 1261L, 1144L, 1161L, 1129L, 1876L, 1199L, 587L, 1299L, 1772L, 833L, 881L)), .Names = c("died", "admitted", "readmitted", "discharged"), row.names = c(NA, -20L), class = "data.frame")
# 
# df$died.prop <- df$died / df$admitted
# df$readmitted.prop <- df$readmitted / df$discharged
# 
# #    died admitted readmitted discharged  died.prop readmitted.prop
# # 1   141     1229        248       1088 0.11472742       0.2279412
# # 2   166     2071        441       1905 0.08015451       0.2314961
# # 3   134     1243        285       1109 0.10780370       0.2569883
# # 4   122     1076        214        954 0.11338290       0.2243187
# # 5   181     1550        288       1369 0.11677419       0.2103725
# # 6   107      827        128        720 0.12938331       0.1777778
# # 7   386     2917        666       2531 0.13232773       0.2631371
# # 8   197     1456        232       1259 0.13530220       0.1842732
# # 9   111      881        157        770 0.12599319       0.2038961
# # 10  149     1410        311       1261 0.10567376       0.2466297
# # 11  153     1297        258       1144 0.11796453       0.2255245
# # 12  162     1323        192       1161 0.12244898       0.1653747
# # 13  102     1231        262       1129 0.08285946       0.2320638
# # 14  234     2110        424       1876 0.11090047       0.2260128
# # 15  190     1389        203       1199 0.13678906       0.1693078
# # 16   94      681        111        587 0.13803231       0.1890971
# # 17  139     1438        328       1299 0.09666203       0.2525019
# # 18  212     1984        438       1772 0.10685484       0.2471783
# # 19   99      932        163        833 0.10622318       0.1956783
# # 20  167     1048        171        881 0.15935115       0.1940976
# 
# cor(df$died.prop, df$readmitted.prop) # -0.55

death.and.readmissions.cor <- 
  sapply(disease.results.table, function(x) cor(x$died.prop,x$prop))
           
save(disease.results.table, bound.expression, death.and.readmissions.cor, file=output.file)
