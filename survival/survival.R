# long form function.
expand.to.long <- function(censor,A,W) {
}

# It's dump.time.varying2... it's just soooo fast.
dump.time.varying2 <- function(big.mat, censor, discharge, tte) {
  
  # starts<-df[x,'discharge.num']+1
  start <- as.numeric(discharge) + 1
  # ends<-df[x,'readmit.num']
  ends <- start + tte - 1
  
  # Convert to numeric for a faster seq function.
  list.of.sequences<-mcmapply(seq.int, starts, ends, by=1, SIMPLIFY=FALSE)
  time.varying.date<-as.Date(do.call(c,list.of.sequences),origin='1970-01-01')
  dow.by.day<-strftime(time.varying.date,'%A')
  fast.seq<-mclapply(ends-starts+1,seq.int)
  stop<-do.call(c,fast.seq) 
  start<-stop-1
  stacked.times<-data.frame(start,stop,dow.by.day)
  
  # Event calculator
  event<-!df$censor[x]
  event.sequence<-unlist(mclapply(seq_along(event),function(i) if(event[i]) c(rep.int(0,length(list.of.sequences[[i]])-1),1) else rep.int(0,length(list.of.sequences[[i]])) ))
  
  const<-data.table(df[,c('id','admit','discharge','readmit','dow','dow.admit','dob','sex','age','num_drugs','admit.diag.mdc','prev_readmissions','holiday.admit','holiday.discharge','specific.holiday.admit','specific.holiday.discharge','post.holiday.admit','post.holiday.discharge','dvt','pe','dvt_pe','dvt_history','pe_history','dvt_pe_history')])
  lengths<-sapply(list.of.sequences,length)
  row.selector<-unlist(mcmapply(rep.int,x,lengths))
  data.frame(event=event.sequence==1,const[row.selector,],stacked.times)
}

132,210,887

36,450,462,263

# Make a fast GLM fitter.
# Take advantage of two main things. 
# A sparse matrix multiplier, and compute the effect of time sanely.

library(Matrix)
?Matrix
m <- Matrix(disease.big.matrix,sparse=TRUE)

g = m %*% rnorm(ncol(m))

likelihood <- function(coefs, time.coefs){
  
}


explode.data <- function(data){
  cov <- data[rep(seq.int(nrow(data)),times=data$tte),
              setdiff(colnames(data),c('tte','censor'))]  
  
  cov$time <- do.call(c,lapply(data$tte,seq.int))
  cov$tte <- rep(data$tte,times=data$tte)
  cov$censor <- rep(data$censor,times=data$tte)
  cov$event <- (cov$tte == cov$time) & !cov$censor
  cov
}

log.likelihood <- function(coefs, data, outcome){
  p=plogis(data %*% coefs)
  sum(log(ifelse(outcome==1,p,1-p)))
}

set.seed(1)
n <- 1e3
m <- 10
mat <- matrix(sample(c(1,0),n*m, replace=TRUE),nrow=n)
true.coef <- rnorm(10)/m
# actual <- rbinom(n=n,prob=plogis(mat %*% true.coef), size=1)
actual <- c(runif(n) < plogis(mat %*% true.coef))

summary(glm(actual~., family=binomial, data=data.frame(mat)))

baseline.time <- 10
tte <- ceiling(exp(mat %*% true.coef) * baseline.time)


class(exp.param)
debugonce(predict)
predict(exp.param)

getMethod(predict,signature=class(exp.param))
methods(predict)
?survival:::predict.survreg

exp.param <- survreg(Surv(tte)~.,data=data.frame(tte,mat), dist='exponential')

a = (c(1,mat[1,]) %*% coef(exp.param))
lambda <- exp(-a)

plot(x= seq(1,100,by=0.1), y=exp(-lambda*(seq(1,100,by=0.1))), type='l', lwd=3, col='blue')


points(x=1/lambda,y=exp(-lambda*(1/lambda)),col='red')


exp(-1)

1/lambda


predict(exp.param)[1:6]

exp(-a*3) * a



summary(exp.param)

exp(c(1,mat[1,]) %*% coef(exp.param))

exp(- (c(1,mat[1,]) %*% coef(exp.param)))


cumsum(exp(a * (1:10)))


(1 - exp(-a)[1]) ^(1:20)

exp(a)

1/exp(a)

head(tte)

summary(exp.param)
# Scale was fixed at 1




data.frame(data$tte,predict(exp.param))

data <- data.frame(tte, censor=FALSE, mat)
exploded.data <- explode.data(data)
exploded.data$time <- as.factor(exploded.data$time)

head(exploded.data, 50)

log.reg<-glm(event~.,
             data=exploded.data[setdiff(colnames(exploded.data),
                                        c('tte','censor'))],
             family=binomial,
             maxit=500)

table(data$tte)


plogis(3.047e+03)

summary(log.reg)

predict(log.reg,type='response')[7:9]



log.reg

coef(log.reg)

survreg(Surv(tte)~.,)

?coxph

true.coef



head(mat)



true.coef

exp(1)*20


exp(0)


exp(sum(true.coef))*10

ifelse(outcome == 1)

abs(log.likelihood(rep(0,m),mat,actual))
abs(log.likelihood(true.coef,mat,actual))


o = optim(par=rep(0, m), 
          function(coefs, data, outcome) 
            abs(log.likelihood(coefs, data, outcome)), 
          data=mat, outcome=actual)


o

log(seq(0,1))

?optim


head(as.vector(g))

help(package='Matrix')


dump.time.varying2(big.mat=disease.big.matrix, 
                   censor=disease.df$censor, 
                   discharge=disease.df$discharge, 
                   tte=disease.df$tte)



head(disease.df$discharge)
head(disease.df$readmit)

head(disease.df)

head(disease.df$tte)
head(disease.df$censor)


head(colnames(disease.big.matrix),50)

nrow(disease.big.matrix)
nrow(disease.df)



# Estimate Q 
estimate.Q <- function(Y,A,W) {
}

# Estimate g_20 (treatment mechanism)
estimate.censoring <- function(A,W){
  
}

# Estimate g_20 (censoring mechanism)
estimate.censoring <- function(censor,A,W){
  
}

