n<-1e4 
m<-100 # variables
num.hosps <- 20
set.seed(1)
hosp <- as.factor(sample(1:num.hosps,n, replace=TRUE))
hosp.mat <- (model.matrix(~hosp-1,hosp))
hosp.effects <- rnorm(num.hosps,mean=1.2,sd=0.2)

# Very weakly correlated variables.
# Baseline probability of readmission. (But somehow not quite right)
set.seed(12)
intercept.effect <- qlogis(0.20) 
var.effects <- c(intercept.effect, rnorm(m)/sqrt(m))
var.mat <- cbind(1, matrix(sample(0:1,n*m,replace=TRUE),nrow=n))

# Build the big mat and truth
info.mat <- cbind(var.mat,hosp.mat)
truth <- c(var.effects, hosp.effects)
p <- plogis(info.mat %*% truth)
y <- runif(n,0,1) <= p


full.mat<-cbind(y,info.mat[,-1]) # Strip the intercept
colnames(full.mat)[1]<-'y'
colnames(full.mat)[2:(m+1)]<-paste0('x',1:m)

# Fit a big fat model.
model <- glm(y~.,data.frame(full.mat), family=binomial(link='logit'))

# Seems okay. Now let's fit a random forest just to see what works.
library(bigrf)
library(doParallel)
registerDoParallel(cores=12) # Register a parallel backend

file.backing <- 'my.big.fat.greek.matrix'
desc.file <- paste(file.backing, 'desc', sep='.')
big.x.mat <- big.matrix(full.mat[,-1], ncol=m+num.hosps, nrow=n, type='integer',
                        backingfile=file.backing, descriptorfile=desc.file)


rf.model <- bigrfc(y=as.factor(full.mat[,'y']), x=big.x.mat, 
             varnlevels=rep(as.integer(2),ncol(big.x.mat)))



varimp(rf.model, x=big.x.mat)
