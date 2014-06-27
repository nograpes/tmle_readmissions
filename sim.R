# FUNCTIONS... straight out of build_rf_Q_star_model.R should be moved to improve generalizability.
# This will calculate vote proportions for any rf model
# For a given tree, get predicted class for data
bump.zeroes <- function(x){
  min.baseline<-min(x[x!=0])
  ifelse(x==0, min.baseline, x)
}

predict.by.tree<-function(tree.num, forest, data, cachepath=matrix.cache, n){
  xtype <- as.integer(.Call("CGetType", data@address, PACKAGE = "bigmemory"))
  tree<-forest[[tree.num]]
  treepredict.result <- .Call('treepredictC',
                              data@address, xtype, n, forest, tree, 
                              PACKAGE = "bigrf")
  sapply(seq_along(levels(forest@y)),function(c){
    w <- treepredict.result$testpredclass == c
    ifelse(w,tree@nodewt[treepredict.result$testprednode],0)
  })  
}

rf.probs <- function(rf.model, data, oob.only=TRUE){
  data.pointer <- bigrf:::makex(data, "xtest", cachepath=matrix.cache)
  # Use foreach instead of 
  prediction.by.tree<-foreach(tree=seq_along(rf.model)) %dopar% 
    predict.by.tree(tree.num=tree, forest=rf.model, data=data.pointer, n=nrow(data))

  # Now set all in-bag to NA
  if(oob.only){
    # Assume provided data is the first n of training sample.
    in.bag.mat<-sapply(rf.model, function(tree) tree@insamp!=0)[1:nrow(data),]
    for(col in 1:ncol(in.bag.mat)){
      in.bag <- in.bag.mat[,col]
      prediction.by.tree[[col]][in.bag,] <- 0
    }
    # prediction.by.tree[in.bag.mat]<-NA
  }
  votes <- Reduce("+", prediction.by.tree)
  classes <- names(rf.model@ytable)
  colnames(votes) <- classes
  votes / rowSums(votes)
}

# This returns a disease matrix with the exposure fixed to a certain hospital.
fixed.hosp.data <- function(hosp, data=disease.big.matrix){
  set.to.hosp <- data # Fix A to the hospital.
  set.to.hosp[,grep('^hosp',colnames(data))] <- 0
  var.name <- paste0('hosp',hosp)
  if(var.name %in% colnames(set.to.hosp)) set.to.hosp[,var.name] <- 1
  set.to.hosp
}

rf.prob.by.hosp <- function(rf.model, hosp, oob.only=TRUE, data=disease.big.matrix)
  rf.probs(rf.model=rf.model, data=fixed.hosp.data(hosp,data), oob.only=oob.only)
####### END FUNCTIONS  
  

n<-1e4 
m<-100 # variables
num.hosps <- 20
set.seed(1)
hosp <- sample(1:num.hosps,n, replace=TRUE)

hosp.mat <- sapply(as.character(sort(unique(hosp))), 
                   function(x) as.numeric(hosp==as.numeric(x))) # Don't use model.matrix. You'll regret it.
# hosp.mat <- (model.matrix(~hosp-1,hosp))
colnames(hosp.mat) <- paste0('hosp',colnames(hosp.mat))
hosp.effects <- log(rnorm(num.hosps,mean=1.2,sd=0.2))

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
suppressPackageStartupMessages(library(bigrf))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores=12) # Register a parallel backend

file.backing <- 'my.big.fat.greek.matrix'
desc.file <- paste(file.backing, 'desc', sep='.')
big.x.mat <- big.matrix(full.mat[,-1], ncol=m+num.hosps, nrow=n, type='integer',
                        backingfile=file.backing, descriptorfile=desc.file)

big.x.df<- data.frame(lapply(data.frame(full.mat[,-1]),as.factor))

# Build Q.						
Q.model <- bigrfc(y=as.factor(full.mat[,'y']), 
				   x=big.x.df, 
				   ntrees=1200, 
				   yclasswts=1 / prop.table(table(as.factor(full.mat[,'y']))))

# Build G.
G.model <- bigrfc(y=as.factor(hosp), 
				   x=big.x.df[!grepl('^hosp',colnames(big.x.df))], 
				   ntrees=1200, 
				   yclasswts=1 / prop.table(table(as.factor(hosp))))


				   
				   
# G model - calibrated RF 
votes <- G.model@oobvotes
prob <- votes/rowSums(votes)
rf.prob.of.exposure <- prob[cbind(1:length(hosp),
                            match(hosp,gsub('hosp','',colnames(prob))))]

# Q model - calibrated RF - as observed (A=a)
votes <- Q.model@oobvotes
prob <- votes/rowSums(votes)
vote.prop <- prob[,"1"]

# Very important to scale the vote proportions.
platt.scaler <- 
  glm(y ~ vote.prop, family=binomial(link='logit'))
rf.prob.of.outcome <-  predict(platt.scaler, type='response')

# Q model - calibrated RF - manipulating exposure to each of twenty levels.
# (A=1, A=2..,A=20)


matrix.cache='~'
unscaled.all.rf.Q.by.hosp<-sapply(sort(unique(hosp)),
                                  function(x,rf.model,data) 
                                    rf.prob.by.hosp(rf.model,x, data=data)[,"1"],
                                  rf.model=Q.model,
								  data=big.x.df)

all.rf.Q.by.hosp <- apply(unscaled.all.rf.Q.by.hosp,2, 
                           function(x) predict(platt.scaler,
                                               newdata=data.frame(prob=x),
                                               type='response'))


											   
# Create a matrix of indicator variables.
exposure.mat <- hosp.mat

# Exposure
modified.rf.prob.of.exposure <- bump.zeroes(rf.prob.of.exposure)
# Outcome
modified.rf.prob.of.outcome <- bump.zeroes(rf.prob.of.outcome)

# I want to find one epsilon for each hospital.
# Each column becomes (A==a)/g
iptw <- exposure.mat  / modified.rf.prob.of.exposure 


glm.BFGS<-function(x,y,offset=rep(plogis(0),length(Y))) {
  likelihood<-function(betas) {
    preds<-plogis(((x %*% betas) + qlogis(offset)))
    sum((y * log(preds)) + ((1-y)*log(1-preds)))
  }
  setNames(optim(rep(0,ncol(x)), 
                 function(betas) abs(likelihood(betas)), method='BFGS')$par,
           colnames(x))
}

epsilons<-function(offset, iptw, y=as.numeric(disease.df$day_30_readmit)) {
  glm.BFGS(x=iptw,
           y=y,
           offset=offset)
}

# I have chosen to use only the calibrated IPTW, because it is more theoretically sound. (I really don't care about accuracy here.)
rf.epsilons <- epsilons(modified.rf.prob.of.outcome, iptw, y=y)

Q.star<-function(Q, iptw, epsilons)
				plogis(qlogis(Q) + ((1/iptw) %*% t(epsilons)))

Q.star(all.rf.Q.by.hosp, iptw, rf.epsilons)
				
rf.Q.star <- Q.star(all.rf.Q.by.hosp, 
                    iptw, 
					rf.epsilons)


							  
colMeans(rf.Q.star) [1] / colMeans(rf.Q.star) 
# Truth

hosp.effects[2] 
hosp.effects[1] / hosp.effects