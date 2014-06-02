# This is some preliminary exploration of the problem of no convergence for one of the hospitals in the pneumonia subset.
debugonce(glm.fit)

library(parallel)
options(mc.cores=12)
Y=as.factor(disease.df$day_30_readmit)

rm(list=setdiff(ls(), c('Y','modified.rf.prob.of.outcome','iptw')))
gc()

coefs.by.iteration <- mclapply(90:100, function(x) suppressWarnings(coef(glm(Y ~ .-1, offset=qlogis(modified.rf.prob.of.outcome), data=data.frame(Y=Y,iptw),family=binomial(link=logit),maxit=x))))




r<-do.call(rbind,coefs.by.iteration)
fixed.coefs<-r[11,]

flip<-r[10,"Hôpital.Charles.Lemoyne"]
flop<-r[11,"Hôpital.Charles.Lemoyne"]

eval.points<-seq(min(flip,flop) - abs(flip-flop),max(flip,flop) + abs(flip-flop),by=abs(flip-flop)/100)

# It's Charles Lemoyne that is all flippy floppy.
likelihood<-function(beta) {
  betas<-fixed.coefs
  betas["Hôpital.Charles.Lemoyne"]<-beta
  preds<-plogis(((iptw %*% betas) + qlogis(modified.rf.prob.of.outcome)))
  sum((Y * log(preds)) + ((1-Y)*log(1-preds)))
}

optim(0, function(beta) abs(likelihood(beta)), method='BFGS')$par



plot(eval.points,sapply(eval.points, likelihood), col='red', type='l')
points(flip, likelihood(flip))
points(flop, likelihood(flop))

names(fixed.coefs)

# Calculate *just* for the current coefficient.
betas<-rf.epsilon.model$coef

class(iptw)

# TEst ground for new fitting method.
likelihood<-function(x, y, offset, betas) {
  preds<-plogis(((x %*% betas) + offset))
  sum((y * log(preds)) + ((1-y)*log(1-preds)))
}

# Why even bother trying to push it through GLM?
# It is a little slow but who cares?
# Plus the names don't get mucked up because of rowname character restrictions.
glm.BFGS<-function(x,y,offset=rep(plogis(0),length(Y))) {
  likelihood<-function(betas) {
    preds<-plogis(((x %*% betas) + qlogis(offset)))
    sum((y * log(preds)) + ((1-y)*log(1-preds)))
  }
  setNames(optim(rep(0,ncol(x)), function(betas) abs(likelihood(betas)), method='BFGS')$par,
           colnames(x))
}

glm.BFGS(x=iptw, y=Y, offset=modified.rf.prob.of.outcome)

qlogis(plogis(0))

qlogis(0)

glm.BFGS(x=iptw, y=Y)


likelihood(r[11,])
likelihood(o)

o-r[11,]



optim(rep(0,ncol(x)), 
      function(betas) 
        abs(likelihood(x,y,offset,betas)), method='BFGS')$par


optim.bfgs<-function(x,y,offset, weights,start,etastart,mustart,family,link,linkfun,linkinv,control,intercept) {
  
}

glm(Y ~ .-1, 
    offset=qlogis(modified.rf.prob.of.outcome), 
    data=data.frame(Y=Y,iptw),
    family=binomial(link=logit),
    method=optim.bfgs
)




