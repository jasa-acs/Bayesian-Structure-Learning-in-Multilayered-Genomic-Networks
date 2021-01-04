rm(list=ls())
gc()
r = 1 # replication number

#### Simulation Settings
p  = 20; q = 6; pE = 0.3; n = 200
#p  = 100; q = 6; pE = 0.03; n = 200
#p  = 100; q = 10; pE = 0.03; n = 200
#p  = 100; q = 10; pE = 0.03; n = 200
#p  = 200; q = 10; pE = 0.03; n = 100

library(mvtnorm)
scaledMat <- function(x){
    newx=x/sqrt(diag(x) %*% t(diag(x)))
    return(newx)
}

load(paste("../simul_results/fitlist_p",p,"q",q,"pE",pE,"n",n,"r",r,".rda",sep=""))
ch.addr = sim$ch.addr
Y = dat$Y

burnin.S = 10000
inf.S = 20000
chlist = sim$ch.addr

palist = vector("list",q)
for (i in 2:q) {
	for (j in 1:(i-1)) {
      palist[[i]] = c(chlist[[j]],palist[[i]])
    }
}

G = sim$A

lambda = 5
delta = 2


fit = sapply(1:q,ch.chaingraph.str(v.ch=chlist[[t]],v.pa=palist[[t]],Y=Y,G=G,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S),simplify=F) ## This can be replaced by parSapply in the snow package for parallel computation

save(dat,sim,lambda,delta,fit,file=paste("../simul_results/strfitlist_p",p,"q",q,"pE",pE,"n",n,"r",r,".rda",sep=""))



