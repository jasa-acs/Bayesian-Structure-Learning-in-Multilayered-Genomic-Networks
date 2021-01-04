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
#library(devtools)
#install_github("MinJinHa/BANS")
library(BANS)

scaledMat <- function(x){
    newx=x/sqrt(diag(x) %*% t(diag(x)))
    return(newx)
}

simul.Chaingraph <- function(p,q,pE) {
# ----------------------------------------
# generate ER model for the skeleton of the chain graph
# p : number of vertices
# q : number of chain components
# pE : connection probability of ER model
# ----------------------------------------
	
	stopifnot(pE>0,pE<1,q>1)
	addr = cbind(1:p,1:p%%q+1)
	ch.addr = sapply(1:q,function(x)addr[addr[,2]==x,1],simplify=F)
	ov = do.call(c,ch.addr)
	
### edges in each chain 
	A = matrix(0,p,p)
	for (i in 1:q) {
		qp = length(ch.addr[[i]])
		Atmp = matrix(0,qp,qp)
		w = which(lower.tri(Atmp))
		Atmp[w] = rbinom(n=length(w),size=1,prob=pE) 
		A[ch.addr[[i]],ch.addr[[i]]] = Atmp + t(Atmp)
	}
    
# edges between chains #
	for (i in 1:(q-1)){
		qp = length(ch.addr[[i+1]])
		qpp = length(ch.addr[[i]])
		tmp = matrix(rbinom(n=qp*qpp,size=1,prob=pE/2) ,ncol=qpp,nrow=qp)
		A[ch.addr[[i+1]],ch.addr[[i]]] = tmp
	}
# Replace undirected edges to directed edges
	return(list(A=A,ch.addr=ch.addr))
}


rmvnorm.Chaingraph <- function(A,ch.addr,n) {
# ----------------------------------------
# generate data given a chain graph
# A : adjacency matrix 
# ch.addr : a list for chain components
# n : sample size
# ----------------------------------------
	p = ncol(A)
	ov = do.call(c,ch.addr)
	stopifnot(sum(1:p%in%ov)==p)
	
	und.edges = which(A!=0 & t(A)!=0,arr.ind=T)
	und.edges = und.edges[!duplicated(t(apply(und.edges,1,sort))),]
	dir.edges = which(A!=0 & t(A)==0,arr.ind=T)
	
# Set undirected part
	K = matrix(0,p,p)
	K[und.edges] = sample(c(-1,1),nrow(und.edges),replace=T)*runif(n=nrow(und.edges),min=0.5,max=1.5)
	K = K + t(K)
	diag(K) = colSums(abs(K)) + 0.1
	R = scaledMat(K)
# Set
    B = matrix(0,p,p)
	B[dir.edges] = sample(c(-1,1),nrow(dir.edges),replace=T)*runif(n=nrow(dir.edges),min=0.5,max=1.5)
# Generate data 
	Ip = diag(p)
	Omega = t(Ip-B)%*%K%*%(Ip-B)
	Sigma = solve(Omega)
	Y = rmvnorm(n,mean=rep(0,p),sigma=Sigma)
	return(list(K=K,B=B,Y=Y,R=R))
}

#### Graph generation
set.seed(77100)
sim  = simul.Chaingraph(p,q,pE)
#plot(graph_from_adjacency_matrix(sim$A))
#dev.off()

## draw graph ##
ch.addr = sim$ch.addr
ov = do.call(c,sim$ch.addr)
A = sim$A[ov,ov]
colnames(A) = rownames(A) = ov

####### Data generation #######
set.seed(101 * (r+77177))
dat = rmvnorm.Chaingraph(A=sim$A,ch.addr=sim$ch.addr,n=n)
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

lambda = 5
delta = 2
fit = sapply(1:q,function(t) ch.chaingraph(v.ch=chlist[[t]],v.pa=palist[[t]],Y=Y,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S),simplify=F) ## This can be replaced by parSapply in the snow package for parallel computation


save(dat,sim,lambda,delta,fit,file=paste("../simul_results/fitlist_p",p,"q",q,"pE",pE,"n",n,"r",r,".rda",sep=""))



