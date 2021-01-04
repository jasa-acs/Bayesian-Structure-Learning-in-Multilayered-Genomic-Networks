rm(list=ls())
gc()
r = 1

p  = 200
q = 10
pE = 0.03
n = 100

library(mvtnorm)
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



simul.Chaingraph.prev <- function(p,q,pE) {
# ----------------------------------------
# generate ER model for the skeleton of the chain graph
# p : number of vertices
# q : number of chain components
# pE : connection probability of ER model
# ----------------------------------------
	
	stopifnot(pE>0,pE<1)
# Construct skeleton
	A = matrix(0,p,p)
	w = which(lower.tri(A))
	A[w] = rbinom(n=length(w),size=1,prob=pE) 
	A = A + t(A)
    
# Define dependence chain #
	addr = cbind(1:p,1:p%%q+1)
	ch.addr = sapply(1:q,function(x)addr[addr[,2]==x,1],simplify=F)
	
# Replace undirected edges to directed edges
	for (i in 1:(q-1)) {		
		for (j in 1:i) {
			A[as.matrix(expand.grid(ch.addr[[j]],ch.addr[[i+1]]))] =0 ### Column-> row
		}
	}
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
set.seed(7168)
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
chlist = sim$ch.addr

palist = vector("list",q)
for (i in 2:q) {
	for (j in 1:(i-1)) {
      palist[[i]] = c(chlist[[j]],palist[[i]])
    }
}

###### Fitting ######

library("capme")
library("clime")
lams=rev(10^seq(from=-1, to=0.5, length=20))

fitlist = vector("list",length(lams))
for (j in 1:length(lams)) {
	lam = lams[j]
	G = matrix(0,p,p)
	for (i in 1:q){
		Y = dat$Y[,chlist[[i]]]
		X = NULL
		if (i ==1){
			G[chlist[[i]],chlist[[i]]] = clime(x=Y,lambda=lam)$Omegalist[[1]]
		}else{
			X = dat$Y[,palist[[i]]]
			fit=capme(y=Y, x=X, lambda=lam, tau=lam)
            G[as.matrix(expand.grid(chlist[[i]],palist[[i]]))] = t(fit$Gammalist[[1]])
			G[chlist[[i]],chlist[[i]]] = fit$Omegalist[[1]]
		}
	}
	fitlist[[j]] = G
}

save(dat,sim,fitlist,file=paste("../simul_results/capme_fitlist_p",p,"q",q,"pE",pE,"n",n,"r",r,".rda",sep=""))



