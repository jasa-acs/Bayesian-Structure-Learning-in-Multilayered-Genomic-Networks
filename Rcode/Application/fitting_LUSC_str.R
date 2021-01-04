rm(list=ls())
gc()
pp = 1

load("../data/StringV10_ppi_pathway.Rdata")
padat = ppi$pathwaydat
padat[,2] = gsub("-","",padat[,2])
pname = unique(padat[,1])[pp]
proteins = toupper(padat[padat[,1] ==pname,2])
genes =unlist(strsplit(padat[padat[,1] ==pname,3],split=", "))

load("../data/TCGA_LUSC_CNA_Methylation_mRNA_RPPA.rda")
sum(!colnames(RPPA)%in% toupper(padat[,2]))
RPPA = RPPA[,colnames(RPPA)%in%proteins]
mRNA = mRNA[,colnames(mRNA)%in%genes]
CNA = CNA[,colnames(CNA)%in%genes]
Methylation = Methylation[,colnames(Methylation)%in%genes]

Y = cbind(Methylation,CNA,mRNA,RPPA)
colnames(Y) = c(paste("Methyl_",colnames(Methylation),sep=""),paste("CNA_",colnames(CNA),sep=""),paste("mRNA_",colnames(mRNA),sep="")
,paste("RPPA_",colnames(RPPA),sep=""))


p1 = ncol(Methylation)
p2 = ncol(CNA)
p3 = ncol(mRNA)
p4 = ncol(RPPA)
addr = cumsum(c(p1,p2,p3,p4))

q=4

chlist = vector("list",q)
chlist[[1]] = 1:addr[1]
chlist[[2]] = (addr[1]+1):addr[2]
chlist[[3]] = (addr[2]+1):addr[3]
chlist[[4]] = (addr[3]+1):addr[4]

palist = vector("list",q)
palist[[3]] = 1:addr[2]
palist[[4]] = 1:addr[3]

load(paste("../appl_results/LUSC_fit_pp",pp,".rda",sep=""))
p = addr[q]

B = matrix(0,p,p)
cB = numeric(0)
for (j in 1:q) {
    tmp = apply(fit[[j]]$eta[,,10001:30000],c(1,2),sum)/20000
    B[chlist[[j]],chlist[[j]]] = tmp
    cB = c(cB,tmp[lower.tri(tmp)])
    if (j%in%c(3,4)) {
        tmp = c(apply(fit[[j]]$Gamma[,,10001:30000],c(1,2),sum)/20000)
        B[as.matrix(expand.grid(chlist[[j]],palist[[j]]))] = tmp
        cB = c(cB,tmp)
    }
}
o = order(cB,decreasing=T)
ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB))<0.2))  ## FDR at 0.2
cut = cB[o][ii]

G= matrix(0,p,p)
G[B>=cut] = 1


burnin.S = 0
inf.S = 30000

lambda = 5
delta = 2

library(snow)
cl = makeCluster(rep('localhost', q), 'SOCK') 
clusterEvalQ(cl, source("../Rpackage/chaingraph.R")) 
clusterExport(cl, list("Y","chlist","palist","burnin.S","inf.S","lambda","delta","G"))
fit = parSapply(cl,1:q,function(t) ch.chaingraph.str(v.ch=chlist[[t]],v.pa=palist[[t]],Y=Y,G=G,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S),simplify=F)
stopCluster(cl)

save(Y,lambda,delta,a,b,fit,file=paste("../appl_results/LUSC_fit_str_pp",pp,".rda",sep=""))



