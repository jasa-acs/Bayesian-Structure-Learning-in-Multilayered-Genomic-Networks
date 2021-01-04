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

burnin.S = 0
inf.S = 30000

chlist = vector("list",q)
chlist[[1]] = 1:addr[1]
chlist[[2]] = (addr[1]+1):addr[2]
chlist[[3]] = (addr[2]+1):addr[3]
chlist[[4]] = (addr[3]+1):addr[4]


palist = vector("list",q)
palist[[3]] = 1:addr[2]
palist[[4]] = 1:addr[3]


lambda = 5
delta = 2

library(snow)
cl = makeCluster(rep('localhost', q), 'SOCK') 
clusterEvalQ(cl, source("../Rpackage/chaingraph.R")) 
clusterExport(cl, list("Y","chlist","palist","burnin.S","inf.S","lambda","delta"))
a = date()
fit = parSapply(cl,1:q,function(t) ch.chaingraph(v.ch=chlist[[t]],v.pa=palist[[t]],Y=Y,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S),simplify=F)
b = date()
stopCluster(cl)
save(Y,lambda,delta,a,b,fit,file=paste("../appl_results/LUSC_fit_pp",pp,".rda",sep=""))

