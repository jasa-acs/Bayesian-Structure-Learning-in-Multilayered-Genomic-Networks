rm(list=ls())
gc()
library(UpSetR)

setwd("/Users/mjha/MDanderson/ChainGraphModel/Rcode_application/TCGARPPA")
dn.array = c("UCEC","SKCM","OV","READ","COAD","LUAD","LUSC")
ll = length(dn.array)
dn = dn.array[1]
edgeattr = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]

for (dd in 2:ll) {
  dn = dn.array[dd]
  edgeattr = rbind(edgeattr,read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2])
}

NodeA = unlist(lapply(strsplit(as.character(edgeattr[,1]),split="_") ,function(x)x[1]))
NodeB = unlist(lapply(strsplit(as.character(edgeattr[,2]),split="_") ,function(x)x[1]))

# CNA-CNA #
w = intersect(which(NodeA=="CNA"),which(NodeB=="CNA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_CNA_CNA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# Methyl-Methyl #
w = intersect(which(NodeA=="Methyl"),which(NodeB=="Methyl"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_Methyl_Methyl.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# mRNA-mRNA #
w = intersect(which(NodeA=="mRNA"),which(NodeB=="mRNA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_mRNA_mRNA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# Protein-Protein #
w = intersect(which(NodeA=="RPPA"),which(NodeB=="RPPA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_RPPA_RPPA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# CNA-mRNA #
w = intersect(which(NodeA=="CNA"),which(NodeB=="mRNA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_CNA_mRNA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# CNA-Protein #
w = intersect(which(NodeA=="CNA"),which(NodeB=="RPPA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_CNA_RPPA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# Methyl-mRNA #
w = intersect(which(NodeA=="Methyl"),which(NodeB=="mRNA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_Methyl_mRNA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# Methyl-Protein #
w = intersect(which(NodeA=="Methyl"),which(NodeB=="RPPA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_Methyl_RPPA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()

# mRNA-Protein #
w = intersect(which(NodeA=="mRNA"),which(NodeB=="RPPA"))

edgen = unique(paste("(",edgeattr[w,1],",",edgeattr[w,2],")",sep=""))
mat = matrix(0,nrow=length(edgen),ncol=length(dn.array))
for (dd in 1:ll) {
  dn = dn.array[dd]
  tmp = read.table(paste("results/edgeattr_",dn,".txt",sep=""))[,1:2]
  edgetmp = unique(paste("(",tmp[,1],",",tmp[,2],")",sep=""))
  edgetmp1 = intersect(edgetmp,edgen)
  w = match(edgetmp1,edgen)
  mat[w,dd]=1
}
rownames(mat) = paste("sample",1:nrow(mat))
colnames(mat) = dn.array
pdf("UpSet_mRNA_RPPA.pdf",width=6,height=7)
upset(data.frame(mat),nsets=7,text.scale=c(2,2.3,2,2,2,1),shade.alpha=0.5,sets.bar.color="green",point.size=1.5,order.by="freq",mb.ratio = c(0.6, 0.4),show.numbers=F)
dev.off()
